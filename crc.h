#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <concepts>
#include <limits>
#include <span>
#include <array>

template<typename T>
concept CRCDef = requires {
	typename T::ValueType;
		requires ::std::unsigned_integral<typename T::ValueType>;
	{ T::is_msb_first } -> ::std::same_as<const bool&>;
	{ T::modulus } -> ::std::same_as<const typename T::ValueType&>;
	{ T::initial } -> ::std::same_as<const typename T::ValueType&>;
	{ T::final_xor } -> ::std::same_as<const typename T::ValueType&>;
};

template<CRCDef T_CRCDef>
struct CRC {
	using DefType = T_CRCDef;
	using ValueType = typename T_CRCDef::ValueType;
	using ByteUIntType = ::std::underlying_type_t<::std::byte>;
	static_assert(::std::numeric_limits<ValueType>::radix == 2, "not binary");
	static constexpr bool is_msb_first = DefType::is_msb_first;
	static constexpr ValueType modulus = DefType::modulus;
	static constexpr ValueType initial = DefType::initial;
	static constexpr ValueType final_xor = DefType::final_xor;
private:
	static constexpr int bitwidth_value = ::std::numeric_limits<ValueType>::digits;
	static constexpr int bitwidth_byte = ::std::numeric_limits<ByteUIntType>::digits;
	static constexpr ::std::size_t size_table_byte = static_cast<::std::size_t>(::std::numeric_limits<ByteUIntType>::max()) + 1;
public:
	/*
	* To describe in GF(2) arithmetic (using the CRC polynomial as the divisor),
	* this lookup table for Sarwate algorithm corresponds to
	* the results (remainders) of the table indices multiplied by 2**n. That is:
	*
	* table[i] = (i * (2**bitwidth_value)) % polynomial
	*/
	static constexpr ::std::array<ValueType, size_table_byte> GenerateLookupTable() noexcept {
		/*
		* This table generation algorithm is equivalent to calculating for each individual bit
		* by iterative multiplication, then summing up the bits to yield the result:
		*
		* table[0] = 0
		* table[2**0] = (2**bitwidth_value) % polynomial
		* for i from 1 to (bitwidth_byte - 1):
		* 	table[2**i] = (table[2**(i - 1)] * 2) % polynomial
		* 	for j from 0 to (i - 1):
		* 		table[2**i + j] = table[i] + table[j]
		*
		* For LSB-first, both the indices and the values are bit-swapped.
		*/
		::std::array<ValueType, size_table_byte> arr_table{};
		if constexpr (is_msb_first) {
			ValueType crc_singlebit = static_cast<ValueType>(1u) << (bitwidth_value - 1);
			for (::std::size_t i = 1; i < size_table_byte; i <<= 1) {
				crc_singlebit = (crc_singlebit << 1) ^ ((crc_singlebit & (static_cast<ValueType>(1u) << (bitwidth_value - 1))) ? modulus : ValueType{});
				for (::std::size_t j = 0; j < i; ++j) arr_table[i + j] = crc_singlebit ^ arr_table[j];
			}
		} else {
			ValueType crc_singlebit = 1u;
			for (::std::size_t i = size_table_byte / 2; i; i >>= 1) {
				crc_singlebit = (crc_singlebit >> 1) ^ ((crc_singlebit & static_cast<ValueType>(1u)) ? modulus : ValueType{});
				for (::std::size_t j = 0; j < size_table_byte; j += (i << 1)) arr_table[i + j] = crc_singlebit ^ arr_table[j];
			}
			return arr_table;
		}
	}
	inline static constexpr ::std::array<ValueType, size_table_byte> table = GenerateLookupTable();
	/*
	* To describe in GF(2) arithmetic (using the CRC polynomial as the divisor),
	* this lookup table corresponds to powers of 2**bitwidth_byte. That is:
	*
	* table[i] = (2**(i*bitwidth_byte)) % polynomial
	*/
	//static constexpr ::std::array<ValueType, size_table_byte> GenerateLookupTable_Exp2() noexcept {
	//	::std::array<ValueType, size_table_byte> arr_table_exp2{};
	//	if constexpr (is_msb_first) {
	//		arr_table_exp2[0] = table[1u];
	//		for (::std::size_t i = 0; i < size_table_byte; ++i) {
	//			arr_table_exp2[i];
	//		}
	//	} else {
	//		ValueType crc_singlebit = 1u;
	//		for (::std::size_t i = size_table_byte / 2; i; i >>= 1) {
	//			crc_singlebit = (crc_singlebit >> 1) ^ ((crc_singlebit & static_cast<ValueType>(1u)) ? modulus : ValueType{});
	//			for (::std::size_t j = 0; j < size_table_byte; j += (i << 1)) arr_table[i + j] = crc_singlebit ^ arr_table[j];
	//		}
	//		return arr_table;
	//	}
	//}
	// Sarwate algorithm for updating the CRC with appended bytes.
	static constexpr void UpdateCRC(ValueType& crc, ::std::span<const ::std::byte> span_byte) noexcept {
		/*
		* In GF(2) arithmetic:
		*
		* CRC(M1..M)
		* = CRC((2**bitwidth_byte) * M1 + M)
		* = ((2**bitwidth_byte) * CRC(M1) + CRC(M)) % polynomial
		* = ((2**bitwidth_byte) * CRC(M1) + (2**bitwidth_value) * M) % polynomial
		* = ((2**bitwidth_byte) * (
		* 		CRC(M1)
		* 			/ (2**(bitwidth_value-bitwidth_byte))
		* 			* (2**(bitwidth_value-bitwidth_byte))
		* 		+ CRC(M1)
		* 			% (2**(bitwidth_value-bitwidth_byte))
		* 	) + (2**bitwidth_byte) * (
		* 		M
		* 			* (2**(bitwidth_value-bitwidth_byte))
		* 	)) % polynomial
		* = (2**bitwidth_value)
		* 		* (CRC(M1) / (2**(bitwidth_value-bitwidth_byte)) + M)
		* 		% polynomial
		* 	+ (2**bitwidth_byte)
		* 		* (CRC(M1) % (2**(bitwidth_value-bitwidth_byte)))
		* 		% polynomial
		* = table[CRC(M1) / (2**(bitwidth_value-bitwidth_byte)) + M] + (2**bitwidth_byte) * (CRC(M1) % (2**(bitwidth_value-bitwidth_byte)))
		*
		* Again, for LSB-first, both the inputs and the CRC values are bit-swapped.
		* That is, The LSB of each input byte becomes the highest-order term in polynomial representation.
		*/
		if constexpr (is_msb_first) {
			for (::std::byte byte : span_byte) {
				crc = table[static_cast<ByteUIntType>(byte) ^ (crc >> (bitwidth_value - bitwidth_byte))] ^ (crc << bitwidth_byte);
			}
		} else {
			for (::std::byte byte : span_byte) {
				crc = table[static_cast<ByteUIntType>(byte) ^ (crc & (static_cast<ValueType>(1u) << (bitwidth_byte - 1)))] ^ (crc >> bitwidth_byte);
			}
		}
	}
	static constexpr ValueType InitializeCRC() noexcept { return initial; }
	static constexpr void FinalizeCRC(ValueType& crc) { crc ^= final_xor; }
	static constexpr ValueType CalculateCRC(::std::span<const ::std::byte> span_byte) noexcept {
		ValueType crc = InitializeCRC();
		UpdateCRC(crc, ::std::move(span_byte));
		FinalizeCRC(crc);
		return crc;
	}
};

struct CRCDef_CRC32 {
	using ValueType = ::std::uint32_t;
	static constexpr bool is_msb_first = false;
	static constexpr ValueType modulus = 0xEDB88320u;
	static constexpr ValueType initial = 0xFFFFFFFFu;
	static constexpr ValueType final_xor = 0xFFFFFFFFu;
};

extern template struct CRC<CRCDef_CRC32>;
using CRC32 = CRC<CRCDef_CRC32>;
