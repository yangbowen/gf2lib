#pragma once

#include <cstddef>
#include <type_traits>
#include <concepts>
#include <bit>
#include <limits>

template<::std::unsigned_integral T_Result, ::std::unsigned_integral T_Value, bool is_reversed = false>
struct Clmul_Naive {
	using ResultType = T_Result;
	using ValueType = T_Value;
	static_assert(::std::numeric_limits<ResultType>::radix == 2, "not binary");
	static_assert(::std::numeric_limits<ValueType>::radix == 2, "not binary");
	static constexpr int bitwidth_result = ::std::numeric_limits<ResultType>::digits;
	static constexpr int bitwidth_value = ::std::numeric_limits<ValueType>::digits;
private:
	static constexpr ValueType bitmask_value_lsb = is_reversed ? static_cast<ValueType>(1u) << (bitwidth_value - 1) : 1u;
public:
	static constexpr ResultType CalculateClmul(ValueType l, ValueType r) noexcept {
		ResultType result{};
		ResultType tmp_l = is_reversed ? static_cast<ResultType>(l) << (bitwidth_result - bitwidth_value) : l;
		for (int i = 0; i < bitwidth_value; ++i) {
			result ^= ((r & bitmask_value_lsb) ? tmp_l : 0u);
			if constexpr (is_reversed) {
				tmp_l >>= 1;
				r <<= 1;
			} else {
				tmp_l <<= 1;
				r >>= 1;
			}
		}
		return result;
	}
};

template<::std::unsigned_integral T_Result, ::std::unsigned_integral T_Value, bool is_reversed = false>
struct Clmul_Strided {
	using ResultType = T_Result;
	using ValueType = T_Value;
	using StrideType = ::std::underlying_type_t<::std::byte>;
	static_assert(::std::numeric_limits<ResultType>::radix == 2, "not binary");
	static_assert(::std::numeric_limits<ValueType>::radix == 2, "not binary");
	static_assert(::std::numeric_limits<StrideType>::radix == 2, "not binary");
	static constexpr int bitwidth_result = ::std::numeric_limits<ResultType>::digits;
	static constexpr int bitwidth_value = ::std::numeric_limits<ValueType>::digits;
	static constexpr int bitwidth_stride = ::std::numeric_limits<StrideType>::digits;
	static_assert(!(bitwidth_value % bitwidth_stride), "bit width of stride does not evenly divide bit width of value");
	static constexpr int count_stride = bitwidth_value / bitwidth_stride;
private:
	static constexpr StrideType bitmask_stride_lsb = is_reversed ? static_cast<StrideType>(1u) << (bitwidth_stride - 1) : 1u;
	static constexpr void IterateStride(ResultType& stride_result, const ResultType& tmp_l, const StrideType& stride_r) noexcept {
		stride_result ^= ((stride_r & bitmask_stride_lsb) ? tmp_l : 0u);
	}
	template<int... index_stride>
	static constexpr void IterateRound(ResultType(&arr_result)[count_stride], ResultType& tmp_l, ValueType& r, ::std::integer_sequence<int, index_stride...>) noexcept {
		static constexpr int direction_endian = []() consteval {
			int ret = 0;
			if constexpr (::std::endian::native == ::std::endian::little) {
				ret = 1;
			} else if constexpr (::std::endian::native == ::std::endian::big) {
				ret = -1;
			}
			ret *= is_reversed ? -1 : 1;
			return ret;
			}
		();
		if constexpr (direction_endian == 1) {
			(..., IterateStride(arr_result[index_stride], tmp_l, reinterpret_cast<const StrideType*>(&r)[index_stride]));
		} else if constexpr (direction_endian == -1) {
			(..., IterateStride(arr_result[index_stride], tmp_l, reinterpret_cast<const StrideType*>(&r)[count_stride - 1 - index_stride]));
		} else {
			if constexpr (is_reversed) {
				(..., IterateStride(arr_result[index_stride], tmp_l, static_cast<ValueType>((r << (index_stride * bitwidth_stride)) & ~StrideType{})));
			} else {
				(..., IterateStride(arr_result[index_stride], tmp_l, static_cast<ValueType>((r >> (index_stride * bitwidth_stride)) & ~StrideType{})));
			}
		}
		if constexpr (is_reversed) {
			tmp_l >>= 1;
			r <<= 1;
		} else {
			tmp_l <<= 1;
			r >>= 1;
		}
	}
public:
	static constexpr ResultType CalculateClmul(ValueType l, ValueType r) noexcept {
		ResultType arr_result[count_stride]{};
		ResultType tmp_l = is_reversed ? static_cast<ResultType>(l) << (bitwidth_result - bitwidth_value) : l;
		for (int i = 0; i < bitwidth_stride; ++i) IterateRound(arr_result, tmp_l, r, ::std::make_integer_sequence<int, count_stride>{});
		ResultType result{};
		for (int index_stride = count_stride; index_stride > 0; --index_stride) {
			if constexpr (is_reversed) {
				result >>= bitwidth_stride;
				result ^= arr_result[index_stride - 1];
			} else {
				result <<= bitwidth_stride;
				result ^= arr_result[index_stride - 1];
			}
		}
		return result;
	}
};
