#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <concepts>
#include <bit>
#include <limits>
#include <iostream>
#include <random>

#include "clmul.h"

int main(int, const char**) {
	::std::default_random_engine randengine(::std::random_device{}());
	const unsigned int count_round = 1000;
	unsigned int i = 0u;
	for (i = 0; i < count_round; ++i) {
		unsigned int l = randengine();
		unsigned int r = randengine();
		::std::uintmax_t result_strided = Clmul_Strided<::std::uintmax_t, unsigned int>::CalculateClmul(l, r);
		::std::uintmax_t result_naive = Clmul_Naive<::std::uintmax_t, unsigned int>::CalculateClmul(l, r);
		if (result_strided != result_naive) {
			::std::cout << "clmul-test[" << i << "/" << count_round << "] l=" << l << " r=" << r;
			::std::cout << " result_strided=" << result_strided << " result_naive=" << result_naive << ::std::endl;
			break;
		}
	}
	if (i < count_round) {
		::std::cerr << "FAILED: clmul-test" << ::std::endl;
	} else {
		::std::cerr << "PASSED: clmul-test" << ::std::endl;
	}
	for (i = 0; i < count_round; ++i) {
		unsigned int l = randengine();
		unsigned int r = randengine();
		::std::uintmax_t result_strided = Clmul_Strided<::std::uintmax_t, unsigned int, true>::CalculateClmul(l, r);
		::std::uintmax_t result_naive = Clmul_Naive<::std::uintmax_t, unsigned int, true>::CalculateClmul(l, r);
		if (result_strided != result_naive) {
			::std::cout << "clmul-test-reversed[" << i << "/" << count_round << "] l=" << l << " r=" << r;
			::std::cout << " result_strided=" << result_strided << " result_naive=" << result_naive << ::std::endl;
			break;
		}
	}
	if (i < count_round) {
		::std::cerr << "FAILED: clmul-test-reversed" << ::std::endl;
	} else {
		::std::cerr << "PASSED: clmul-test-reversed" << ::std::endl;
	}
}
