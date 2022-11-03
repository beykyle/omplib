#ifndef __FACTORIAL_HEADER__
#define __FACTORIAL_HEADER__

#include <cstdint>

namespace omplib {

/// @brief Computes the factorial for a positive integer
const uint64_t factorial(const uint32_t n);

/// @brief Computes the doulbe factorial for a positive integer
///
/// See http://mathworld.wolfram.com/DoubleFactorial.html
const uint64_t doubleFactorial(const uint32_t n);

} // namespace util

#endif
