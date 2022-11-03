#include "util/factorial.hpp"

using namespace omplib;

// It might be worth implementing a cache of precomputed factorial and double
// factorial values to improve performace in the future. For right now, the
// problems being run do not need to compute large factorials.

const uint64_t omplib::factorial(const uint32_t n) {
  return n <= 1 ? 1 : n * factorial(n - 1);
}

const uint64_t omplib::doubleFactorial(const uint32_t n) {
  return n <= 1 ? 1 : n * doubleFactorial(n - 2);
}
