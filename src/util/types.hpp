#ifndef TYPES_HEADER
#define TYPES_HEADER

#include <complex>

#include "include/boost/rational.hpp"

namespace omplib {
  using real = double;
  using cmpl = std::complex<real>;
  using frac = boost::rational<int>;

  template<typename T> 
  T cast(frac f) {
    return boost::rational_cast<T>(f);
  }

}

#endif
