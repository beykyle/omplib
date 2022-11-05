
#ifndef ASYMPTOTICS_HEADER
#define ASYMPTOTICS_HEADER

#include <cmath>
#include <complex>

#include "util/constants.hpp"

namespace omplib {

///@brief Spherical Bessel function of 1st kind 
struct j {
  int l;
  j(int l) : l(l) {};
  std::complex<double> operator () (double kr) {
    return std::sph_bessel(l,kr);
  }
};

///@brief Spherical Bessel function double
struct eta {
  int l;
  eta(int l) : l(l) {};
  std::complex<double> operator () (double kr) {
    return std::sph_neumann(l,kr);
  }
};

///@brief Spherical Hankel function double
struct h_out {
  int l;
  h_out(int l) : l(l) {};
  std::complex<double> operator () (double kr) {
    using constants::i;
    return j{l}(kr) + i * eta{l}(kr);
  }
};

///@brief Spherical Hankel function double
struct h_in {
  int l;
  h_in(int l) : l(l) {};
  std::complex<double> operator () (double kr) {
    using constants::i;
    return j{l}(kr) - i * eta{l}(kr);
  }
};

}
#endif
