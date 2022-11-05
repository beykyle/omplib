
#ifndef RMATRIX_HEADER
#define RMATRIX_HEADER

#include <cmath>
#include <complex>

#include "util/constants.hpp"

namespace omplib {
///@brief Spherical Bessel function of 1st kind 
struct j {
  int l;
  j(int l) : l(l) {};
  std::complex<double> operator () (std::complex<double> kr) {
    return std::sph_bessel(l,kr);
  }
};

///@brief Spherical Bessel function of 2nd kind 
struct eta {
  int l;
  eta(int l) : l(l) {};
  std::complex<double> operator () (std::complex<double> kr) {
    return std::sph_neumann(l,kr);
  }
};

///@brief Spherical Hankel function of 1st kind 
struct h_out {
  int l;
  h_out(int l) : l(l) {};
  std::complex<double> operator () (std::complex<double> kr) {
    using constants::i;
    return j{l}(kr) + i * eta{l}(kr);
  }
};

///@brief Spherical Hankel function of 2nd kind 
struct h_in {
  int l;
  h_in(int l) : l(l) {};
  std::complex<double> operator () (std::complex<double> kr) {
    using constants::i;
    return j{l}(kr) - i * eta{l}(kr);
  }
};

///@brief derivative of spherical Bessel function of 1st kind 
struct jp {
  int l;
  jp(int l) : l(l) {};
  std::complex<std::complex<std::complex<double>>> operator () (std::complex<double> kr) {
    //TODO
  }
};

///@brief derivative of spherical Bessel function of 2nd kind 
struct etap {
  int l;
  etap(int l) : l(l) {};
  std::complex<std::complex<std::complex<double>>> operator () (std::complex<double> kr) {
    //TODO
  }
};

///@brief derivative of spherical Hankel function of 1st kind 
struct hp_out {
  int l;
  hp_out(int l) : l(l) {};
  std::complex<std::complex<std::complex<double>>> operator () (std::complex<double> kr) {
    //TODO
  }
};

///@brief derivative of spherical Hankel function of 2nd kind 
struct hp_in {
  int l;
  hp_in(int l) : l(l) {};
  std::complex<std::complex<std::complex<double>>> operator () (std::complex<double> kr) {
    //TODO
  }
};

}
#endif
