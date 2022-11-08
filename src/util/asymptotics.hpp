
#ifndef ASYMPTOTICS_HEADER
#define ASYMPTOTICS_HEADER

#include <cmath>
#include <complex>

#include "util/constants.hpp"
#include "util/types.hpp"

namespace omplib {

///@brief Spherical Bessel function of 1st kind 
struct j {
  int l;
  j(int l) : l(l) {};
  cmpl operator () (real kr) {
    return std::sph_bessel(l,kr);
  }
};

///@brief Spherical Bessel function of 2nd kind
struct eta {
  int l;
  eta(int l) : l(l) {};
  cmpl operator () (real kr) {
    return std::sph_neumann(l,kr);
  }
};

///@brief Spherical outgoing Hankel function
struct h_out {
  int l;
  h_out(int l) : l(l) {};
  cmpl operator () (real kr) {
    using constants::i;
    return j{l}(kr) + i * eta{l}(kr);
  }
};

///@brief Spherical incoming Hankel function
struct h_in {
  int l;
  h_in(int l) : l(l) {};
  cmpl operator () (real kr) {
    using constants::i;
    return j{l}(kr) - i * eta{l}(kr);
  }
};

}
#endif
