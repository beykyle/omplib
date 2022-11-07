#include "potential/params_base.hpp"

using namespace omplib;

real OMParams<Proj::proton>::asym(int Z, int A) {
  const real n = static_cast<real>(A - Z);
  const real a = static_cast<real>(A);
  const real z = static_cast<real>(Z);
  const real alpha = (n - z)/a; 
  return alpha;
}

template<>
real OMParams<Proj::neutron>::asym(int Z, int A) {
  return -OMParams<Proj::proton>::asym(Z,A);
}
