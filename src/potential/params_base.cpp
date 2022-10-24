#include "potential/params_base.hpp"

using namespace omplib;

double OMParams<Proj::proton>::asym(int Z, int A) {
  const double n = static_cast<double>(A - Z);
  const double a = static_cast<double>(A);
  const double z = static_cast<double>(Z);
  const double alpha = (n - z)/a; 
  return alpha;
}

template<>
double OMParams<Proj::neutron>::asym(int Z, int A) {
  const double n = static_cast<double>(A - Z);
  const double a = static_cast<double>(A);
  const double z = static_cast<double>(Z);
  const double alpha = (n - z)/a; 
  return -alpha;
}
