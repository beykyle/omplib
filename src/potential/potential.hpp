#ifndef POTENTIAL_HEADER
#define POTENTIAL_HEADER

#include "potential/params_base.hpp"
#include "potential/params.hpp"

namespace omplib {

///@brief Base class for an instance of a potential between
/// some target (A,Z) and a projectile at a given center of mass
/// energy
struct Potential {
  virtual double eval(double r) const = 0;
  double operator () (double r) const { return eval(r); }
};

/// @brief Common phenomenological potential form used for central potentials
class WoodsSaxon : public Potential {
private:
  double V, R, a;
public:
  WoodsSaxon(double V, double R, double a)
    : V(V), R(R), a(a) {};
  
  double eval(double r) const final {
    return V/(1. + exp((r-R)/a));
  };
};

/// @brief Common phenomenological potential form used for surface peakes potentials.
/// The derivative in r of a Woods-Saxon
class DerivWoodsSaxon : public Potential {
private:
  double V, R, a;
public:
  DerivWoodsSaxon(double V, double R, double a)
    : V(V), R(R), a(a) {};
  
  double eval(double r) const final {
    const auto y = exp((r-R))/a;
    return - V / a * ( y / ( (1. + y) * (1+y) ));
  };
};

/// @brief Common phenomenological potential form used for spin orbit potentials
/// the derivative in r of a Woods-Saxon, times 1/r
class Thomas : Potential {
private:
  double V, R, a;
public:
  Thomas(double V, double R, double a)
    : V(V), R(R), a(a) {};
  
  double eval(double r) const final {
    return Thomas(V,R,a).eval(r)/r;
  };
};

class Gaussian : Potential {
private:
  double V, R, sigma;
public:
  Gaussian(double V, double R, double sigma)
    : V(V), R(R), sigma(sigma) {};
  
  double eval(double r) const final {
    return exp( (r-R)*(r-R)/(sigma*sigma) );
  };
};

template<size_t N>
class NGaussian : Potential {
private:
  std::array<Gaussian,N> gaussians;

public:
  NGaussian(std::array<Gaussian,N> gaussians):
    gaussians(gaussians) {};
  
  double eval(double r) const final {
    double v = 0;
    for (int i = 0; i < N; ++i) {
      v += gaussians[i].eval(r);
    }
  };
};

using DoubleGaussian = NGaussian<2>;

class KoningDelaroche03 : Potential {};
class ChapelHill89 : Potential {};
class WLH21 : Potential {};

}

#endif 
