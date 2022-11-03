#ifndef POTENTIAL_HEADER
#define POTENTIAL_HEADER

#include <memory>
#include <complex>


#include "potential/params_base.hpp"
#include "potential/params.hpp"

namespace omplib {

///@brief Base class for an instance of a potential between
/// some target (A,Z) and a projectile at a given cms energy
struct Potential {
  using ptr = std::unique_ptr<Potential>;
  using cptr = std::unique_ptr<const Potential>;

  virtual std::complex<double> eval(double r, double erg) const = 0;
  std::complex<double> operator () (double r, double erg) const { return eval(r, erg); }
};

/// @brief Common phenomenological potential form used for central potentials
class WoodsSaxon : public Potential {
private:
  std::complex<double> V;
  double R, a;
public:
  WoodsSaxon(double V, double R, double a)
    : V(V), R(R), a(a) {};
  WoodsSaxon(std::complex<double> V, double R, double a)
    : V(V), R(R), a(a) {};
  
  std::complex<double> eval(double r, double erg) const final {
    return V/(1. + exp((r-R)/a));
  };
};

/// @brief Common phenomenological potential form used for surface peakes potentials.
/// The derivative in r of a Woods-Saxon
class DerivWoodsSaxon : public Potential {
private:
  std::complex<double> V;
  double R, a;
public:
  DerivWoodsSaxon(double V, double R, double a)
    : V(V), R(R), a(a) {};
  DerivWoodsSaxon(std::complex<double> V, double R, double a)
    : V(V), R(R), a(a) {};
  
  std::complex<double> eval(double r, double erg) const final {
    const auto y = exp((r-R))/a;
    return - V / a * ( y / ( (1. + y) * (1+y) ));
  };
};

/// @brief Common phenomenological potential form used for spin orbit potentials
/// the derivative in r of a Woods-Saxon, times 1/r
class Thomas : Potential {
private:
  std::complex<double> V;
  double R, a;
public:
  Thomas(double V, double R, double a)
    : V(V), R(R), a(a) {};
  Thomas(std::complex<double> V, double R, double a)
    : V(V), R(R), a(a) {};
  
  std::complex<double> eval(double r, double erg) const final {
    return DerivWoodsSaxon(V,R,a).eval(r,erg)/r;
  };
};

class Gaussian : Potential {
private:
  std::complex<double> V;
  double R, sigma;
public:
  Gaussian(double V, double R, double sigma)
    : V(V), R(R), sigma(sigma) {};
  Gaussian(std::complex<double> V, double R, double sigma)
    : V(V), R(R), sigma(sigma) {};
  
  std::complex<double> eval(double r, double erg) const final {
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
  
  std::complex<double> eval(double r, double erg) const final {
    double v = 0;
    for (int i = 0; i < N; ++i) {
      v += gaussians[i].eval(r);
    }
  };
};

using DoubleGaussian = NGaussian<2>;

/// @brief Common global phenomenological OMP form for a given isotope, A,Z.
template<class Params>
class OMP : Potential {
private:
  /// @brief Determines the parameterization used for calculating 
  /// term depths, radii, and diffusivities
  Params params;
  int A;
  int Z;
public:
  OMP(int A, int Z)
    : A(A), Z(Z), params(Params()) {};
  OMP(int A, int Z, Params params)
    : A(A), Z(Z), params(params) {};
  OMP(int A, int Z, Params&& params)
    : A(A), Z(Z), params(params) {};

  void reset_target(int An, int Zn) {
    A = An;
    Z = An;
  }

  std::complex<double> eval(double r, double erg) const final {
    
    const auto V = WoodsSaxon{ 
      std::complex<double>{params.real_cent_V(Z,A,erg),0},
      params.real_cent_R(Z,A,erg),
      params.real_cent_A(Z,A,erg) };
    const auto Vs = DerivWoodsSaxon{ 
      std::complex<double>{params.real_surf_V(Z,A,erg),0},
      params.real_surf_R(Z,A,erg),
      params.real_surf_A(Z,A,erg) };
    const auto Vso = Thomas{ 
      std::complex<double>{params.real_spin_V(Z,A,erg),0},
      params.real_spin_R(Z,A,erg),
      params.real_spin_A(Z,A,erg) };
    
    const auto W = WoodsSaxon{ 
      std::complex<double>{0,params.cmpl_cent_V(Z,A,erg)},
      params.cmpl_cent_R(Z,A,erg),
      params.cmpl_cent_A(Z,A,erg) };
    const auto Ws = DerivWoodsSaxon{ 
      std::complex<double>{0,params.cmpl_surf_V(Z,A,erg)},
      params.cmpl_surf_R(Z,A,erg),
      params.cmpl_surf_A(Z,A,erg) };
    const auto Wso = Thomas{ 
      std::complex<double>{0,params.cmpl_spin_V(Z,A,erg)},
      params.cmpl_spin_R(Z,A,erg),
      params.cmpl_spin_A(Z,A,erg) };

    return V.eval(r,erg) + Vs.eval(r,erg) + Vso.eval(r,erg)  // Re
         + W.eval(r,erg) + Ws.eval(r,erg) + Wso.eval(r,erg); // Im
  };
};

}

#endif 
