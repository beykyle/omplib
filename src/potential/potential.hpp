#ifndef POTENTIAL_HEADER
#define POTENTIAL_HEADER

#include <memory>
#include <complex>

#include "potential/params_base.hpp"
#include "potential/params.hpp"

namespace omplib {

///@brief Base class for an instance of a local potential between
/// some target (A,Z) and a projectile at a given cms energy
struct Potential {
  using ptr = std::unique_ptr<Potential>;
  using cptr = std::unique_ptr<const Potential>;

  /// @brief r is the distance between the projectile and target
  virtual std::complex<double> eval(double r, double erg) const = 0;
  std::complex<double> operator () (double r, double erg) const { return eval(r, erg); }
};

///@brief Base class for an instance of a non-local potential between
/// some target (A,Z) and a projectile at a given cms energy
struct NonLocalPotential {
  using ptr = std::unique_ptr<NonLocalPotential>;
  using cptr = std::unique_ptr<const NonLocalPotential>;

  /// @brief r is a distance between the projectile and target
  /// @brief rp is another distance between the projectile and target
  virtual std::complex<double> eval(double r, double rp, double erg) const = 0;
  std::complex<double> operator () (double r, double rp, double erg) const { 
    return eval(r,rp, erg); }
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
  
  WoodsSaxon(const WoodsSaxon& rhs) = default;
  WoodsSaxon(WoodsSaxon&& rhs) = default;
  WoodsSaxon() = default;
  
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
  
  DerivWoodsSaxon(const DerivWoodsSaxon& rhs) = default;
  DerivWoodsSaxon(DerivWoodsSaxon&& rhs) = default;
  DerivWoodsSaxon() = default;
  
  std::complex<double> eval(double r, double erg) const final {
    const auto y = exp((r-R))/a;
    return - V / a * ( y / ( (1. + y) * (1+y) ));
  };
};

/// @brief Common phenomenological potential form used for spin orbit potentials
/// the derivative in r of a Woods-Saxon, times 1/r
class Thomas : public Potential {
private:
  std::complex<double> V;
  double R, a;
public:
  Thomas(double V, double R, double a)
    : V(V), R(R), a(a) {};
  Thomas(std::complex<double> V, double R, double a)
    : V(V), R(R), a(a) {};
  
  Thomas(const Thomas& rhs) = default;
  Thomas(Thomas&& rhs) = default;
  Thomas() = default;
  
  std::complex<double> eval(double r, double erg) const final {
    return DerivWoodsSaxon(V,R,a).eval(r,erg)/r;
  };
};

class Gaussian : public Potential {
private:
  std::complex<double> V;
  double R, sigma;
public:
  Gaussian(double V, double R, double sigma)
    : V(V), R(R), sigma(sigma) {};
  Gaussian(std::complex<double> V, double R, double sigma)
    : V(V), R(R), sigma(sigma) {};
  
  Gaussian(const Gaussian& rhs) = default;
  Gaussian(Gaussian&& rhs) = default;
  Gaussian() = default;
  
  std::complex<double> eval(double r, double erg) const final {
    return exp( (r-R)*(r-R)/(sigma*sigma) );
  };
};

template<size_t N>
class NGaussian : public Potential {
private:
  std::array<Gaussian,N> gaussians;

public:
  NGaussian(std::array<Gaussian,N> gaussians):
    gaussians(gaussians) {};
  
  NGaussian(const NGaussian& rhs) = default;
  NGaussian(NGaussian&& rhs) = default;
  
  
  std::complex<double> eval(double r, double erg) const final {
    std::complex<double> v = 0;
    for (int i = 0; i < N; ++i) {
      v += gaussians[i].eval(r, erg);
    }
    return v;
  };
};

using DoubleGaussian = NGaussian<2>;

///@brief Thompson, D. R., LeMere, M., & Tang, Y. C. (1977). 
/// Systematic investigation of scattering problems with the resonating-group method. 
/// Nuclear Physics A, 286(1), 53-66.
class Minnesota : public DoubleGaussian {
public:
  static Minnesota build_1S0() {
    return Minnesota(200, -98.15, 1.487, 0.465);
  }
  static Minnesota build_3S1() {
    return Minnesota(200, -178, 1.487, 0.639);
  }
  
  Minnesota(const Minnesota& rhs) = default;
  Minnesota(Minnesota&& rhs) = default;

  Minnesota(double V01, double V02, double k1, double k2):
    DoubleGaussian({Gaussian(V01, 0., 1/sqrt(k1) ), Gaussian(V02, 0., 1/sqrt(k2) )}) {}
};


/// @brief Common global phenomenological OMP form for a given isotope, A,Z.
template<class Params>
class OMP : public Potential {
private:
  /// @brief Determines the parameterization used for calculating 
  /// term depths, radii, and diffusivities
  Params params;


  /// @brief target
  int A, Z;

public:
  OMP(int A, int Z)
    : A(A), Z(Z), params(Params()) {};
  OMP(int A, int Z, Params params)
    : A(A), Z(Z), params(params) {};
  OMP(int A, int Z, Params&& params)
    : A(A), Z(Z), params(params) {};

  OMP(const OMP<Params>& rhs) = default;
  OMP(OMP<Params>&& rhs) = default;

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

/// @brief An arbitrary local potential in r smeared into the off-diagonal 
/// by a Gaussian factor in (r-rp), from:
/// Perey, F., and B. Buck. 
/// "A non-local potential model for the scattering of neutrons by nuclei." 
/// Nuclear Physics 32 (1962): 353-380.
class PereyBuck : public NonLocalPotential {
private:
  Potential::cptr local_potential;
  Gaussian        non_local_factor;
public:
  PereyBuck(Potential::cptr&& potential, double beta)
    : local_potential(std::move(potential))
    , non_local_factor(1/(pow(constants::pi, 3./2.)) * beta*beta*beta, 0, beta) {}

  std::complex<double> eval(double r, double rp, double erg) const final {
    return local_potential->eval(r, erg) * non_local_factor( (r-rp), erg);
  };
};

}

#endif 
