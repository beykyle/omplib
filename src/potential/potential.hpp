#ifndef POTENTIAL_HEADER
#define POTENTIAL_HEADER

#include <memory>
#include <complex>

#include "util/types.hpp"
#include "solver/channel.hpp"

#include "potential/params_base.hpp"
#include "potential/params.hpp"

namespace omplib {

struct GeneralPotential {
  using ptr = std::unique_ptr<GeneralPotential>;
  using cptr = std::unique_ptr<const GeneralPotential>;
  
  virtual cmpl eval(real r, real, const Channel& ch) const = 0;
  virtual cmpl eval_reduced(real r, real, const Channel& ch) const  = 0;
};


///@brief Base class for an instance of a local potential between
/// some target (A,Z) and a projectile at a given cms ench.Tlaby
struct Potential : public GeneralPotential {
  using ptr = std::unique_ptr<Potential>;
  using cptr = std::unique_ptr<const Potential>;

  /// @brief r is the distance between the projectile and target
  virtual cmpl eval(real r, const Channel& ch) const = 0;
  cmpl eval_reduced(real r, const Channel& ch) const { return eval(r, ch) * r; }
  
  cmpl eval(real r, real rp, const Channel& ch) const final { 
    if (r!=rp) return 0;
    return eval(r, ch); 
  }
  cmpl eval_reduced(real r, real rp, const Channel& ch) const final { 
    if (r!=rp) return 0;
    return eval_reduced(r,ch); 
  }
  
  cmpl operator () (real r, const Channel& ch) const { return eval(r, ch); }
};

///@brief Base class for an instance of a non-local potential between
/// some target (A,Z) and a projectile at a given cms ench.Tlaby
struct NonLocalPotential : public GeneralPotential {
  using ptr = std::unique_ptr<NonLocalPotential>;
  using cptr = std::unique_ptr<const NonLocalPotential>;

  /// @brief r is a distance between the projectile and target
  /// @brief rp is another distance between the projectile and target
  cmpl eval_reduced(real r, real rp, const Channel& ch) const final { 
    return eval(r, rp, ch) * r * rp; 
  }
  cmpl operator () (real r, real rp, const Channel& ch) const { 
    return eval(r, rp, ch); }
};

/// @brief Common phenomenological potential form used for central potentials
class WoodsSaxon : public Potential {
private:
  cmpl V;
  real R, a;
public:
  WoodsSaxon(real V, real R, real a)
    : V(V), R(R), a(a) {};
  WoodsSaxon(cmpl V, real R, real a)
    : V(V), R(R), a(a) {};
  
  WoodsSaxon(const WoodsSaxon& rhs) = default;
  WoodsSaxon(WoodsSaxon&& rhs) = default;
  WoodsSaxon() = default;
  
  cmpl eval(real r, const Channel&) const final {
    return V/(1. + exp((r-R)/a));
  };
};

/// @brief Common phenomenological potential form used for surface peakes potentials.
/// The derivative in r of a Woods-Saxon
class DerivWoodsSaxon : public Potential {
private:
  cmpl V;
  real R, a;
public:
  DerivWoodsSaxon(real V, real R, real a)
    : V(V), R(R), a(a) {};
  DerivWoodsSaxon(cmpl V, real R, real a)
    : V(V), R(R), a(a) {};
  
  DerivWoodsSaxon(const DerivWoodsSaxon& rhs) = default;
  DerivWoodsSaxon(DerivWoodsSaxon&& rhs) = default;
  DerivWoodsSaxon() = default;
  
  cmpl eval(real r, const Channel&) const final {
    const auto y = exp((r-R))/a;
    return - V / a * ( y / ( (1. + y) * (1+y) ));
  };
};

/// @brief Common phenomenological potential form used for spin orbit potentials
/// the derivative in r of a Woods-Saxon, times 1/r
class Thomas : public Potential {
private:
  cmpl V;
  real R, a;
public:
  Thomas(real V, real R, real a)
    : V(V), R(R), a(a) {};
  Thomas(cmpl V, real R, real a)
    : V(V), R(R), a(a) {};
  
  Thomas(const Thomas& rhs) = default;
  Thomas(Thomas&& rhs) = default;
  Thomas() = default;
  
  cmpl eval(real r, const Channel& ch) const final {
    return DerivWoodsSaxon(V,R,a).eval(r,ch)/r;
  };
};

class Gaussian : public Potential {
private:
  cmpl V;
  real R, sigma;
public:
  Gaussian(real V, real R, real sigma)
    : V(V), R(R), sigma(sigma) {};
  Gaussian(cmpl V, real R, real sigma)
    : V(V), R(R), sigma(sigma) {};
  
  Gaussian(const Gaussian& rhs) = default;
  Gaussian(Gaussian&& rhs) = default;
  Gaussian() = default;
  
  cmpl eval(real r, const Channel&) const final {
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
  
  
  cmpl eval(real r, const Channel& ch) const final {
    cmpl v = 0;
    for (unsigned int i = 0; i < N; ++i) {
      v += gaussians[i].eval(r, ch);
    }
    return v;
  };
};

using realGaussian = NGaussian<2>;

///@brief Thompson, D. R., LeMere, M., & Tang, Y. C. (1977). 
/// Systematic investigation of scattering problems with the resonating-group method. 
/// Nuclear Physics A, 286(1), 53-66.
class Minnesota : public realGaussian {
public:
  static Minnesota build_1S0() {
    return Minnesota(200, -98.15, 1.487, 0.465);
  }
  static Minnesota build_3S1() {
    return Minnesota(200, -178, 1.487, 0.639);
  }
  
  Minnesota(const Minnesota& rhs) = default;
  Minnesota(Minnesota&& rhs) = default;

  Minnesota(real V01, real V02, real k1, real k2):
    realGaussian({Gaussian(V01, 0., 1/sqrt(k1) ), Gaussian(V02, 0., 1/sqrt(k2) )}) {}
};


/// @brief Common global phenomenological OMP form for a given isotope, A,Z.
template<class Params>
class OMP : public Potential {
private:
  /// @brief target
  int A, Z;

  
  /// @brief Determines the parameterization used for calculating 
  /// term depths, radii, and diffusivities
  Params params;


public:
  OMP(int A, int Z)
    : A(A), Z(Z), params(Params()) {};
  OMP(int A, int Z,  Params params)
    : A(A), Z(Z),  params(params) {};
  OMP(int A, int Z,  Params&& params)
    : A(A), Z(Z), params(params) {};

  OMP(const OMP<Params>& rhs) = default;
  OMP(OMP<Params>&& rhs)      = default;

  void reset_target(int An, int Zn) {
    A = An;
    Z = Zn;
  }

  cmpl eval(real r, const Channel& ch) const final {
    
    const auto V = WoodsSaxon{ 
      cmpl{params.real_cent_V(Z,A,ch.Tlab),0},
      params.real_cent_r(Z,A,ch.Tlab),
      params.real_cent_a(Z,A,ch.Tlab) };
    const auto Vs = DerivWoodsSaxon{ 
      cmpl{params.real_surf_V(Z,A,ch.Tlab),0},
      params.real_surf_r(Z,A,ch.Tlab),
      params.real_surf_a(Z,A,ch.Tlab) };
    const auto Vso = Thomas{ 
      cmpl{params.real_spin_V(Z,A,ch.Tlab),0},
      params.real_spin_r(Z,A,ch.Tlab),
      params.real_spin_a(Z,A,ch.Tlab) };
    
    const auto W = WoodsSaxon{ 
      cmpl{0,params.cmpl_cent_V(Z,A,ch.Tlab)},
      params.cmpl_cent_r(Z,A,ch.Tlab),
      params.cmpl_cent_a(Z,A,ch.Tlab) };
    const auto Ws = DerivWoodsSaxon{ 
      cmpl{0,params.cmpl_surf_V(Z,A,ch.Tlab)},
      params.cmpl_surf_r(Z,A,ch.Tlab),
      params.cmpl_surf_a(Z,A,ch.Tlab) };
    const auto Wso = Thomas{ 
      cmpl{0,params.cmpl_spin_V(Z,A,ch.Tlab)},
      params.cmpl_spin_r(Z,A,ch.Tlab),
      params.cmpl_spin_a(Z,A,ch.Tlab) };

    return V.eval(r,ch) + Vs.eval(r,ch) + Vso.eval(r,ch) * ch.spin_orbit()  // R
         + W.eval(r,ch) + Ws.eval(r,ch) + Wso.eval(r,ch) * ch.spin_orbit(); // Im
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
  PereyBuck(Potential::cptr&& potential, real beta)
    : local_potential(std::move(potential))
    , non_local_factor(1/(pow(constants::pi, 3./2.)) * beta*beta*beta, 0, beta) {}

  cmpl eval(real r, real rp, const Channel& ch) const final {
    return local_potential->eval(r, ch) * non_local_factor( (r-rp), ch);
  };
};

}

#endif 
