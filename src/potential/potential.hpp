#ifndef POTENTIAL_HEADER
#define POTENTIAL_HEADER

#include <memory>
#include <complex>

#include "util/types.hpp"
#include "solver/channel.hpp"

#include "potential/params_base.hpp"
#include "potential/params.hpp"

namespace omplib {

/// @brief interface for general local or non-local potential
struct GeneralPotential {
  using ptr = std::unique_ptr<GeneralPotential>;
  using cptr = std::unique_ptr<const GeneralPotential>;
  
  /// @brief  evaluate the potential at r, rp
  virtual cmpl eval(real r, real rp, const PData& data) const = 0;
  /// @brief  evaluate the reduced potential at r, rp
  /// @returns V(r) * r * delta(r - rp) if potetial is local, 
  /// r * V(r,rp) * rp if potential is non-local
  virtual cmpl eval_reduced(real r, real rp, const PData& data) const = 0;
};

///@brief Base class for an instance of a local potential between
/// some target (A,Z) and a projectile
struct Potential : public GeneralPotential {
  using ptr = std::unique_ptr<Potential>;
  using cptr = std::unique_ptr<const Potential>;

  virtual cmpl eval(real r, const PData& d) const = 0;
  cmpl eval_reduced(real r, const PData& d) const { return eval(r, d) * r; }
  
  cmpl eval(real r, real rp, const PData& d) const final { 
    if (r!=rp) return {0,0};
    return eval(r, d); 
  }
  cmpl eval_reduced(real r, real rp, const PData& d) const final { 
    if (r!=rp) return {0,0};
    return eval_reduced(r,d); 
  }
  
  cmpl operator () (real r, const PData& d) const { return eval(r, d); }
  cmpl operator () (real r, real rp, const PData& d) const { return eval(r, rp, d); }
};

///@brief Base class for an instance of a non-local potential between
/// some target (A,Z) and a projectile
struct NonLocalPotential : public GeneralPotential {
  using ptr = std::unique_ptr<NonLocalPotential>;
  using cptr = std::unique_ptr<const NonLocalPotential>;

  /// @brief r is a distance between the projectile and target
  /// @brief rp is another distance between the projectile and target
  cmpl eval_reduced(real r, real rp, const PData& d) const final { 
    return eval(r, rp, d) * r * rp; 
  }
  cmpl operator () (real r, real rp, const PData& d) const { 
    return eval(r, rp, d); }
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
  
  cmpl eval(real r, const PData&) const final {
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
  
  cmpl eval(real r, const PData&) const final {
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
  
  cmpl eval(real r, const PData& d) const final {
    return DerivWoodsSaxon(V,R,a).eval(r,d)/r;
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
  
  cmpl eval(real r, const PData&) const final {
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
  
  
  cmpl eval(real r, const PData& d) const final {
    cmpl v = 0;
    for (unsigned int i = 0; i < N; ++i) {
      v += gaussians[i].eval(r, d);
    }
    return v;
  };
};

using doubleGaussian = NGaussian<2>;

///@brief Thompson, D. R., LeMere, M., & Tang, Y. C. (1977). 
/// Systematic investigation of scattering problems with the resonating-group method. 
/// Nuclear Physics A, 286(1), 53-66.
class Minnesota : public doubleGaussian {
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
    doubleGaussian({Gaussian(V01, 0., 1/sqrt(k1) ), Gaussian(V02, 0., 1/sqrt(k2) )}) {}
};

class Yukawa : public Potential {
private:
  real mass_coupling;
  real force_coupling;
public:
  cmpl eval(real r, const PData&) const final {
    return - force_coupling * exp(-r/mass_coupling) / r;
  };
  Yukawa(real mass_coupling, real force_coupling)
    : mass_coupling(mass_coupling), force_coupling(force_coupling) {};
  Yukawa(real alpha, real m, real g)
    : mass_coupling(alpha * m), force_coupling(g*g) {};
};

class SphereWell : public Potential {
private: 
  real R;
  cmpl depth;
public:
  cmpl eval(real r, const PData&) const final {
    if (r < R) return -depth;
    return 0;
  };

  SphereWell(real R, cmpl depth)
    : R(R), depth(depth) {};
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
  OMP(Isotope isotope)
    : A(isotope.A), Z(isotope.Z), params(Params()) {};
  OMP(Isotope isotope,  Params params)
    : A(isotope.A), Z(isotope.Z),  params(params) {};
  OMP(Isotope isotope,  Params&& params)
    : A(isotope.A), Z(isotope.Z), params(params) {};

  OMP(const OMP<Params>& rhs) = default;
  OMP(OMP<Params>&& rhs)      = default;

  void reset_target(int An, int Zn) {
    A = An;
    Z = Zn;
  }

  cmpl eval(real r, const PData& d) const final {

    const auto erg_lab = d.e.erg_lab;
    
    const auto V = WoodsSaxon{ 
      cmpl{params.real_cent_V(Z,A,erg_lab),0},
      params.real_cent_r(Z,A,erg_lab),
      params.real_cent_a(Z,A,erg_lab) };
    const auto Vs = DerivWoodsSaxon{ 
      cmpl{params.real_surf_V(Z,A,erg_lab),0},
      params.real_surf_r(Z,A,erg_lab),
      params.real_surf_a(Z,A,erg_lab) };
    const auto Vso = Thomas{ 
      cmpl{params.real_spin_V(Z,A,erg_lab),0},
      params.real_spin_r(Z,A,erg_lab),
      params.real_spin_a(Z,A,erg_lab) };
    
    const auto W = WoodsSaxon{ 
      cmpl{0,params.cmpl_cent_V(Z,A,erg_lab)},
      params.cmpl_cent_r(Z,A,erg_lab),
      params.cmpl_cent_a(Z,A,erg_lab) };
    const auto Ws = DerivWoodsSaxon{ 
      cmpl{0,params.cmpl_surf_V(Z,A,erg_lab)},
      params.cmpl_surf_r(Z,A,erg_lab),
      params.cmpl_surf_a(Z,A,erg_lab) };
    const auto Wso = Thomas{ 
      cmpl{0,params.cmpl_spin_V(Z,A,erg_lab)},
      params.cmpl_spin_r(Z,A,erg_lab),
      params.cmpl_spin_a(Z,A,erg_lab) };

    return V.eval(r,d) + Vs.eval(r,d) + Vso.eval(r,d) * d.am.spin_orbit()  // R
         + W.eval(r,d) + Ws.eval(r,d) + Wso.eval(r,d) * d.am.spin_orbit(); // Im
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

  cmpl eval(real r, real rp, const PData& d) const final {
    return local_potential->eval(r, d) * non_local_factor( (r-rp), d);
  };
};


/// @brief Yamaguchi, Yoshio. 
/// "Two-nucleon problem when the potential is nonlocal but separable. I." 
/// Physical Review 95.6 (1954): 1628.
class Yamaguchi : public NonLocalPotential {
private:
  real alpha;
  real beta;
  real f;
public:
  /// @param mu reduced mass [amu]
  /// @param alpha [fm]^-1
  /// @param beta [fm]^-1
  Yamaguchi(real mu, real alpha, real beta)
    : alpha(alpha)
    , beta(beta)
    , f(
          constants::hbar * constants::hbar * constants::c * constants::c 
        / (mu * constants::MeV_per_amu) 
      ) 
  {};
  
  ///@brief parameters chosen to reproduce bound state of deuteron and 
  /// and neutron-proton triplet scattering length 
  Yamaguchi(): alpha(0.2316053), beta(1.3918324), f(41.472) {};
  cmpl eval(real r, real rp, const PData&) const final {
    return 2 * f * beta * (alpha + beta)*(alpha + beta) * exp(-beta * (r + rp));
  };

  real analytic_swave_kmatrix(real k) const {
    const auto a = alpha;
    const auto b = beta;
    const auto d = 2*(a + b)*(a+b);
    real cot_delta = a*b*(a+2*b)/(d*k) 
                   + (a*a + 2*a*b + 3 * b *b)/(b*d) * k
                   + k*k*k/(b*d);
    return 1./cot_delta;
  }
};

}

#endif 
