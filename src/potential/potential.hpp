#ifndef POTENTIAL_HEADER
#define POTENTIAL_HEADER

#include <memory>
#include <complex>

#include "util/types.hpp"

#include "potential/params_base.hpp"
#include "potential/params.hpp"

namespace omplib {

///@brief Base class for an instance of a local potential between
/// some target (A,Z) and a projectile at a given cms energy
struct Potential {
  using ptr = std::unique_ptr<Potential>;
  using cptr = std::unique_ptr<const Potential>;

  /// @brief r is the distance between the projectile and target
  virtual cmpl eval(real r, real erg) const = 0;
  cmpl operator () (real r, real erg) const { return eval(r, erg); }
};

///@brief Base class for an instance of a non-local potential between
/// some target (A,Z) and a projectile at a given cms energy
struct NonLocalPotential {
  using ptr = std::unique_ptr<NonLocalPotential>;
  using cptr = std::unique_ptr<const NonLocalPotential>;

  /// @brief r is a distance between the projectile and target
  /// @brief rp is another distance between the projectile and target
  virtual cmpl eval(real r, real rp, real erg) const = 0;
  cmpl operator () (real r, real rp, real erg) const { 
    return eval(r,rp, erg); }
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
  
  cmpl eval(real r, real erg) const final {
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
  
  cmpl eval(real r, real erg) const final {
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
  
  cmpl eval(real r, real erg) const final {
    return DerivWoodsSaxon(V,R,a).eval(r,erg)/r;
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
  
  cmpl eval(real r, real erg) const final {
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
  
  
  cmpl eval(real r, real erg) const final {
    cmpl v = 0;
    for (unsigned int i = 0; i < N; ++i) {
      v += gaussians[i].eval(r, erg);
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

  /// @brief dimension of SU(2) repr for orbital angular momentum, 
  /// projectile intrinsic spin, and total projectile angular momentum.
  /// J = L + S, respectvely. 
  /// Dimensions of SU(2) representation w/ spin j is j2 = 2*j+1
  int l2, s2, j2;
  
  /// @brief Determines the parameterization used for calculating 
  /// term depths, radii, and diffusivities
  Params params;

  /// @returns projection of spin along axis of orbital ang. mom.; e.g. L * S
  real spin_orbit() const {
    const real L = (l2 - 1.)/2.;
    const real S = (s2 - 1.)/2.;
    const real J = (l2 - 1.)/2.;

    return 0.5*(J*(J+1) - L*(L+1) - S*(S+1));
  }

public:
  OMP(int A, int Z, int l2, int s2, int j2)
    : A(A), Z(Z), l2(l2), s2(s2), j2(j2), params(Params()) {};
  OMP(int A, int Z, int l2, int s2, int j2, Params params)
    : A(A), Z(Z), l2(l2), s2(s2), j2(j2), params(params) {};
  OMP(int A, int Z, int l2, int s2, int j2, Params&& params)
    : A(A), Z(Z), l2(l2), s2(s2), j2(j2), params(params) {};

  OMP(const OMP<Params>& rhs) = default;
  OMP(OMP<Params>&& rhs)      = default;

  void reset_target(int An, int Zn) {
    A = An;
    Z = An;
  }

  void set_proj_spin(int s2n)   { s2 = s2n; }
  void set_orb_am(int l2n)      { l2 = l2n; }
  void set_proj_tot_am(int j2n) { j2 = j2n; }

  cmpl eval(real r, real erg) const final {
    
    const auto V = WoodsSaxon{ 
      cmpl{params.real_cent_V(Z,A,erg),0},
      params.real_cent_r(Z,A,erg),
      params.real_cent_a(Z,A,erg) };
    const auto Vs = DerivWoodsSaxon{ 
      cmpl{params.real_surf_V(Z,A,erg),0},
      params.real_surf_r(Z,A,erg),
      params.real_surf_a(Z,A,erg) };
    const auto Vso = Thomas{ 
      cmpl{params.real_spin_V(Z,A,erg),0},
      params.real_spin_r(Z,A,erg),
      params.real_spin_a(Z,A,erg) };
    
    const auto W = WoodsSaxon{ 
      cmpl{0,params.cmpl_cent_V(Z,A,erg)},
      params.cmpl_cent_r(Z,A,erg),
      params.cmpl_cent_a(Z,A,erg) };
    const auto Ws = DerivWoodsSaxon{ 
      cmpl{0,params.cmpl_surf_V(Z,A,erg)},
      params.cmpl_surf_r(Z,A,erg),
      params.cmpl_surf_a(Z,A,erg) };
    const auto Wso = Thomas{ 
      cmpl{0,params.cmpl_spin_V(Z,A,erg)},
      params.cmpl_spin_r(Z,A,erg),
      params.cmpl_spin_a(Z,A,erg) };

    return V.eval(r,erg) + Vs.eval(r,erg) + Vso.eval(r,erg) * spin_orbit()  // R
         + W.eval(r,erg) + Ws.eval(r,erg) + Wso.eval(r,erg) * spin_orbit(); // Im
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

  cmpl eval(real r, real rp, real erg) const final {
    return local_potential->eval(r, erg) * non_local_factor( (r-rp), erg);
  };
};

}

#endif 
