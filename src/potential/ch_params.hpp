#ifndef CH_PARAMS_HEADER
#define CH_PARAMS_HEADER

#include "util/types.hpp"

#include "potential/params_base.hpp"

namespace omplib {

/// @brief  Phenomenological global OM potential parameterization from:
/// R. Varner, W. Thompson, T. McAbee, E. Ludwig, and T. Clegg, 
/// Physics Reports 201, 57 (1991), ISSN 0370-1573, 
/// URL https://www.sciencedirect.com/science/article/pii/037015739190039O
template <Proj projectile>
class CH89Params : public OMParams<projectile> {
protected:
  
  // real central shape
  real r_0, r_A;
  real a0;
  
  // complex central and surface shape
  real rw_0, rw_A;
  real aw;
  
  // real spin orbit shape
  real rso_0, rso_A;
  real aso;
  
  // real central depth
  real v_0, v_e, v_asym;
  
  // complex central depth
  real wv_0, wve_0, wv_ew;
  
  // complex surface depth
  real ws_0, ws_asym, ws_e0, ws_ew;

  // real spin orbit depth
  real vso_0;

  real Ec(int, int, real) const { return 0; }
  
  using OMParams<projectile>::asym;

public:
  
  real real_cent_r(int Z, int A, real erg) const final;
  real cmpl_cent_r(int Z, int A, real erg) const final;
  real cmpl_surf_r(int Z, int A, real erg) const final;
  real real_spin_r(int Z, int A, real erg) const final;
  
  real real_cent_a(int Z, int A, real erg) const final;
  real cmpl_cent_a(int Z, int A, real erg) const final;
  real cmpl_surf_a(int Z, int A, real erg) const final;
  real real_spin_a(int Z, int A, real erg) const final;

  real real_cent_V(int Z, int A, real erg) const final;
  real cmpl_cent_V(int Z, int A, real erg) const final;
  real cmpl_surf_V(int Z, int A, real erg) const final;
  real real_spin_V(int Z, int A, real erg) const final;

  // CH89 does not have real surface or complex spin terms
  real cmpl_spin_r(int, int, real) const final { return 0; }
  real cmpl_spin_V(int, int, real) const final { return 0; }
  real cmpl_spin_a(int, int, real) const final { return 0; }
  
  real real_surf_a(int, int, real) const override { return 0; }
  real real_surf_V(int, int, real) const override { return 0; }
  real real_surf_r(int, int, real) const override { return 0; }
  
  CH89Params( const CH89Params<projectile>& rhs ) = default;
  
  // construct using default CH89 params
  CH89Params();

  CH89Params(json p):
    r_0(     p["CH89RealCentral_r_o_0"] ) , 
    r_A(     p["CH89RealCentral_r_o"]   ) ,
    a0(      p["CH89RealCentral_a_0"]   ) , 
    rw_0(    p["CH89ImagCentral_r_w0"]  ),
    rw_A(    p["CH89ImagCentral_r_w"]   ),
    aw(      p["CH89ImagCentral_a_w"]   ),
    rso_0(   p["CH89SpinOrbit_r_so_0"]   ),
    rso_A(   p["CH89SpinOrbit_r_so"]     ),
    aso(     p["CH89SpinOrbit_a_so"]     ),
    
    v_0(     p["CH89RealCentral_V_0"]   ) ,
    v_e(     p["CH89RealCentral_V_e"]   ) ,
    v_asym(  p["CH89RealCentral_V_t"]   ) ,
    wv_0(    p["CH89ImagCentral_W_v0"]  ),
    wve_0(   p["CH89ImagCentral_W_ve0"] ),
    wv_ew(   p["CH89ImagCentral_W_vew"] ),
    ws_0(    p["CH89ImagCentral_W_s0"]  ),
    ws_asym( p["CH89ImagCentral_W_st"]  ),
    ws_e0(   p["CH89ImagCentral_W_se0"] ),
    ws_ew(   p["CH89ImagCentral_W_sew"] ),
    vso_0(   p["CH89SpinOrbit_V_so"]     )
    
  {}

  /// @brief constructs a CH89Params\<p\> with params refit w/ MCMC; from
  /// Pruitt, C. D. et al, 
  /// “Uncertainty-Quantified Phenomenological Optical Potentials 
  /// for Single-Nucleon Scattering”, 
  /// LLNL release number LLNL-JRNL-835671-DRAFT (to be published).
  static CH89Params<projectile> build_CHUQ()
  {
    auto p = CH89Params<projectile>{};
    
    p.v_0     = 56.19;
    p.v_asym  = 13.82;
    p.v_e     = -0.36;
    p.r_0     = -0.20;
    p.r_A     = 1.20;
    p.a0      = 0.73;
    p.vso_0   = 5.58;
    p.rso_0   = -1.12;
    p.rso_A   = 1.29;
    p.aso     = 0.61;
    p.wv_0    = 9.92;
    p.wve_0   = 33.15;
    p.wv_ew   = 24.0;
    p.ws_0    = 10.59;
    p.ws_asym = 27.09;
    p.ws_e0   = 20.00;
    p.ws_ew   = 36.38;
    p.rw_0    = -0.41;
    p.rw_A    = 1.32;
    p.aw      = 0.69;

    if constexpr (projectile == Proj::proton) {
      p.rc_0 = 0.13;
      p.rc_A = 1.25;
    }

    return p;
  };
};

template<>
class CH89Params<Proj::proton> : 
  public CH89Params<Proj::neutron> , OMParams<Proj::proton> {
protected: 
  real rc_0, rc_A;
  real Ec(int Z, int A, real erg) const;

public:
  constexpr static Proj projectile = Proj::proton;
  real real_coul_r(int Z, int A, real erg) const final;
  
  CH89Params( 
      const CH89Params<Proj::proton>& rhs) = default;
  CH89Params();
  CH89Params(json p);
};


template<Proj proj>
real CH89Params<proj>::real_cent_r(int, int A, real) const {
  const real a = static_cast<real>(A);
  return r_0 + r_A * pow(a, 1./3.);
}

template<Proj proj>
real CH89Params<proj>::cmpl_cent_r(int, int A, real) const {
  const real a = static_cast<real>(A);
  return rw_0 + rw_A * pow(a, 1./3.);
}

template<Proj proj>
real CH89Params<proj>::cmpl_surf_r(int Z, int A, real erg) const {
  return cmpl_cent_r(Z,A,erg);
}

template<Proj proj>
real CH89Params<proj>::real_spin_r(int, int A, real) const {
  const real a = static_cast<real>(A);
  return rso_0 + rso_A * pow(a, 1./3.);
}

template<Proj proj>
real CH89Params<proj>::real_cent_a(int, int, real) const {
  return a0;
}

template<Proj proj>
real CH89Params<proj>::cmpl_cent_a(int, int, real) const {
  return aw;
}

template<Proj proj>
real CH89Params<proj>::cmpl_surf_a(int, int, real) const {
  return aw;
}

template<Proj proj>
real CH89Params<proj>::real_spin_a(int, int, real) const {
  return aso;
}

template<Proj proj>
real CH89Params<proj>::real_cent_V(int Z, int A, real erg) const {
  const real dE = erg - Ec(Z,A,erg);
  return -(v_0 + v_e * dE + asym(Z,A) * v_asym);
}

template<Proj proj>
real CH89Params<proj>::cmpl_cent_V(int Z, int A, real erg) const {
  const real dE = erg - Ec(Z,A,erg);
  return -wv_0 / (1 + exp( (wve_0 - dE)/wv_ew ));
}

template<Proj proj>
real CH89Params<proj>::cmpl_surf_V(int Z, int A, real erg) const {
  const real dE = erg - Ec(Z,A,erg);
  const real Ws = (ws_0  + asym(Z,A) * ws_asym) / (1 + exp( (dE - ws_e0)/ws_ew ));
  return 4 * cmpl_surf_a(Z,A,erg) * Ws;
}

template<Proj proj>
real CH89Params<proj>::real_spin_V(int, int, real) const {
  return 2. * vso_0;
}

}
#endif
