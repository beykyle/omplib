#ifndef CH_PARAMS_HEADER
#define CH_PARAMS_HEADER

#include "potential/params_base.hpp"

namespace omplib {

/// @brief  Phenomenological global OM potential parameterization from:
/// R. Varner, W. Thompson, T. McAbee, E. Ludwig, and T. Clegg, 
/// Physics Reports 201, 57 (1991), ISSN 0370-1573, 
/// URL https://www.sciencedirect.com/science/article/pii/037015739190039O
template <Proj projectile>
class ChapelHill89 : public OMParams<projectile> {
protected:
  
  // real central shape
  double r_0, r_A;
  double a0;
  
  // complex central and surface shape
  double rw_0, rw_A;
  double aw;
  
  // real spin orbit shape
  double rso_0, rso_A;
  double aso;
  
  // real central depth
  double v_0, v_e, v_asym;
  
  // complex central depth
  double wv_0, wve_0, wv_ew;
  
  // complex surface depth
  double ws_0, ws_asym, ws_e0, ws_ew;

  // real spin orbit depth
  double vso_0;

  double Ec(int Z, int A, double erg) const { return 0; }
  
  using OMParams<projectile>::asym;

public:
  
  double real_cent_r(int Z, int A, double erg) const final;
  double cmpl_cent_r(int Z, int A, double erg) const final;
  double cmpl_surf_r(int Z, int A, double erg) const final;
  double real_spin_r(int Z, int A, double erg) const final;
  
  double real_cent_a(int Z, int A, double erg) const final;
  double cmpl_cent_a(int Z, int A, double erg) const final;
  double cmpl_surf_a(int Z, int A, double erg) const final;
  double real_spin_a(int Z, int A, double erg) const final;

  double real_cent_V(int Z, int A, double erg) const final;
  double cmpl_cent_V(int Z, int A, double erg) const final;
  double cmpl_surf_V(int Z, int A, double erg) const final;
  double real_spin_V(int Z, int A, double erg) const final;

  // CH89 does not have real surface or complex spin terms
  double cmpl_spin_r(int Z, int A, double erg) const final { return 0; }
  double cmpl_spin_V(int Z, int A, double erg) const final { return 0; }
  double cmpl_spin_a(int Z, int A, double erg) const final { return 0; }
  
  double real_surf_a(int Z, int A, double erg) const override { return 0; }
  double real_surf_V(int Z, int A, double erg) const override { return 0; }
  double real_surf_r(int Z, int A, double erg) const override { return 0; }
  
  ChapelHill89( const ChapelHill89<projectile>& rhs ) = default;
  
  // construct using default CH89 params
  ChapelHill89();

  ChapelHill89(json p):
    v_0(     p["CH89RealCentral"]["V_0"] ) ,
    v_e(     p["CH89RealCentral"]["V_e"] ) ,
    v_asym(  p["CH89RealCentral"]["V_t"] ) ,
    r_0(     p["KDHartreeFock"]["r_o_0"] ) , 
    r_A(     p["KDImagVolume"]["r_o"]    ) ,
    a0(      p["KDImagSurface"]["a_0"]   ) , 
    
    wv_0(    p["CH89ImageCentral"]["W_v0"]  ),
    wve_0(   p["CH89ImageCentral"]["W_ve0"] ),
    wv_ew(   p["CH89ImageCentral"]["W_vew"] ),
    ws_0(    p["CH89ImageCentral"]["W_s0"]  ),
    ws_e0(   p["CH89ImageCentral"]["W_se0"] ),
    ws_ew(   p["CH89ImageCentral"]["W_sew"] ),
    ws_asym( p["CH89ImageCentral"]["W_st"]  ),
    rw_0(    p["CH89ImageCentral"]["r_w0"]  ),
    rw_A(    p["CH89ImageCentral"]["r_w"]   ),
    aw(      p["CH89ImageCentral"]["a_w"]   ),
    
    vso_0(   p["CH89SpinOrbit"]["V_so"]     ),
    rso_0(   p["CH89SpinOrbit"]["r_so"]     ),
    rso_A(   p["CH89SpinOrbit"]["r_so_0"]   ),
    aso(     p["CH89SpinOrbit"]["a_so"]     )
    
  {}

  /// @brief constructs a ChapelHill89\<p\> with params refit w/ MCMC; from
  /// Pruitt, C. D. et al, 
  /// “Uncertainty-Quantified Phenomenological Optical Potentials 
  /// for Single-Nucleon Scattering”, 
  /// LLNL release number LLNL-JRNL-835671-DRAFT (to be published).
  static ChapelHill89<projectile> build_CHUQ()
  {
    auto p = ChapelHill89<projectile>{};
    
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
class ChapelHill89<Proj::proton> : 
  public ChapelHill89<Proj::neutron> , OMParams<Proj::proton> {
protected: 
  double rc_0, rc_A;
  double Ec(int Z, int A, double erg) const;

public:
  constexpr static Proj projectile = Proj::proton;
  double real_coul_r(int Z, int A, double erg) const final;
  double real_coul_V(int Z, int A, double erg) const final;
  
  ChapelHill89( 
      const ChapelHill89<Proj::proton>& rhs) = default;
  ChapelHill89();
  ChapelHill89(json p);
};


template<Proj proj>
double ChapelHill89<proj>::real_cent_r(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return r_0 + a * r_A;
}

template<Proj proj>
double ChapelHill89<proj>::cmpl_cent_r(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return rw_0 + a * rw_A;
}

template<Proj proj>
double ChapelHill89<proj>::cmpl_surf_r(int Z, int A, double erg) const {
  return cmpl_cent_r(Z,A,erg);
}

template<Proj proj>
double ChapelHill89<proj>::real_spin_r(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return rso_0 + a * rso_A;
}

template<Proj proj>
double ChapelHill89<proj>::real_cent_a(int Z, int A, double erg) const {
  return a0;
}

template<Proj proj>
double ChapelHill89<proj>::cmpl_cent_a(int Z, int A, double erg) const {
  return aw;
}

template<Proj proj>
double ChapelHill89<proj>::cmpl_surf_a(int Z, int A, double erg) const {
  return aw;
}

template<Proj proj>
double ChapelHill89<proj>::real_spin_a(int Z, int A, double erg) const {
  return aso;
}

template<Proj proj>
double ChapelHill89<proj>::real_cent_V(int Z, int A, double erg) const {
  const double dE = erg - Ec(Z,A,erg);
  return v_0 + v_e * dE + asym(Z,A) * v_asym;
}

template<Proj proj>
double ChapelHill89<proj>::cmpl_cent_V(int Z, int A, double erg) const {
  const double dE = erg - Ec(Z,A,erg);
  return wv_0 / (1 + exp( (wve_0 - dE)/wv_ew ));
}

template<Proj proj>
double ChapelHill89<proj>::cmpl_surf_V(int Z, int A, double erg) const {
  const double dE = erg - Ec(Z,A,erg);
  return (ws_0  + asym(Z,A) * ws_asym) / (1 + exp( (dE - ws_e0)/ws_ew ));
}

template<Proj proj>
double ChapelHill89<proj>::real_spin_V(int Z, int A, double erg) const {
  return vso_0;
}

}
#endif
