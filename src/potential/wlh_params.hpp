#ifndef WLH_PARAMS_HEADER
#define WLH_PARAMS_HEADER

#include "potential/params_base.hpp"

namespace omplib {

/// @brief  Microscopic global OM potential parameterization using XEFT, from
/// T. R. Whitehead, Y. Lim, and J. W. Holt, 
/// Phys. Rev. Lett. 127, 182502 (2021), 
/// URL https://link.aps.org/doi/10.1103/PhysRevLett.127.182502.
template <Proj projectile>
class WLH21 : public OMParams<projectile> {
protected:
  // real central depth
  double v0, v1, v2, v3, v4, v5, v6;

  // real central shape
  double r0, r1, r2, r3;
  double a0, a1, a2, a3, a4;
  
  // complex central depth
  double w0, w1, w2, w3, w4;

  // complex central shape
  double rw0, rw1, rw2, rw3, rw4, rw5;
  double aw0, aw1, aw2, aw3, aw4;
  
  // complex surface depth
  double d0, d1, d2, d3;

  // complex surface shape
  double rs0, rs1, rs2;
  double as0;
  
  // real spn-orbit depth
  double vso_0, vso_1;

  // real spin-orbit shape
  double rso_0, rso_1;
  double aso_0, aso_1;

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
  
  // WLH does not have a real surface or complex SO term
  double real_surf_V(int Z, int A, double erg) const final { return 0; }
  double real_surf_a(int Z, int A, double erg) const final { return 0; }
  double real_surf_r(int Z, int A, double erg) const final { return 0; }

  double cmpl_spin_V(int Z, int A, double erg) const final { return 0; }
  double cmpl_spin_a(int Z, int A, double erg) const final { return 0; }
  double cmpl_spin_r(int Z, int A, double erg) const final { return 0; }

  WLH21(json p);
  WLH21();

  using OMParams<projectile>::asym;
};

template<>
class WLH21<Proj::proton> : public WLH21<Proj::neutron> , OMParams<Proj::proton> {
public:
  WLH21(json p);
  WLH21();
};


template<Proj proj>
double WLH21<proj>::real_cent_r(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return r0 - r1 * erg + r2 * erg * erg - r3 * pow(a,-1./3.);
}

template<Proj proj>
double WLH21<proj>::cmpl_cent_r(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return rw0 + (rw1 + rw2 * a)/(rw3 + a + rw4 * erg) + rw5 * erg * erg;
}

template<Proj proj>
double WLH21<proj>::cmpl_surf_r(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return rs0 - rs1 * erg - rs2 * pow(a, -1./3.);
}

template<Proj proj>
double WLH21<proj>::real_spin_r(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return rso_0 - rso_1 * pow(a,-1./3.);
}

template<Proj proj>
double WLH21<proj>::real_cent_a(int Z, int A, double erg) const {
  // delta = +(N-Z)/A
  const double delta = OMParams<Proj::proton>::asym(Z,A);
  const double a_np =  a0 - a2 * erg * erg - (a3 - a4 * delta ) * delta;
  if constexpr (proj == Proj::neutron)
    return a_np - a1 * erg;
  else if constexpr (proj == Proj::proton)
    return a_np + a1 * erg;
}

template<Proj proj>
double WLH21<proj>::cmpl_cent_a(int Z, int A, double erg) const {
  // delta = +(N-Z)/A
  const double delta = OMParams<Proj::proton>::asym(Z,A);
  return aw0 + aw1 * erg / (aw2 + erg) + (aw3 - aw4 * erg) * delta;
}

template<Proj proj>
double WLH21<proj>::cmpl_surf_a(int Z, int A, double erg) const {
  return as0;
}

template<Proj proj>
double WLH21<proj>::real_spin_a(int Z, int A, double erg) const {
  return aso_0 - aso_1 * static_cast<double>(A);
}

template<Proj proj>
double WLH21<proj>::real_cent_V(int Z, int A, double erg) const {
  // delta = +/- (N-Z)/A ; - for neutron, + for proton
  const double delta  = asym(Z,A);
  const double v_erg  = v0 - v1 * erg + v2 * erg * erg + v3 * erg * erg * erg;
  const double v_asym = ( v4 - v5  * erg + v6 * erg * erg ) * delta;
  return v_erg + v_asym;
}

template<Proj proj>
double WLH21<proj>::cmpl_cent_V(int Z, int A, double erg) const {
  // delta = +/- (N-Z)/A ; - for neutron, + for proton
  const double delta  = asym(Z,A);
  const double v_erg = w0 + w1 * erg - w2 * erg * erg;
  if constexpr (proj == Proj::neutron)
    return v_erg + (w3 - w4 * erg ) * delta;
  else if constexpr (proj == Proj::proton)
    return v_erg + (-w3 - w4 * erg ) * delta;
}

template<Proj proj>
double WLH21<proj>::cmpl_surf_V(int Z, int A, double erg) const {
  // delta = +(N-Z)/A
  const double delta = OMParams<Proj::proton>::asym(Z,A);
  return d0 - d1 * erg - (d2  - d3 * erg) * delta;
}

template<Proj proj>
double WLH21<proj>::real_spin_V(int Z, int A, double erg) const {
  return vso_0 - vso_1 * static_cast<double>(A);
}

}

#endif 
