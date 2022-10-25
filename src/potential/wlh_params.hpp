#ifndef WLH_PARAMS_HEADER
#define WLH_PARAMS_HEADER

#include "potential/params_base.hpp"
#include "util/constants.hpp"

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
  WLH21():
  v0(    52.6912521913),  
  v1(    0.2849592984),   
  v2(    -0.0002654968),  
  v3(    2.6234e-06),     
  v4(    21.0895801061),  
  v5(    0.2847889774),   
  v6(    0.0010745253),   
  r0(    1.2978610209),   
  r1(    0.312492324),    
  r2(    0.0008999731),   
  r3(    8.6727e-06),     
  a0(    0.741279032),    
  a1(    0.0008329878),   
  a2(    1.02475e-05),    
  a3(    0.2792868326),   
  a4(    0.8764603676),   
  w0(    2.8722965665), 
  w1(    0.2480075276), 
  w2(    0.0004858014), 
  w3(    8.7072302176), 
  w4(    0.0246183387), 
  rw0(   0.5629863841),
  rw1(   79.8221560535),
  rw2(   0.7264601343),
  rw3(   86.5142396423),
  rw4(   0.9175065953),
  rw5(   3.6419e-06),
  aw0(   0.2391198977),
  aw1(   0.5849807386),
  aw2(   6.7877182692),
  aw3(   0.0529392097),
  aw4(   0.0005803378),
  d0(    1.6408748071),
  d1(    0.0395001751),
  d2(    2.355450576),
  d3(    0.1418364594),
  rs0(   1.2882811754),
  rs1(   1.4774733457),
  rs2(   0.0030997663),
  as0(   0.8506076645),
  vso_0( 9.6220335001),
  vso_1( 0.0085454732),
  rso_0( 1.2794000707),
  rso_1( 0.8734769907),
  aso_0( 0.8060570111),
  aso_1( 0.0003509748)
  {}

  using OMParams<projectile>::asym;
};

template<>
class WLH21<Proj::proton> : public WLH21<Proj::neutron> , OMParams<Proj::proton> {
public:
  //TODO what is the Coulomb contribution to WLH
  double real_coul_r(int Z, int A, double erg) const final { return 0; }
  WLH21(json p);
  WLH21();
};


template<Proj proj>
double WLH21<proj>::real_cent_r(int Z, int A, double erg) const {
  const double a  = static_cast<double>(A);
  const double a3 = pow(a,1./3.); 
  return (r0 - r1 * erg + r2 * erg * erg)*a3  - r3;
}

template<Proj proj>
double WLH21<proj>::cmpl_cent_r(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  const double a3 = pow(a,1./3.); 
  return a3*(rw0 + (rw1 + rw2 * a)/(rw3 + a + rw4 * erg) + rw5 * erg * erg);
}

template<Proj proj>
double WLH21<proj>::cmpl_surf_r(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  const double a3 = pow(a,1./3.); 
  return a3*(rs0 - rs1 * erg) - rs2;
}

template<Proj proj>
double WLH21<proj>::real_spin_r(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  const double a3 = pow(a,1./3.); 
  return rso_0 * a3 - rso_1;
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
  if constexpr (proj ==Proj::proton) {
    if (erg > 20) return 0;
  }
  if constexpr (proj ==Proj::neutron) {
    if (erg > 40) return 0;
  }
  // delta = +(N-Z)/A
  const double delta = OMParams<Proj::proton>::asym(Z,A);
  const double Uso = d0 - d1 * erg - (d2  - d3 * erg) * delta;
  return Uso;
}

template<Proj proj>
double WLH21<proj>::real_spin_V(int Z, int A, double erg) const {
  return vso_0 - vso_1 * static_cast<double>(A);
}

}

#endif 
