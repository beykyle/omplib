#ifndef PARAMS_HPP
#define PARAMS_HPP

#include "nlohmann/json.hpp"
using nlohmann::json;

namespace omplib {

enum class Projectile : bool {
  proton,
  neutron,
};
 
template<Projectile p>
struct OMParams {
  constexpr static Projectile projectile = p;

  // Woods-Saxon term radii 
  virtual double real_volu_r(int Z, int A, double erg) const = 0;
  virtual double cmpl_volu_r(int Z, int A, double erg) const = 0;
  virtual double real_surf_r(int Z, int A, double erg) const = 0;
  virtual double cmpl_surf_r(int Z, int A, double erg) const = 0;
  virtual double real_spin_r(int Z, int A, double erg) const = 0;
  virtual double cmpl_spin_r(int Z, int A, double erg) const = 0;
  
  // Woods-Saxon term diffusivity
  virtual double real_volu_a(int Z, int A, double erg) const = 0;
  virtual double cmpl_volu_a(int Z, int A, double erg) const = 0;
  virtual double real_surf_a(int Z, int A, double erg) const = 0;
  virtual double cmpl_surf_a(int Z, int A, double erg) const = 0;
  virtual double real_spin_a(int Z, int A, double erg) const = 0;
  virtual double cmpl_spin_a(int Z, int A, double erg) const = 0;

  // Woods-Saxon term depth
  virtual double real_volu_V(int Z, int A, double erg) const = 0;
  virtual double cmpl_volu_V(int Z, int A, double erg) const = 0;
  virtual double real_surf_V(int Z, int A, double erg) const = 0;
  virtual double cmpl_surf_V(int Z, int A, double erg) const = 0;
  virtual double real_spin_V(int Z, int A, double erg) const = 0;
  virtual double cmpl_spin_V(int Z, int A, double erg) const = 0;
};

template<>
struct OMParams<Projectile::proton> : public OMParams<Projectile::neutron>{

  virtual double real_coul_r(int Z, int A, double erg) const = 0;
  virtual double real_coul_V(int Z, int A, double erg) const = 0;
};


// A. Koning and J. Delaroche, 
// Nuclear Physics A 713, 231 (2003), ISSN 0375-9474, 
// URL https://www.sciencedirect.com/science/article/pii/S0375947402013210.
template <Projectile projectile>
class KoningDelaroche03 : public OMParams<projectile> {
  // fermi energy
  double e_fermi_0, e_fermi_A;
  
  // real central
  double v1_0, v1_asym, v1_A, v2_0, v2_A, v3_0, v3_A, v4_0; // depth
  double rv_0, rv_A, av_0, av_A;                                // shape

  // cmplex central
  double w1_0, w1_A, w2_0, w2_A;
  
  // cmplex surface
  double d1_0, d1_asym, d2_0, d2_A, d2_A2, d2_A3, d3_0;
  double rd_0, rd_A, ad_0, ad_A;

  // real spin orbit
  double vso1_0, vso1_A, vso2_0;
  double rso_0, rso_A, aso_0;
  
  // cmplex spin orbit
  double wso1, wso2;

  // structure and energy factors
  double Ef(int A) const;
  double asym(int Z, int A) const;

public:
  
  double real_volu_r(int Z, int A, double erg) const final;
  double cmpl_volu_r(int Z, int A, double erg) const final;
  double real_surf_r(int Z, int A, double erg) const final;
  double cmpl_surf_r(int Z, int A, double erg) const final;
  double real_spin_r(int Z, int A, double erg) const final;
  double cmpl_spin_r(int Z, int A, double erg) const final;
  
  double real_volu_a(int Z, int A, double erg) const final;
  double cmpl_volu_a(int Z, int A, double erg) const final;
  double real_surf_a(int Z, int A, double erg) const final;
  double cmpl_surf_a(int Z, int A, double erg) const final;
  double real_spin_a(int Z, int A, double erg) const final;
  double cmpl_spin_a(int Z, int A, double erg) const final;

  double real_volu_V(int Z, int A, double erg) const final;
  double cmpl_volu_V(int Z, int A, double erg) const final;
  double real_surf_V(int Z, int A, double erg) const final;
  double cmpl_surf_V(int Z, int A, double erg) const final;
  double real_spin_V(int Z, int A, double erg) const final;
  double cmpl_spin_V(int Z, int A, double erg) const final;

  KoningDelaroche03<projectile>(json param_file);
  KoningDelaroche03<projectile>();
};

template<>
class KoningDelaroche03<Projectile::proton> : public OMParams<Projectile::proton>{
  double rc_0, rc_A, rc_A2;

  double real_coul_r(int Z, int A, double erg) const final;
  double real_coul_V(int Z, int A, double erg) const final;
  
  KoningDelaroche03(json param_file);
  KoningDelaroche03();
};

// KD03 - term params that are the same between projectiles
template<Projectile p>
double omplib::KoningDelaroche03<p>::Ef(int A) const {
  return e_fermi_0 + e_fermi_A * static_cast<double>(A);
};

template<Projectile p>
double omplib::KoningDelaroche03<p>::asym(int Z, int A) const {
  const double a = static_cast<double>(A);
  const double z = static_cast<double>(Z);
  const double n = a - z;
  return (n - z)/a;
}

template<Projectile p>
double omplib::KoningDelaroche03<p>::real_volu_r(
    int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return rv_0 - rv_A / pow(A , -1./3.);
};

template<Projectile p>
double omplib::KoningDelaroche03<p>::cmpl_volu_r(
    int Z, int A, double erg) const {
  return real_volu_r(Z,A,erg);
};

template<Projectile p>
double omplib::KoningDelaroche03<p>::real_surf_r(
    int Z, int A, double erg) const {
  return 0;
};

template<Projectile p>
double omplib::KoningDelaroche03<p>::cmpl_surf_r(
    int Z, int A, double erg) const {
  return rd_0 - rd_A * pow(static_cast<double>(A), -1./3.);;
};

template<Projectile p>
double omplib::KoningDelaroche03<p>::real_spin_r(
    int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return rso_0 - rso_A * pow(a, -1./3.);
};

template<Projectile p>
double omplib::KoningDelaroche03<p>::cmpl_spin_r(
    int Z, int A, double erg) const {
  return real_spin_r(Z,A,erg);
};

template<Projectile p>
double omplib::KoningDelaroche03<p>::real_volu_a(
    int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return av_0 - av_A / a;
};

template<Projectile p>
double omplib::KoningDelaroche03<p>::cmpl_volu_a(
    int Z, int A, double erg) const {
  return real_volu_a(Z,A,erg);
};

template<Projectile p>
double omplib::KoningDelaroche03<p>::real_surf_a(
    int Z, int A, double erg) const {
  return 0;
};

template<Projectile p>
double omplib::KoningDelaroche03<p>::cmpl_surf_a(
    int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  if constexpr (p == Projectile::neutron)
    return ad_0 - ad_A * a;
  else if constexpr(p == Projectile::proton)
    return ad_0 + ad_A * a;
};

template<Projectile p>
double omplib::KoningDelaroche03<p>::real_spin_a(
    int Z, int A, double erg) const {
  return aso_0;
};

template<Projectile p>
double omplib::KoningDelaroche03<p>::cmpl_spin_a(
    int Z, int A, double erg) const {
  return real_spin_a(Z,A,erg);
};

// R. Varner, W. Thompson, T. McAbee, E. Ludwig, and T. Clegg, 
// Physics Reports 201, 57 (1991), ISSN 0370-1573, 
// URL https://www.sciencedirect.com/science/article/pii/037015739190039O
template <Projectile projectile>
class ChapelHill89 : public OMParams<projectile> {
public:
  
  double real_volu_r(int Z, int A, double erg) const final;
  double cmpl_volu_r(int Z, int A, double erg) const final;
  double real_surf_r(int Z, int A, double erg) const final;
  double cmpl_surf_r(int Z, int A, double erg) const final;
  double real_spin_r(int Z, int A, double erg) const final;
  double cmpl_spin_r(int Z, int A, double erg) const final;
  
  double real_volu_a(int Z, int A, double erg) const final;
  double cmpl_volu_a(int Z, int A, double erg) const final;
  double real_surf_a(int Z, int A, double erg) const final;
  double cmpl_surf_a(int Z, int A, double erg) const final;
  double real_spin_a(int Z, int A, double erg) const final;
  double cmpl_spin_a(int Z, int A, double erg) const final;

  double real_volu_V(int Z, int A, double erg) const final;
  double cmpl_volu_V(int Z, int A, double erg) const final;
  double real_surf_V(int Z, int A, double erg) const final;
  double cmpl_surf_V(int Z, int A, double erg) const final;
  double real_spin_V(int Z, int A, double erg) const final;

  ChapelHill89<projectile>(json param_file);
  ChapelHill89<projectile>();
};


// T. R. Whitehead, Y. Lim, and J. W. Holt, 
// Phys. Rev. Lett. 127, 182502 (2021), 
// URL https://link.aps.org/doi/10.1103/PhysRevLett.127.182502.
template <Projectile projectile>
class WLH21 : public OMParams<projectile> {
public:
  
  double real_volu_r(int Z, int A, double erg) const final;
  double cmpl_volu_r(int Z, int A, double erg) const final;
  double real_surf_r(int Z, int A, double erg) const final;
  double cmpl_surf_r(int Z, int A, double erg) const final;
  double real_spin_r(int Z, int A, double erg) const final;
  double cmpl_spin_r(int Z, int A, double erg) const final;
  
  double real_volu_a(int Z, int A, double erg) const final;
  double cmpl_volu_a(int Z, int A, double erg) const final;
  double real_surf_a(int Z, int A, double erg) const final;
  double cmpl_surf_a(int Z, int A, double erg) const final;
  double real_spin_a(int Z, int A, double erg) const final;
  double cmpl_spin_a(int Z, int A, double erg) const final;

  double real_volu_V(int Z, int A, double erg) const final;
  double cmpl_volu_V(int Z, int A, double erg) const final;
  double real_surf_V(int Z, int A, double erg) const final;
  double cmpl_surf_V(int Z, int A, double erg) const final;
  double real_spin_V(int Z, int A, double erg) const final;

  WLH21<projectile>(json param_file);
  WLH21<projectile>();
};
};


#endif
