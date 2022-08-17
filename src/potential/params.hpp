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
  double e_fermi_0;
  double e_fermi_A;
  
  // real central
  double v1_0, v1_asym, v1_A, v2_0, v2_A, v3_0, v3_A, v4_0; // depth
  double r_0, r_A, a_0, a_A;                                // shape

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
  double Ef(int at) const;
  double asym(int zt, int at) const;

public:
  
  double real_volu_r(int Z, int A, double erg) const final { return 0; }
  double cmpl_volu_r(int Z, int A, double erg) const final { return 0; }
  double real_surf_r(int Z, int A, double erg) const final { return 0; }
  double cmpl_surf_r(int Z, int A, double erg) const final { return 0; }
  double real_spin_r(int Z, int A, double erg) const final { return 0; }
  double cmpl_spin_r(int Z, int A, double erg) const final { return 0; }
  
  double real_volu_a(int Z, int A, double erg) const final { return 0; }
  double cmpl_volu_a(int Z, int A, double erg) const final { return 0; }
  double real_surf_a(int Z, int A, double erg) const final { return 0; }
  double cmpl_surf_a(int Z, int A, double erg) const final { return 0; }
  double real_spin_a(int Z, int A, double erg) const final { return 0; }
  double cmpl_spin_a(int Z, int A, double erg) const final { return 0; }

  double real_volu_V(int Z, int A, double erg) const final { return 0; }
  double cmpl_volu_V(int Z, int A, double erg) const final { return 0; }
  double real_surf_V(int Z, int A, double erg) const final { return 0; }
  double cmpl_surf_V(int Z, int A, double erg) const final { return 0; }
  double real_spin_V(int Z, int A, double erg) const final { return 0; }
  double cmpl_spin_V(int Z, int A, double erg) const final { return 0; }

  KoningDelaroche03<projectile>(json param_file);
  KoningDelaroche03<projectile>();
};
/*

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
*/
};


#endif
