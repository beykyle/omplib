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
public:
  
  double real_cent_r(int Z, int A, double erg) const final;
  double cmpl_cent_r(int Z, int A, double erg) const final;
  double real_surf_r(int Z, int A, double erg) const final;
  double cmpl_surf_r(int Z, int A, double erg) const final;
  double real_spin_r(int Z, int A, double erg) const final;
  double cmpl_spin_r(int Z, int A, double erg) const final;
  
  double real_cent_a(int Z, int A, double erg) const final;
  double cmpl_cent_a(int Z, int A, double erg) const final;
  double real_surf_a(int Z, int A, double erg) const final;
  double cmpl_surf_a(int Z, int A, double erg) const final;
  double real_spin_a(int Z, int A, double erg) const final;
  double cmpl_spin_a(int Z, int A, double erg) const final;

  double real_cent_V(int Z, int A, double erg) const final;
  double cmpl_cent_V(int Z, int A, double erg) const final;
  double real_surf_V(int Z, int A, double erg) const final;
  double cmpl_surf_V(int Z, int A, double erg) const final;
  double real_spin_V(int Z, int A, double erg) const final;

  WLH21<projectile>(json param_file);
  WLH21<projectile>();
};

}

#endif 
