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
public:
  
  double real_cent_r(int Z, int A, double erg) const final;
  double cmpl_cent_r(int Z, int A, double erg) const final;
  double cmpl_surf_r(int Z, int A, double erg) const final;
  double cmpl_spin_r(int Z, int A, double erg) const final;
  
  double real_cent_a(int Z, int A, double erg) const final;
  double cmpl_cent_a(int Z, int A, double erg) const final;
  double cmpl_surf_a(int Z, int A, double erg) const final;
  double cmpl_spin_a(int Z, int A, double erg) const final;

  double real_cent_V(int Z, int A, double erg) const final;
  double cmpl_cent_V(int Z, int A, double erg) const final;
  double cmpl_surf_V(int Z, int A, double erg) const final;

  // CH89 does not have real surface or spin terms
  double real_surf_r(int Z, int A, double erg) const final;
  double real_spin_r(int Z, int A, double erg) const final;
  double real_surf_V(int Z, int A, double erg) const final;
  double real_spin_V(int Z, int A, double erg) const final;
  double real_surf_a(int Z, int A, double erg) const final;
  double real_spin_a(int Z, int A, double erg) const final;
  
  ChapelHill89<projectile>(json param_file);
  ChapelHill89<projectile>();
};

template <Proj p>
/// @brief constructs a ChapelHill89\<p\> with params refit w/ MCMC; from
/// Pruitt, C. D. et al, 
/// “Uncertainty-Quantified Phenomenological Optical Potentials 
/// for Single-Nucleon Scattering”, 
/// LLNL release number LLNL-JRNL-835671-DRAFT (to be published).
class CHUQ : public ChapelHill89<p> { public: CHUQ(); };

template <>
class CHUQ<Proj::proton> { public: CHUQ(); };

}

#endif
