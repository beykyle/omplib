#ifndef PARAMS_HPP
#define PARAMS_HPP

#include "util/constants.hpp"

#include "nlohmann/json.hpp"
using nlohmann::json;

namespace omplib {

/// @brief incident particle type
enum class Proj : bool {
  proton,
  neutron,
};
 
/// @tparam  incident neutron or proton 
template<Proj p>
/// @brief  Pure abstract base class for OM potential parameters
struct OMParams {
  constexpr static Proj projectile = p;
  
  // Woods-Saxon term radii 
  virtual double real_cent_r(int Z, int A, double erg) const = 0;
  virtual double cmpl_cent_r(int Z, int A, double erg) const = 0;
  virtual double real_surf_r(int Z, int A, double erg) const = 0;
  virtual double cmpl_surf_r(int Z, int A, double erg) const = 0;
  virtual double real_spin_r(int Z, int A, double erg) const = 0;
  virtual double cmpl_spin_r(int Z, int A, double erg) const = 0;
  
  // Woods-Saxon term diffusivity
  virtual double real_cent_a(int Z, int A, double erg) const = 0;
  virtual double cmpl_cent_a(int Z, int A, double erg) const = 0;
  virtual double real_surf_a(int Z, int A, double erg) const = 0;
  virtual double cmpl_surf_a(int Z, int A, double erg) const = 0;
  virtual double real_spin_a(int Z, int A, double erg) const = 0;
  virtual double cmpl_spin_a(int Z, int A, double erg) const = 0;

  // Woods-Saxon term depth
  virtual double real_cent_V(int Z, int A, double erg) const = 0;
  virtual double cmpl_cent_V(int Z, int A, double erg) const = 0;
  virtual double real_surf_V(int Z, int A, double erg) const = 0;
  virtual double cmpl_surf_V(int Z, int A, double erg) const = 0;
  virtual double real_spin_V(int Z, int A, double erg) const = 0;
  virtual double cmpl_spin_V(int Z, int A, double erg) const = 0;

protected:

  /// @returns -(N-Z)/A  when p == neutron and +(N-Z)/A when p == proton
  static double asym(int Z, int A);

};

template<>
struct OMParams<Proj::proton> : public OMParams<Proj::neutron> {
  virtual double real_coul_r(int Z, int A, double erg) const = 0;

  /// @brief Coulomb potential w/in uniformly charge sphere or radius R
  /// v(r) = q^2/(2*R)(3 - r/R)
  virtual double real_coul_V_outer(int Z, int A, double erg) const {
    return static_cast<double>(Z) * e_sqr;
  }
  /// @brief Coulomb potential outside uniformly charge sphere or radius R
  /// v(r) = q^2/r
  virtual double real_coul_V_inner(int Z, int A, double erg) const { 
    const double a = static_cast<double>(A);
    const double z = static_cast<double>(Z);
    const double RC = real_coul_r(Z,A,erg) * pow(a,1./3.);
    return z * e_sqr / (2 * RC)  ;
  };
};

template<Proj proj>
double OMParams<proj>::asym(int Z, int A) {
  const double n = static_cast<double>(A - Z);
  const double a = static_cast<double>(A);
  const double z = static_cast<double>(Z);
  const double alpha = (n - z)/a; 
  if constexpr      (proj == Proj::neutron) return  alpha;
  else if constexpr (proj == Proj::proton ) return -alpha;
}

}

#endif
