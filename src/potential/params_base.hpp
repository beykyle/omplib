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
  double asym(int Z, int A) const;

};

template<>
struct OMParams<Proj::proton> {
  virtual double real_coul_r(int Z, int A, double erg) const = 0;
  virtual double real_coul_V(int Z, int A, double erg) const = 0;
};

template<Proj proj>
double OMParams<proj>::asym(int Z, int A) const {
  const double n = static_cast<double>(A - Z);
  const double a = static_cast<double>(A);
  const double z = static_cast<double>(Z);
  return (n - z)/a;
}


}

#endif
