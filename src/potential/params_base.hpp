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
 
/// @tparam p incident neutron or proton 
///
template<Proj p>
/// @brief  Pure abstract base class for OM potential parameters, where the potential 
/// takes the common global phenomenological form:
/// V(r,E) + V_s(r,E) + V_so(r,E) * L.S + i*W_v(r,E) + i*W_s(r,E) + i W_so(r,E) * L.S
///
/// Here, V,W_v refer to real and an imaginary volume terms, which are typically 
/// negative/attractive, and take Woods-Saxon form factors f(r,R,a) = 1/(1+exp((r-R)/a))
/// V(r,erg) = real_cent_V(Z,A,erg) * f(r,real_cent_R(Z,A,erg), real_cent_a(Z,A,erg))
/// W(r,erg) = cmpl_cent_V(Z,A,erg) * f(r,cmpl_cent_R(Z,A,erg), cmpl_cent_a(Z,A,erg))
///
/// V_s, and W_s refer to surface peaked potentials, typically positive/repulsive, and take
/// derivative of Woods-Saxon form factors.
/// V_s(r,erg) = real_surf_V(Z,A,erg) * d/dr f(r,real_surf_R(Z,A,erg), real_surf_a(Z,A,erg))
/// W_s(r,erg) = cmpl_surf_V(Z,A,erg) * d/dr f(r,cmpl_surf_R(Z,A,erg), cmpl_surf_a(Z,A,erg))
/// The typical factor of -4*a_s is absorbed into *_surf_V.
///
/// V_so, and W_so refer to spin orbit potentials, typically positive/repulsive, and take
/// derivative of Woods-Saxon form factors, modulated by 1/e.
/// V_so(r,erg) = real_spin_V(Z,A,erg) * 1/r d/dr f(r,real_spin_R(Z,A,erg), real_spin_a(Z,A,erg))
/// W_so(r,erg) = cmpl_spin_V(Z,A,erg) * 1/r  d/dr f(r,cmpl_spin_R(Z,A,erg), cmpl_spin_a(Z,A,erg))
///
/// The typical factor of inverse pion mass squared is absorbed into *_spin_V.
/// All depth terms have units of MeV, and r and a terms of fm.
///
/// See the function omplib::OMP::eval for how these parameters are used.
///
/// NOTE: radii are NOT reduced - they are the exact radii R that get plugged into f(r,R,a).
/// To get reduced radii, divide by A^(1/3)
struct OMParams {
  constexpr static Proj projectile = p;
  
  // Woods-Saxon term radii 
  // NOTE these are not reduced radii (r/a^(1/3))
  // they are the exact radii that get plugged into a Wood-Saxon term
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
struct OMParams<Proj::proton>  {
  virtual double real_coul_r(int Z, int A, double erg) const = 0;

  /// @brief Coulomb potential w/in uniformly charge sphere or radius R
  /// v(r) = q^2/(2*R)(3 - r/R)
  virtual double real_coul_V_outer(int Z, int A, double erg) const {
    return static_cast<double>(Z) * constants::e_sqr;
  }
  /// @brief Coulomb potential outside uniformly charge sphere or radius R
  /// v(r) = q^2/r
  virtual double real_coul_V_inner(int Z, int A, double erg) const { 
    const double a = static_cast<double>(A);
    const double z = static_cast<double>(Z);
    const double RC = real_coul_r(Z,A,erg) * pow(a,1./3.);
    return z * constants::e_sqr / (2 * RC)  ;
  };
  
  /// @returns +(N-Z)/A 
  static double asym(int Z, int A);
};

}

#endif
