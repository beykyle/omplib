#ifndef CHANNEL_HEADER
#define CHANNEL_HEADER

#include <cmath>
#include <complex>
#include <cassert>
#include <iterator>

#include "util/asymptotics.hpp"
#include "util/constants.hpp"
#include "util/types.hpp"

namespace omplib {

enum class Parity : bool {
  odd,
  even
};

/// @brief  all information relevant to a specific initial/final state 
struct Channel {
  /// @brief channel energy threshold [MeV]
  real threshold;
  /// @brief CMS energy [MeV]
  real erg_cms;
  /// @brief Bombarding kinetic energy in the lab frame [MeV]
  real erg_lab;
  /// @brief channel radius [fm]
  real radius;
  
  /// @brief scattering system reduced mass [MeV]
  real reduced_mass;
  /// @brief target mass [amu]
  real target_mass;
  /// @brief projectile mass [amu]
  real proj_mass;
  /// @brief hbar^2 * c^2/(2 * reduced mass  * radius) [Mev fm]
  real h2ma;

  /// @brief product of proton # of target and projectile
  real Zz;
  
  /// @brief CMS momentum  w/ relativistic correction according to 
  /// Eq. 17 of Ingemarsson, A. 
  /// "Some notes on optical model calculations at medium energies." 
  /// Physica Scripta 9.3 (1974): 156.
  real k;

  /// @brief Sommerfield parameter: sgn(Z_t * Z_p)/(a_B * k) [dimensionless]
  /// ab is the Bohr radius, Z_t and Z_p are the target and product charges
  real sommerfield_param;
  
  /// @brief value of the outgoing asymptotic wavefunction at 
  /// the channel radius
  cmpl asymptotic_wvfxn_out;
  /// @brief value of the derivative of the outgoing asymptotic 
  /// wavefunction at the channel radius
  cmpl asymptotic_wvfxn_deriv_out;
  /// @brief value of the incoming asymptotic wavefunction at 
  /// the channel radius
  cmpl asymptotic_wvfxn_in;
  /// @brief value of the derivative of the incoming asymptotic 
  /// wavefunction at the channel radius
  cmpl asymptotic_wvfxn_deriv_in;
  
  /// @brief orbital angular momentum of scattering system [hbar]
  int l;
  /// @brief (2J+1), where J is the total angular 
  /// momentum of the channel [hbar]
  int J2;
  /// @brief (2S+1), where S is the intrinsic angular momentum of
  /// the proprojectile
  int S2;
  /// @brief spatial parity of channel state
  Parity pi;
  
  /// @returns projection of spin along axis of orbital ang. mom.; e.g. L * S
  real spin_orbit() const {
    const real L = l;
    const real S = (S2 - 1.)/2.;
    const real J = (J2 - 1.)/2.;

    return 0.5*(J*(J+1) - L*(L+1) - S*(S+1));
  }

  /// @param threshold [Mev]
  /// @param erg_cms [Mev]
  /// @param radius [fm]
  /// @param proj_mass [amu] 
  /// @param targ_mass [amu]
  /// @param Zt proton number of target 
  /// @param Zp proton number of projectile 
  /// @param l orbital angular momentum [hbar]
  /// @param J2 = 2*j+1, where j is total angular momentum  [hbar]
  /// @param S2 = 2*s+1, where s is projectile spin angular momentum [hbar]
  Channel(real threshold, real erg_cms, real radius, 
          real proj_mass, int Zp,
          real targ_mass, int Zt, 
          int l, int J2, int S2)
    : threshold(threshold)
    , erg_cms(erg_cms)
    , radius(radius)
    , target_mass(targ_mass)
    , proj_mass(proj_mass)
    , Zz(Zp*Zt)
    , l(l)
    , J2(J2)
    , S2(S2)
    , pi( l % 2 == 0 ? Parity::even : Parity::odd)
  {
      reset_erg_cms(erg_cms);
      reset_l(l);
  }

  /// @brief resets orbital AM, and re-calculates asymptotics
  void reset_l(int ln) {
    
    //TODO use confluent hypergeometrics for asymptotics
    // for now, assume neutral projectile and use Hankel fxns
    // if asym wavefunction = r*h_n(kr)
    // derivative = d/dr (r * h_n(kr)) = h_n(kr) + kr d/d(kr) h_n(kr)
    assert(Zz == 0);
    
    const auto kr = k*radius;
    l = ln;
    asymptotic_wvfxn_out       = radius * h_out{l}(kr);
    asymptotic_wvfxn_deriv_out = (1.+l) * h_out{l}(kr) -  kr * h_out{l+1}(kr);
    
    asymptotic_wvfxn_in        = radius * h_in{l}(k*radius);
    asymptotic_wvfxn_deriv_in  = (1.+l) * h_in{l}(kr) -  kr * h_in{l+1}(kr);
  }
  
  /// @param erg_cmsn new cms frame energy [Mev]
  /// @brief Holds target, projectile, and other params constant
  /// Recalculates erg_lab, k, reduced_mass, h2ma, sommerfield_param
  void reset_erg_cms(real erg_cmsn) {
    
    using constants::hbar;
    using constants::c;
    using constants::e_sqr;
    using constants::MeV_per_amu;
    
    const auto erg_cms = erg_cmsn - threshold;
    const auto m1      = proj_mass * MeV_per_amu;
    const auto m2      = target_mass * MeV_per_amu;
    
    erg_lab = erg_cms * (m1 + m2) / m2; // Eq. 4 Inermarsson, 1974

    // relativistic-corrected wavenumber [1/fm]
    // Eq. 17 of Ingemarsson, 1974
    k = m2 * sqrt(erg_lab * ( erg_lab + 2 * m1 )) 
      / sqrt( (m1 + m2)*(m1 + m2)  + 2 * m2 * erg_lab)
      / (hbar * c);

    // relativistic-corrected reduced mass [MeV]
    // Eq. 21 of Ingemarsson, 1974
    const auto Ep = m1 + erg_cms;
    reduced_mass  = hbar * hbar * c * c * k * k * Ep / (Ep*Ep - m1*m1);
    
    // [MeV fm]
    h2ma = hbar * hbar * c * c / (2* reduced_mass * radius); 
    
    // Bohr radius [fm]
    const real ab =  bohr_radius();
    
    // dimensionless
    sommerfield_param = std::copysign(1., Zz) / (ab * k);
  }
  
  /// @param erg_labn new lab frame energy [Mev]
  /// @brief Holds target, projectile, and other params constant
  /// Recalculates erg_lab, k, reduced_mass, h2ma, sommerfield_param
  void reset_erg_lab(real erg_labn) {
    using constants::MeV_per_amu;
    const auto m1 = proj_mass * MeV_per_amu;
    const auto m2 = target_mass * MeV_per_amu;

    erg_cms  = erg_labn * m2 / (m1 + m2); 
    reset_erg_cms(erg_cms);
  }
  
  /// @param threshn new channel threshold [Mev]
  /// @brief Holds target, projectile, and other params constant
  /// Recalculates erg_lab, k, reduced_mass, h2ma, sommerfield_param
  void reset_threshold(real threshn) {
    threshold = threshn;
    reset_erg_cms(erg_cms);
  }

  /// @returns the Bohr radius in [fm]
  real bohr_radius() const {
    using constants::hbar, constants::c, constants::e_sqr;
    return hbar * hbar * c * c / (reduced_mass * e_sqr * fabs(Zz));
  }
};

};

#endif
