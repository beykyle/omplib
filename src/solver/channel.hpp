#ifndef CHANNEL_HEADER
#define CHANNEL_HEADER

#include <cmath>
#include <complex>
#include <cassert>

#include "util/asymptotics.hpp"
#include "util/constants.hpp"
#include "util/types.hpp"

namespace omplib {

enum class Parity : bool {
  odd,
  even
};

struct Channel {
  /// @brief channel energy threshold [MeV]
  real threshold;
  /// @brief CMS energy [MeV]
  real energy;
  /// @brief channel radius [fm]
  real radius;
  /// @brief scattering system reduced mass [MeV]
  real reduced_mass;
  /// @brief hbar^2 * c^2/(2 * reduced mass  * radius) [Mev fm]
  real h2ma;
  /// @brief Bombarding kinetic energy in the lab frame [MeV]
  real Tlab;
  
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
  /// @param energy (cms) [Mev]
  /// @param radius [fm]
  /// @param proj_mass [amu] 
  /// @param targ_mass [amu]
  /// @param Zt proton number of target 
  /// @param Zp proton number of projectile 
  /// @param l orbital angular momentum [hbar]
  /// @param J2 total angular momentum  [hbar]
  /// @param pi parity [hbar]
  Channel(real threshold, real erg_cms, real radius, 
          real proj_mass,  int Zp,
          real targ_mass,  int Zt, 
          int l, int J2, int S2)
    : threshold(threshold)
    , energy(erg_cms)
    , radius(radius)
    , l(l)
    , J2(J2)
    , S2(S2)
    , pi( l % 2 == 0 ? Parity::even : Parity::odd)
  {
    using constants::hbar;
    using constants::c;
    using constants::e_sqr;
    using constants::MeV_per_amu;
    
    const auto energy = erg_cms - threshold;
    const auto m1     = proj_mass * MeV_per_amu;
    const auto m2     = targ_mass * MeV_per_amu;

    // for debugging relativistic corrections
    // non-relativistic wave number [1/fm]
    k = sqrt(2 * reduced_mass * constants::MeV_per_amu * energy) / (hbar * c);
    // non-relativistic reduced mass [MeV]
    reduced_mass = targ_mass * proj_mass / (targ_mass + proj_mass) 
                 * MeV_per_amu;
    
    // [Mev]
    Tlab = energy * (m1 + m2) / m2; // Eq. 4 Inermarsson, 1974
    
    // relativistic-corrected wavenumber [1/fm]
    // Eq. 17 of Ingemarsson, 1974
    k = m2 * sqrt(Tlab * ( Tlab + 2 * m1 )) 
      / sqrt( (m1 + m2)*(m1 + m2)  + 2 * m2 * Tlab)
      / (hbar * c);


    // relativistic-corrected reduced mass [MeV]
    // Eq. 21 of Ingemarsson, 1974
    const auto Ep = m1 + energy;
    reduced_mass  = hbar * hbar * c * c * k * k * Ep / (Ep*Ep - m1*m1);
    
    // [MeV fm]
    h2ma = hbar * hbar * c * c / (2* reduced_mass * radius); 

    // product of charges 
    const real Zz = Zt * Zp;
    
    // Bohr radius [fm]
    const real ab = hbar * hbar * c * c / (reduced_mass * e_sqr * fabs(Zz));
    
    // dimensionless
    sommerfield_param = std::copysign(1., Zz) / (ab * k);

    //TODO use confluent hypergeometrics for asymptotics
    // for now, assume neutral projectile and use Hankel fxns
    // if asym wavefunction = r*h_n(kr)
    // derivative = d/dr (r * h_n(kr)) = h_n(kr) + kr d/d(kr) h_n(kr)
    assert(Zp == 0);

    const auto kr = k*radius;

    asymptotic_wvfxn_out       = radius * h_out{l}(kr);
    asymptotic_wvfxn_deriv_out = (1.+l) * h_out{l}(kr) -  kr * h_out{l+1}(kr);
    
    asymptotic_wvfxn_in        = radius * h_in{l}(k*radius);
    asymptotic_wvfxn_deriv_in  = (1.+l) * h_in{l}(kr) -  kr * h_in{l+1}(kr);
  }
};

};

#endif
