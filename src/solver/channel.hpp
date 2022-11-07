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
  /// @brief scattering system reduced mass [amu]
  real reduced_mass;
  
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
  int    J2;
  /// @brief spatial parity of channel state
  Parity pi;

  Channel(real threshold, real energy, real radius, 
          real proj_mass, real targ_mass,
          int Zt, int Zp, int l, int J2, Parity pi )
    : threshold(threshold)
    , energy(energy)
    , radius(radius)
    , reduced_mass( targ_mass * proj_mass / (targ_mass + proj_mass) )
    , l(l)
    , J2(J2)
    , pi(pi)
  {
    using constants::hbar;
    using constants::c;
    using constants::e_sqr;
    
    const auto m1   = proj_mass * constants::MeV_per_amu / c / c;
    const auto m2   = targ_mass * constants::MeV_per_amu / c / c;

    // non-relativistic
    k = sqrt(2 * reduced_mass * constants::MeV_per_amu * energy) / (hbar * c);
    
    // relativistic correction
    const auto Tlab = energy * (m1 + m2) / m2; // Eq. 4 Inermarsson, 1974
    // Eq. 17 of Ingemarsson, 1974
    k = m2 * c * c 
      * sqrt(Tlab * (Tlab + 2 * m1 * c * c)) 
      / sqrt((m1 + m2)*(m1 + m2) * c * c * c * c + 2 * m2 * c * c * Tlab)
      / (hbar * c);

    const real Zz = Zt * Zp;
    const real ab = hbar * hbar / (reduced_mass * e_sqr * fabs(Zz));

    sommerfield_param = std::copysign(1, Zz) / (ab * k);

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
