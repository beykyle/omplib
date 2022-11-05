#ifndef CHANNEL_HEADER
#define CHANNEL_HEADER

#include <cmath>
#include <complex>
#include <cassert>

#include "util/asymptotics.hpp"
#include "util/constants.hpp"

namespace omplib {

enum class : bool Parity {
  odd,
  even
};

struct Channel {
  /// @brief channel energy threshold [MeV]
  double threshold;
  /// @brief CMS energy [MeV]
  double energy;
  /// @brief channel radius [fm]
  double radius;
  /// @brief scattering system reduced mass [MeV]
  double reduced_mass;
  
  /// @brief CMS momentum sqrt( 2 * E * reduced_mass * c^2)/hbar 
  double k;

  /// @brief Sommerfield parameter: sgn(Z_t * Z_p)/(a_B * k) [dimensionless]
  /// ab is the Bohr radius, Z_t and Z_p are the target and product charges
  double sommerfield_param;
  
  /// @brief value of the outgoing asymptotic wavefunction at 
  /// the channel radius
  std::complex<double> asymptotic_wvfxn_out;
  /// @brief value of the derivative of the outgoing asymptotic 
  /// wavefunction at the channel radius
  std::complex<double> asymptotic_wvfxn_deriv_out;
  /// @brief value of the incoming asymptotic wavefunction at 
  /// the channel radius
  std::complex<double> asymptotic_wvfxn_in;
  /// @brief value of the derivative of the incoming asymptotic 
  /// wavefunction at the channel radius
  std::complex<double> asymptotic_wvfxn_deriv_in;
  
  /// @brief orbital angular momentum of scattering system [hbar]
  int l;
  /// @brief (2J+1), where J is the total angular 
  /// momentum of the channel [hbar]
  int    J2;
  /// @brief spatial parity of channel state
  Parity pi;

  Channel(double threshold, double energy, double radius, double reduced_mass, 
          int Zt, int Zp, int l, int J2, Parity pi )
    : threshold(threshold)
    , energy(energy)
    , radius(radius)
    , reduced_mass(reduced_mass)
    , l(l)
    , J2(J2)
    , pi(pi)
  {
    using constants::hbar;
    using constants::c;
    using constants::e_sqr;

    k = sqrt(2 * reduced_mass * energy) * c / hbar;
    
    const double Zz = Zt * Zp;
    const double ab = hbar * hbar / (reduced_mass * e_sqr * abs(Zz));

    sommerfield_param = std::copysign(1, Zz) / (ab * k);

    //TODO use confluent hypergeometrics for asymptotics
    // for now, assume neutral projectile and use Hankel fxns
    assert(Zp == 0);

    asymptotic_wvfxn_out       = h_out{l}(kr);
    asymptotic_wvfxn_deriv_out = hp_out{l}(kr);
    
    asymptotic_wvfxn_in        = h_in{l}(kr);
    asymptotic_wvfxn_deriv_in  = hp_in{l}(kr);
  }
  
};

};

#endif
