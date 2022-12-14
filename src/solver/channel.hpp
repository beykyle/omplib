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

  /// @brief coupled AM states of the channel. This type is meant to be mutable,
  /// with l and J2 updated according to set_spin_up and set_spin_down
  struct AngularMomentum {
    /// @brief intrinsic spin of projectile
    const frac s;
    /// @brief total angular momentum of channel
    frac j;
    /// @brief orbital angular momentum [hbar]
    int l;
    /// @brief spatial parity of channel state
    Parity pi;


    Parity update_pi(int l) {
      return l % 2 == 0 ? Parity::even : Parity::odd;
    }

    /// @tparam n numerator of projection of projectile spin along j axis
    /// @tparam d denominator of projection of projectile spin along j axis
    template<int n, int d = 1>
    void set_spin(int new_l){
      assert(new_l > 0);
      l = new_l;
      pi = update_pi(l);
      set_spin<n,d>();
    }
    
    /// @tparam n numerator of projection of projectile spin along j axis
    /// @tparam d denominator of projection of projectile spin along j axis
    template<int n, int d = 1>
    void set_spin(){
      j = l + frac(n,d);
    }
    
    AngularMomentum(frac s, frac j, int l)
      : s(s), j(j), l(l), pi(update_pi(l)) 
    {
      assert(l >= 0);
      assert(j >= 0);
      assert(s >= 0);
    }
    
    /// @brief initialize AM to S-Wave state
    AngularMomentum(frac s)
      : s(s), j(s), l(0), pi(update_pi(l)) 
    {
      assert(l >= 0);
      assert(j >= 0);
      assert(s >= 0);
    }
    
    /// @returns projection of spin along axis of orbital ang. mom.; e.g. L * S
    real spin_orbit() const {
      const auto J = cast<real>(j);
      const auto S = cast<real>(s);
      const auto L = static_cast<real>(l);
      return 0.5*(J*(J+1) - L*(L+1) - S*(S+1));
    }
  };
  
  /// @brief things that change w/ energy
  struct Energetics {
    /// @brief CMS energy [MeV]
    real erg_cms;
    /// @brief Bombarding kinetic energy in the lab frame [MeV]
    real erg_lab;
    /// @brief CMS momentum  w/ relativistic correction according to 
    /// Eq. 17 of Ingemarsson, A. 
    /// "Some notes on optical model calculations at medium energies." 
    /// Physica Scripta 9.3 (1974): 156.
    real k;
    /// @brief scattering system reduced mass [MeV]
    real reduced_mass;
    /// @brief hbar^2 * c^2/(2 * reduced mass  * radius) [Mev fm]
    real h2ma;
    /// @brief Sommerfield parameter: sgn(Z_t * Z_p)/(a_B * k) [dimensionless]
    /// ab is the Bohr radius, Z_t and Z_p are the target and product charges
    real sommerfield_param;

    Energetics(real erg_cms, const Channel& ch)
    : erg_cms(erg_cms - ch.threshold) 
    , erg_lab(erg_cms 
        * (ch.proj_mass * constants::MeV_per_amu + ch.targ_mass * constants::MeV_per_amu) 
        / (ch.targ_mass * constants::MeV_per_amu))
    , k([this,&ch]() ->real { 
        
        using constants::MeV_per_amu;
        using constants::hbar;
        using constants::c;
        
        // relativistic-corrected wavenumber [1/fm]
        // Eq. 17 of Ingemarsson, 1974
        const auto m1 = ch.proj_mass * MeV_per_amu;
        const auto m2 = ch.targ_mass * MeV_per_amu;
        return m2 * sqrt(erg_lab * (erg_lab + 2 * m1 )) 
                  / sqrt( (m1 + m2)*(m1 + m2)  + 2 * m2 * erg_lab)
                  / (hbar * c);
        }())
    , reduced_mass([&erg_cms,&ch,this]() ->real { 
        
        using constants::hbar;
        using constants::c;
        using constants::MeV_per_amu;
        
        // relativistic-corrected reduced mass [MeV]
        // Eq. 21 of Ingemarsson, 1974
        const auto m1 = ch.proj_mass * MeV_per_amu;
        const auto Ep = m1 + erg_cms;

        return hbar * hbar * c * c * k * k * Ep / (Ep*Ep - m1*m1);

        }())
    , h2ma([&ch,this]() ->real { 
        
        using constants::hbar;
        using constants::c;

        return hbar * hbar * c * c / (2* reduced_mass * ch.radius); 

        }())
    , sommerfield_param([&ch,this]() ->real { 
        
        using constants::hbar;
        using constants::c;
        using constants::e_sqr;

        // Bohr radius [fm]
        const real ab = hbar * hbar * c * c / (reduced_mass * e_sqr * fabs(ch.Zz));
    
        // dimensionless
        return std::copysign(1., ch.Zz) / (ab * k);

        }())
     {}
  };

  /// @brief things that change w/ orbital angular momentum
  struct Asymptotics {
    /// @brief value of the outgoing asymptotic wavefunction at 
    /// the channel radius
    cmpl wvfxn_out;
    /// @brief value of the incoming asymptotic wavefunction at 
    /// the channel radius
    cmpl wvfxn_in;
    /// @brief value of the derivative of the outgoing asymptotic 
    /// wavefunction at the channel radius
    cmpl wvfxn_deriv_out;
    /// @brief value of the derivative of the incoming asymptotic 
    /// wavefunction at the channel radius
    cmpl wvfxn_deriv_in;

    Asymptotics(const AngularMomentum& am, real k, const Channel& ch) 
      : wvfxn_out( ch.radius * h_out{am.l}(k * ch.radius) )
      , wvfxn_in(  ch.radius * h_in{am.l} (k*ch.radius)   )
      , wvfxn_deriv_out( ReducedDeriv(h_out{am.l})(k*ch.radius) )
      , wvfxn_deriv_in(  ReducedDeriv(h_in{am.l}) (k*ch.radius) )
    {
      //TODO use confluent hypergeometrics for asymptotics
      // for now, assume neutral projectile and use Hankel fxns
      // if asym wavefunction = r*h_n(kr)
      // derivative = d/dr (r * h_n(kr)) = h_n(kr) + kr d/d(kr) h_n(kr)
      assert(ch.Zz == 0);
    }
  };

  /// @brief channel energy threshold [MeV]
  real threshold;
  /// @brief channel radius [fm]
  real radius;
  /// @brief target mass [amu]
  real targ_mass;
  /// @brief projectile mass [amu]
  real proj_mass;
  /// @brief product of proton # of target and projectile
  real Zz;
  /// @brief projectile spin
  frac s;
  

  Energetics set_erg_cms(real erg_cms) const {
    return Energetics(erg_cms, *this);
  }

  Energetics set_erg_lab(real erg_lab) const {
    using constants::MeV_per_amu;
    const auto m1 = proj_mass * MeV_per_amu;
    const auto m2 = targ_mass * MeV_per_amu;

    const auto erg_cms  = erg_lab * m2 / (m1 + m2); 
    return set_erg_cms(erg_cms);
  }
  
  Asymptotics set_angular_momentum(AngularMomentum am, real k) const {
    return Asymptotics(am, k, *this);
  }
  
  /// @param threshold [Mev]
  /// @param radius [fm]
  /// @param proj_mass [amu] 
  /// @param Zp proton number of projectile 
  /// @param s spin of projectile
  /// @param targ_mass [amu]
  /// @param Zt proton number of target 
  Channel(real threshold, real radius, 
          real proj_mass, int Zp, frac s,
          real targ_mass, int Zt)
    : threshold(threshold)
    , radius(radius)
    , targ_mass(targ_mass)
    , proj_mass(proj_mass)
    , Zz(Zp*Zt)
    , s(s) 
  {}
};

/// @brief Problem data for evaluating potentials and solving 
/// Schrodinger's eqn. Most potentials don't need all this, 
/// but some do.
struct PData {
  Channel ch;
  Channel::Energetics e;
  Channel::AngularMomentum am;
  
  PData(Channel ch, Channel::Energetics e, Channel::AngularMomentum am)
    : ch(ch), e(e), am(am) {}
};

};

#endif
