#ifndef SCATTER_HEADER
#define SCATTER_HEADER

#include <vector>
#include <array>

#include "util/constants.hpp"
#include "util/config.hpp"

#include "solver/rmatrix.hpp"
#include "potential/potential.hpp"

namespace omplib {

/// @tparam Nucleon projectile
template<Proj proj>
/// @brief single channel scatter for incident nucleon on a nucleus
class NAScatter {
public:
  
  struct CrossSection {
    real tot;
    real rxn;
  };

  struct Amplitude {
    /// @brief spin-preserving amplitude
    std::array<cmpl, MAXL> A;
    /// @brief spin-flipping amplitude
    std::array<cmpl, MAXL> B;

    /// @brief conversion factor from amplitude to T-Matrix
    real ampl_to_T;
  };
  
  using DiffXS         = std::vector<real>;
  using AnalyzingPower = std::vector<real>;

  struct AngularData {
    DiffXS         dxdu;
    AnalyzingPower Ay;
    AngularData(unsigned int sz): dxdu(sz,0), Ay(sz,0) {};
  };

private:
  using Solver = RMatrixSolverSingleChannel<NBASIS>;
  /// @brief Dimension of spin 1/2 repr of SU(2) = (2*1/2 +1)
  constexpr static int S2 = 2;
  Isotope target;
  
  real ch_radius;
  real ch_threshold;


public:
  NAScatter(Isotope target, real ch_radius=15, real ch_threshold=0)
    : target(target)
    , ch_radius(ch_radius)
    , ch_threshold(ch_threshold) {};

  /// @tparam Potential callable (real r [fm], real rp [fm], Channel ch) -> cmpl [Mev]
  /// @tparam SideEffect function that process/stores the S-Matrix elements for each 
  /// partial wave
  template<class Potential, class SideEffect>
  /// @brief calculates the S-Matrix for each partial wave amplitude,
  /// and passes it into func for processing
  /// @returns The maxmium l value used before convergence 
  /// worst case time complexity O(MAXL * NBASIS^3 )
  int core_solver(
      Potential p, 
      const Channel& ch, 
      const Channel::Energetics& e, 
      SideEffect func,
      real rel_err_crit = REL_ERR_EPS) const 
  {
    // TODO for now just neutrons
    static_assert(charge<proj>() == 0);
  
    auto am = Channel::AngularMomentum{S2}; 
    
    // S-wave can only have spin up states   
    cmpl Sminus = 0;
    cmpl Splus  = Solver(ch, e, am, p).smatrix();
    func(am.l,Sminus,Splus);

    auto sm_old = Sminus;

    auto rel_err = [](auto l, auto r) { return (l-r)/r; };
    
    // higher OAM states must account for spin up and down
    for (am.l = 1; am.l < MAXL; ++am.l) {
      // twice the projections of the projectile spin
      // along the channel ang momentum axis
      constexpr int up   =  1; //  ms =  1/2 => 2ms =  1
      constexpr int down = -1; // -ms = -1/2 => 2ms = -1
      
      // spin down
      am.set_spin<down>(am.l);
      Sminus = Solver(ch, e, am, p).smatrix();
      
      // spin
      am.set_spin<up>(am.l);
      Splus = Solver(ch, e, am , p).smatrix();
      
      func(am.l,Sminus,Splus);
      
      if (
          rel_err(Sminus.real(), sm_old.real()) < rel_err_crit and 
          rel_err(Sminus.imag(), sm_old.imag()) < rel_err_crit 
         ) break;

      sm_old = Sminus;

    }

    return am.l;
  } 

  /// @tparam Potential callable (real r [fm], real rp [fm], Channel ch) -> cmpl [Mev]
  template<class Potential>
  /// @brief  Calculate the partial wave amplitudes (porportional to the T-Matrix)
  /// for spin-flip and spin-preserving processes. Eqns. 8,9 in Ingermarsson, 1947.
  Amplitude amplitudes(real erg_cms, Potential p) const {
    using constants::pi;
    using constants::hbar;
    
    const auto ch = Channel(0, ch_radius, mass<proj>(), charge<proj>(), target.mass , target.Z);
    const auto e  = ch.set_erg_cms(erg_cms);

    Amplitude ampl;
    auto& [A, B, conv] = ampl;
    // set conversion factor from ampl to T-Matrix
    conv = -4 * pi*pi * e.reduced_mass / (hbar*hbar);
    
    const auto f = constants::i / (2 * e.k);

    // tabulte spin-flip and non-spin-flip amplitudes with solver
    core_solver(p, ch, e, 
        [&ampl, f](int l, cmpl Sminus, cmpl Splus) {
          ampl.A[l] = ( (real)(2 *l +1) - (real)(l+1) * Splus - (real)l * Sminus)*f;
          ampl.B[l] = ( Sminus - Splus)*f;
        } 
    );
    
    return ampl;
  } 
  
  /// @tparam Potential callable (real r [fm], real rp [fm], Channel ch) -> cmpl [Mev]
  template<class Potential>
  /// @brief  Calculate the total and reaction cross section for the channel; 
  /// Eqns 46-47 in Pruitt, 2021
  CrossSection xs(real erg_cms, Potential p) const {
  
    const auto ch = Channel(0, ch_radius, mass<proj>(), charge<proj>(), target.mass , target.Z);
    const auto e  = ch.set_erg_cms(erg_cms);

    CrossSection xs{0,0};

    // tally up total and rxn cross sections with solver
    core_solver(p, ch, e, 
        [&xs](int l, cmpl Sminus, cmpl Splus) {
          xs.rxn += (l+1) * (1 - norm(Splus) ) + l * (1 - norm(Sminus) );
          xs.tot += (l+1) * (1 - Splus.real()) + l * (1 - Sminus.real());
        }
    );

    return xs;
  }

  /// @param ampl partial-wave amplitudes
  /// @param mu_grid Grid of mu = cos(theta) on [-1,1] on which to calculate 
  /// analyzing powers and differential cross section
  AngularData angular(const Amplitude& ampl, std::vector<real> mu_grid) const {
    assert(mu_grid.front() >= -1);
    assert(mu_grid.back()  <=  1);

    AngularData data(mu_grid.size());
    auto& [dxdu, Ay] = data;

    for (int i = 0; i < mu_grid.size(); ++i) {
      cmpl a{0,0};
      cmpl b{0,0};
      for (int l = 0; l < MAXL; ++l) {
         a += ampl.A[l] * std::legendre(l,mu_grid[i]);
         b += ampl.B[l] * std::assoc_legendre(l,1,mu_grid[i]);
      }
      dxdu[i] = a * conj(a) + b * conj(b);
      Ay[i]   = a * conj(b) + b * conj(a);
      Ay[i]  /= dxdu[i];
    }
    return data;
  }

};

}
#endif 
