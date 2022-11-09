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


  struct Amplitude {
    /// @brief spin-preserving amplitude
    std::array<cmpl, MAXL> A;
    /// @brief spin-flipping amplitude
    std::array<cmpl, MAXL> B;

    /// @brief conversion factor from amplitude to T-Matrix
    real ampl_to_T;
  };


  struct CrossSection {
    real tot;
    real rxn;
  };
  
  struct Solution {
    CrossSection xs;
    Amplitude T;
  };
  
  /// @tparam Potential callable (real r [fm], real rp [fm], Channel ch) -> cmpl [Mev]
  /// @tparam SideEffect
  template<class Potential, class SideEffect>
  /// @brief calculates the S-Matrix for each partial wave amplitude,
  /// and passes it into func for processing
  /// @returns The maxmium l value used before convergence 
  /// worst case time complexity O(MAXL * NBASIS^3 )
  int core_solver(
      Potential p, const Channel& ch, const Channel::Energetics& e, SideEffect func) const 
  {
    // TODO for now just neutrons
    static_assert(charge<proj>() == 0);
  
    auto am = Channel::AngularMomentum{0, S2, S2}; 
    
    // S-wave can only have spin up states   
    cmpl Sminus = 0;
    cmpl Splus  = Solver(ch, e, am, p).smatrix();
    func(am.l,Sminus,Splus);
    
    // higher OAM states must account for spin up and down
    for (am.l = 1; am.l < MAXL; ++am.l) {
      
      am.set_spin_down(am.l);
      Sminus = Solver(ch, e, am, p).smatrix();
      
      am.set_spin_up(am.l);
      Splus = Solver(ch, e, am , p).smatrix();
      

      func(am.l,Sminus,Splus);
    }

    return am.l;
  } 

  /// @tparam Potential callable (real r [fm], real rp [fm], Channel ch) -> cmpl [Mev]
  /// @brief  Calculate the partial wave amplitudes (porportional to the T-Matrix)
  /// for spin-flip and spin-preserving processes. Eqns. 8,9 in Ingermarsson, 1947.
  /// memory complexity O(NBASIS^2 + MAXL)
  /// time complexity O(MAXL * NBASIS^3 )
  template<class Potential>
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
  /// @brief  Calculate the total and reaction cross section for the channel; 
  /// Eqns 46-47 in Pruitt, 2021
  /// memory complexity O(NBASIS^2)
  /// time complexity O(MAXL * NBASIS^3 )
  template<class Potential>
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
};

}
#endif 
