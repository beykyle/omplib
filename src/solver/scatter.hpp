#ifndef SCATTER_HEADER
#define SCATTER_HEADER

#include "util/constants.hpp"
#include "util/config.hpp"

#include "solver/rmatrix.hpp"
#include "potential/potential.hpp"

#include <vector>
#include <array>

namespace omplib {

struct Isotope {
  int A;
  int Z;
  //TODO load a mass table
  real mass;
};

template<Proj p>
constexpr real mass() {
  if constexpr (p == Proj::neutron) return constants::n_mass_amu;
  if constexpr (p == Proj::proton ) return constants::p_mass_amu;
}

template<Proj p>
constexpr int charge() {
  if constexpr (p == Proj::neutron) return 0;
  if constexpr (p == Proj::proton ) return 1;
}

template<Proj proj>
class NAScatter {
private:
  using Solver = RMatrixSolverSingleChannel<NBASIS>;
  /// @brief Dimension of spin 1/2 repr of SU(2) = (2*1/2 +1)
  constexpr static int S2 = 2;
  Isotope target;
  
  real    ch_radius;
  real    ch_threshold;


public:
  NAScatter(Isotope target, real ch_radius=15, real ch_threshold=0)
    : target(target)
    , ch_radius(ch_radius)
    , ch_threshold(ch_threshold) {};

  struct Matrix {
    std::array<cmpl, MAXL> A;
    std::array<cmpl, MAXL> B;
  };

  /// @tparam Potential callable (real r [fm], real rp [fm], Channel ch) -> cmpl [Mev]
  template<class Potential>
  Matrix solve(real erg_cms, Potential p) const {

    // TODO for now just neutrons
    static_assert(charge<proj>() == 0);
  
    Matrix soln;
    auto& [A, B] = soln;

    int l   = 0;
    auto ch = Channel(0, erg_cms, ch_radius, 
                    mass<proj>(), charge<proj>(),
                    target.mass , target.Z,  
                    0, 0, S2);

    // S-wave can only have spin up states   
    cmpl Sminus = 0;
    cmpl Splus  = Solver(ch, p).smatrix();
    A[0] =  Splus;
    B[0] = -Splus;
    
    // higher OAM states must account for spin up and down
    for (l = 1; l < MAXL; ++l) {
      ch.reset_l(l); // re-calculate asymptotics
      
      // spin up
      ch.J2 = 2*l + S2; // 2*(l+1/2) +1 = 2l + (2*1/2+1)
      Splus = Solver(ch, p).smatrix();
      
      // spin down
      ch.J2 = 2*l; // 2*(l-1/2)  +1 = 2l
      Sminus = Solver(ch, p).smatrix();

      A[l] = (real)(2*l+1) - (real)(l+1) * Splus - (real)l * Sminus;
      B[l] = Sminus - Splus;
    } 
    
    return soln; 
  } 
};

}
#endif 
