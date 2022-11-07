#ifndef CONSTANTS_HEADER
#define CONSTANTS_HEADER

#include "util/types.hpp"

#include <complex>

namespace omplib {
namespace constants {
  using namespace std::complex_literals;
  ///@brief square of elementary charge [ Mev fm ]
  constexpr static real e_sqr = 1.4399764;
  ///@brief pion mass [ Mev ]
  constexpr static real m_pi  = 139.57061;
  /// @brief hbar [Mev s]
  constexpr static real hbar = 6.58212196e-22;
  /// @brief Coulomb's constant 
  constexpr static real q_coulomb = 1.60217733e-19; 
  /// @brief free-space permittivity [MeV fm / c^2 ]
  constexpr static real eps_permitiviyty = 5.60958617e+37;
  /// @brief zoom fast light speed [fm/s]
  constexpr static real c = 2.99792458e+23;
  /// @brief sqrt(-1)
  constexpr static cmpl i = 1i;
  
  // dimensionless
  /// @brief hbar^2/(m_pion * c)^2 
  constexpr static real csp = 2.04553; 
  ///@brief bruh, it's pi
  constexpr static real pi    = 3.1415926536;
  /// @brief fine structure
  constexpr static real alpha = 7.2973525376e-03;

  // [amu]
  constexpr static real n_mass_amu = 1.00866491578;
  constexpr static real e_mass_amu = 5.48579911e-4;
  constexpr static real p_mass_amu = 1.00727646688;

  constexpr static real ev_per_amu  = 9.31494013e+8; 
  constexpr static real MeV_per_amu = 931.4943335;
};
};


#endif
