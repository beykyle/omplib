#ifndef CONFIG_HEADER
#define CONFIG_HEADER

#include "util/types.hpp"

namespace omplib {
  /// @brief  The number of basis functions to use in R-Matrix calculations
  constexpr int NBASIS  = 20;
  /// @brief The maximum number of orbital angular momentum values to consider
  /// for a single channel
  constexpr int MAXL    = 20;
  /// @brief Relative error between iterations indicating convergence
  constexpr real REL_ERR_EPS  = 1.0e-8;
}

#endif 
