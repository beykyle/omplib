#ifndef RMATRIX_HEADER
#define RMATRIX_HEADER

#include "util/constants.hpp"
#include "util/asymptotics.hpp"
#include "util/lagrange_legendre.hpp"

#include "solver/channel.hpp"

#include <Eigen/LU>

namespace omplib {

/// @tparam the number of Lagrange-Legendre basis functions for integration
template<unsigned int N>
/// @brief Solves the scattering system in a single, uncoupled channel, using the R-Matrix 
/// method, following the general procedure of:
/// Descouvemont, Pierre. 
/// "An R-matrix package for coupled-channel problems in nuclear physics." 
/// Computer physics communications 200 (2016): 199-219. 
/// Calculates the on-shell R, S, T and K-matrices, as well as the scattering 
/// wavefunction in the channel.
class RMatrixSolverSingleChannel {
private:
  using Matrix = Eigen::Matrix<std::complex<double>,N,N>;
  Matrix Cinv;
  
  const Channel& channel;
  LagrangeLegendreBasis<N> basis;

public:
  /// @tparam callable (double) -> coomplex<double> evaluating a 
  /// radially weighted potential in the lth partial wave channel: 
  /// r * V(r,r') * r' , for a potential V(r,r') in MeV.
  /// If potential is local, V(r,r') should return 0 when r != r' 
  template<class Potential>
  RMatrixSolverSingleChannel(const Channel& channel, Potential potential)
    : Cinv(), channel(channel), basis(channel.radius)
    { 

      // set forward matrix C
      Matrix C;

      const auto mu = channel.reduced_mass;
      const auto l  = channel.l;
      
      // Eq. (6.10) in Descouvemont, 2016
      for (unsigned int n = 0; n < N; ++n) {
        for (unsigned int m = 0; m < N; ++m) {
          C(n,m) = basis.KE_Bloch_matrix_element(mu,l,n,m) 
                 + basis.non_local_matrix_element(potential,n,m) 
                 - channel.energy;
        }
      }
      
      // invert C to solve the system 
      Cinv = C.inverse();  
    }

  
  RMatrixSolverSingleChannel& operator=(const RMatrixSolverSingleChannel<N>&) {
    return *this;
  };
  RMatrixSolverSingleChannel& operator=(RMatrixSolverSingleChannel<N>&&) {
    return *this;
  };
  RMatrixSolverSingleChannel(const RMatrixSolverSingleChannel<N>&) = default;
  RMatrixSolverSingleChannel(RMatrixSolverSingleChannel<N>&&)      = default;

  struct Solution {
    std::complex<double> R, S, T, K;
    std::vector<std::complex<double>>  wvfxn;
  };

  /// @brief R-Matrix element for channel
  std::complex<double> rmatrix() const {
    using constants::hbar;
    const auto a = channel.radius;
    const auto mu = channel.reduced_mass;
    
    std::complex<double> R = 0;
    
    for (unsigned int n = 0; n < N; ++n) {
      for (unsigned int m = 0; m < N; ++m) {
        R += basis.f(n,a) * Cinv(n,m) * basis.f(m,a);
      }
    }
    return hbar*hbar / (2 * mu * a) * R;
  }
  
  /// @brief S-Matrix element for channel
  std::complex<double> smatrix(std::complex<double> R) const {
    const auto k = channel.k;
    const std::complex<double> in   = channel.asymptotic_wvfxn_in;
    const std::complex<double> inp  = channel.asymptotic_wvfxn_deriv_in;
    const std::complex<double> out  = channel.asymptotic_wvfxn_out;
    const std::complex<double> outp = channel.asymptotic_wvfxn_deriv_out;

    return (in - k * R * inp )/(out - k * R * outp);
  }
  
  /// @brief S-Matrix element for channel
  std::complex<double> smatrix() const {
    return smatrix(rmatrix());
  }
  
  /// @brief T-Matrix element for channel
  std::complex<double> tmatrix(std::complex<double> S) const {
    using constants::i;
    return i*S - i;
  }
  std::complex<double> tmatrix() const {
    return tmatrix(smatrix());
  }

  /// @brief K-matrix element of  channel
  /// defined as Caley transform of S-Matrix
  std::complex<double> kmatrix(std::complex<double> S) const {
    using constants::i;
    return i*(1.-S)/(1.+S);
  }
  std::complex<double> kmatrix() const {
    return kmatrix(smatrix());
  }

  /// @brief calculate the wavefunction within the scattering radius
  /// @param r radial grid, must be strictly increasing, and within [0,a]
  std::vector<std::complex<double>> wvfxn(const std::vector<double>& r ) const {
    using constants::hbar;
    const auto a = channel.radius;
    const auto mu = channel.reduced_mass;
    const auto hm =  hbar*hbar / (2 * mu * a);
    
    assert(r.front() >  0);
    assert(r.back() <= a);
    
    std::vector<std::complex<double>> w;

    for (unsigned int i = 0; i < r.size(); ++ i) {
      for (unsigned int n = 0; n < N; ++n) {
        for (unsigned int m = 0; m < N; ++m) {
          w[i] = hm * basis.f(n,a) * Cinv(n,m) * basis.f(m,r[i]);
        }
      }
    }
    return w;
  }

  Solution solve(const std::vector<double>& r) const {
    const auto R = rmatrix();
    const auto S = smatrix(R);
    return Solution{ R, S, tmatrix(S), kmatrix(S), wvfxn(r)  };
  }

};

}

#endif 
