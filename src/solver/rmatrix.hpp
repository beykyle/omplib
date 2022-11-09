#ifndef RMATRIX_HEADER
#define RMATRIX_HEADER

#include "util/constants.hpp"
#include "util/asymptotics.hpp"
#include "util/lagrange_legendre.hpp"
#include "util/types.hpp"

#include "potential/potential.hpp"

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
/// This type follows RAII; the system is initialized and inverted (solved) during 
/// construction. Once an instance has been constructed, the system has been solved,
/// and derived quantities (e.g. R/S/T/K-Matrices and the wavefunction in the basis)
/// can be calculated at little computational cost, using the calculate() function
class RMatrixSolverSingleChannel {
private:
  using Matrix = Eigen::Matrix<cmpl,N,N>;

  /// @brief Matrix on C^(NXN), define by H-E = (T + L + l*(l+1)/r^2 + V(r,r') - E)^-1 
  /// in the basis, where T is the derivative component of the kinetic energy, 
  /// and L is the Bloch operator
  Matrix Cinv;
  
  /// @brief projected basis in which (H-E) is inverted using Gauss-Legendre quadrature
  LagrangeLegendreBasis<N> basis;
  
  /// @brief channel information
  PData data;
  /// @brief wavefunctions exterior to the channel (asymptotic) and their derivatives
  const Channel::Asymptotics asym;
  

public:
  /// @tparam callable (real) -> cmpl evaluating a 
  /// radially weighted potential in the lth partial wave channel: 
  /// r * V(r,r') * r' , for a potential V(r,r') in MeV.
  /// If potential is local, r * V(r,r') * r'  = r V(r) * delta(r - r') 
  /// @brief Solves the scattering equation, with the provided potential in the 
  /// given channel.
  /// time complexity O(N^3)
  /// memory complexity O(N^2)
  template<class Potential>
  RMatrixSolverSingleChannel(
      const Channel channel, 
      const Channel::Energetics e, 
      const Channel::AngularMomentum am, 
      Potential potential)
    : Cinv()
    , basis( channel.radius )
    , data(  channel, e, am )
    , asym(  channel.set_angular_momentum(am, e.k) )
    { 
      using constants::hbar;
      const auto l    = am.l;
      const auto h2ma = e.h2ma;

      // pass channel info into the potential callable
      auto p = [&potential, this](real r, real rp) -> cmpl {
        return potential(r, rp, data);
      };
      
      // Eq. (6.10) in Baye, Daniel. 
      // "The Lagrange-mesh method." Physics reports 565 (2015): 1-107.
      // e.g. Cinv = (C-I*E)^-1
      for (unsigned int n = 0; n < N; ++n) {
        for (unsigned int m = 0; m < N; ++m) {
          
          Cinv(n,m) = h2ma * basis.KE_Bloch_matrix_element(l,n,m) 
                    + basis.non_local_matrix_element(p,n,m);
          
          if (n == m) Cinv(n,n) -= e.erg_cms;
        }
      }
      
      // invert C to solve the system 
      //TODO should be symmetric - probably more efficient inversion
      //TODO GPU/CUSolve?
      Cinv = Cinv.inverse();  
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
    /// \defgroup Matrices 
    /// matrices solving the scatter problem
    cmpl R, S, T, K;

    /// @brief coefficients of the reduced wavefunction in the basis
    std::array<cmpl,N> wvfxn;
  };

  /// @brief R-Matrix element for channel
  /// O(N^2)
  cmpl rmatrix() const {
    using constants::hbar;
    const auto a    = data.ch.radius;
    const auto h2ma = data.e.h2ma;
    
    cmpl R = 0;
    
    for (unsigned int n = 0; n < N; ++n) {
      for (unsigned int m = 0; m < N; ++m) {
        R += basis.f(n,a) * Cinv(n,m) * basis.f(m,a);
      }
    }
    return h2ma * R;
  }
  
  /// @brief S-Matrix element for channel
  cmpl smatrix(cmpl R) const {
    const auto k = data.e.k;
    const cmpl in   = asym.wvfxn_in;
    const cmpl inp  = asym.wvfxn_deriv_in;
    const cmpl out  = asym.wvfxn_out;
    const cmpl outp = asym.wvfxn_deriv_out;

    return (in - k * R * inp )/(out - k * R * outp);
  }
  
  /// @brief S-Matrix element for channel
  cmpl smatrix() const {
    return smatrix(rmatrix());
  }
  
  /// @brief T-Matrix element for channel
  cmpl tmatrix(cmpl S) const {
    using constants::i;
    return i*S - i;
  }
  cmpl tmatrix() const {
    return tmatrix(smatrix());
  }

  /// @brief K-matrix element of  channel
  /// defined as Caley transform of S-Matrix
  cmpl kmatrix(cmpl S) const {
    using constants::i;
    return i*(1.-S)/(1.+S);
  }
  cmpl kmatrix() const {
    return kmatrix(smatrix());
  }

  /// @brief calculate the reduced wavefunction u(r) = r*R(r), within the scattering radius
  /// @param r radial grid, must be strictly increasing, and within [0,a]
  /// @param S the S-Matrix
  std::vector<cmpl> wvfxn_rbasis(
      const std::vector<real>& r, cmpl S ) const {
    return wvfxn_rbasis(r, wvfxn(S));
  }
  
  /// @brief calculate the reduced wavefunction u(r) = r*R(r), within the scattering radius
  /// @param r radial grid, must be strictly increasing, and within [0,a]
  /// @param c coefficients of the reduced wavefunction in the basis
  std::vector<cmpl> wvfxn_rbasis(
      const std::vector<real>& r, std::array<cmpl,N> c ) const {

    assert(r.front() >= 0);
    assert(r.back()  <= data.ch.radius);
    
    std::vector<cmpl> w (r.size(), cmpl{0,0});

    for (unsigned int i = 0; i < r.size(); ++ i) {
      for (unsigned int m = 0; m < N; ++m) {
        w[i] += c[m] * basis.f(m,r);
      }
    }
    return w;
  }

  /// @brief Determine the N coefficients for the basis functions
  /// that uniquely determine the solution for the reduced wavefunction:
  /// u(r) = sum_m=0^N c_m basis.f(m,r)
  /// @param S the S-Matrix
  std::array<cmpl,N> wvfxn( cmpl S ) const {
    
    std::array<cmpl,N> c{N, cmpl{0,0}};
    
    using constants::hbar;
    using constants::i;
    const auto h2ma = data.e.h2ma;
    const auto a    = data.ch.radius;
    const auto k    = data.e.k;
    
    const auto wvext_deriv = 
      (asym.wvfxn_deriv_in + S * asym.wvfxn_deriv_out ) / (2 * k * i * a);
    
    for (unsigned int m = 0; m < N; ++m) {
      for (unsigned int n = 0; n < N; ++n) {
        c[m] += h2ma * wvext_deriv * basis.f(n,a) * Cinv(n,m);
      }
    }

    return c;
  }

  Solution calculate() const {
    const auto R = rmatrix();
    const auto S = smatrix(R);
    return Solution{ R, S, tmatrix(S), kmatrix(S), wvfxn(S) };
  }

};

}

#endif 
