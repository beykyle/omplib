#ifndef LAGRANGE_LEGENDRE_HEADER
#define LAGRANGE_LEGENDRE_HEADER

#include <cmath>
#include <cassert>

#include "util/constants.hpp"
#include "util/gauss_legendre_weights.hpp"
#include "util/types.hpp"

namespace omplib {

///@tparam N is the number of basis functions
template<unsigned int N>
/// @brief functionality for the R-Matrix method on a Lagrange-Legendre basis, following:
/// Descouvemont, Pierre. 
/// "An R-matrix package for coupled-channel problems in nuclear physics." 
/// Computer physics communications 200 (2016): 199-219.
class LagrangeLegendreBasis {
public:
  using GL = GaussLegendre<N,0,1>;
private:
  ///@brief channel radius
  real a;

  static constexpr GL gl = GL{};

public:
  LagrangeLegendreBasis(int a): a(a) {};

  /// @brief computes the Lagrange-Legendre function of order n
  /// @param n referrs to the nth Lagrange Legendre function
  /// @param r is a radial distance in fm on [0,a]
  real f(unsigned int n, real r) const {
    
    assert(n <  N);
    assert(r <= a);
    
    const real f = (N+n)%2==0 ? 1. : -1.;
    const real xn = gl.abscissa[n];
    return f * r / ( a * xn ) * sqrt( a * xn * ( 1 - xn ) ) 
      * std::legendre(N, 2 * r / a - 1) / (r - a*xn); 
  }

  /// @tparam callable (real) -> complex<real> evaluating a radially weighted potential in 
  /// the lth partial wave channel: r * V(r), for a potential V(r) in MeV
  template<class Potential>
  /// @brief matrix element of an arbitrary local potential in the Lagrange Legendre basis, 
  /// approximated using Gauss-Legendre quadrature
  /// @param n is the order of the Lagrange-Legendre bra
  /// @param m is the order of the Lagrange-Legendre ket
  cmpl local_matrix_element(Potential p, unsigned int n, unsigned int m ) const {
    assert(n < N);
    assert(m < N);
    
    if (n != m) return 0;
    
    const real xn = gl.abscissa[n];
    return p(a*xn);
  }

  /// @tparam callable (real) -> complex<real> evaluating a radially weighted potential in 
  /// the lth partial wave channel: r * V(r,r') * r' , for a potential V(r,r') in MeV
  template<class Potential>
  /// @brief matrix element of an arbitrary non-local potential in the Lagrange Legendre basis, 
  /// approximated using Gauss-Legendre quadrature
  /// @param n is the order of the Lagrange-Legendre bra
  /// @param m is the order of the Lagrange-Legendre ket
  cmpl non_local_matrix_element(Potential p, unsigned int n, unsigned int m) const {
    assert(n < N);
    assert(m < N);

    const real xn       = gl.abscissa[n];
    const real xm       = gl.abscissa[m];
    const real lambda_n = gl.weights[n];
    const real lambda_m = gl.weights[m];

    return a * sqrt(lambda_n * lambda_m) * p(a*xn, a*xm);
  }

  /// @brief matrix elements of the radial kinetic energy for partial wave l
  /// plus the Bloch operator. These are symnmetric in n and m.
  /// @param l partial wave channel
  /// @param n is the order of the Lagrange-Legendre bra
  /// @param m is the order of the Lagrange-Legendre ket
  real KE_Bloch_matrix_element(unsigned int l, unsigned int n, unsigned int m) const {
    assert(n < N);
    assert(m < N);
    
    using constants::hbar;
    
    const real xn = gl.abscissa[n];

    if (n == m) {
      return l*(l+1)/(a*a*xn*xn) // centrifugal
        * (4. * N * N + 4. * N + 3.) * (xn * (1. - xn) - 6. * xn + 1.)  // kinetic
        / (3. * a * a * xn * xn * (1. - xn) * (1. - xn) );
    }
    
    const real xm = gl.abscissa[m];
    const real f  = (n+m)%2==0 ? 1. : -1.;
    
    // centrifugal only on diag
    return f / (xn * xm  * (1. - xn) *(1. - xm)) 
      * (
          N*N + N + 1. + (xn + xm - 2 * xn * xm)/( ( xn - xm ) * ( xn - xm ) ) 
          - 1. / ( 1 - xn ) - 1 / ( 1 - xm ) 
        );
  }

};

}

#endif 
