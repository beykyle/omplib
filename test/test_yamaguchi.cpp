#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "solver/rmatrix.hpp"
#include "potential/potential.hpp"

using namespace omplib;
using Catch::Approx;
using Solver = RMatrixSolverSingleChannel<20>;


TEST_CASE("Yamaguchi analytic s-wave phase shift") {
 // test of R-Matrix solver against analytic potential
 // see Table 16 of 
 // Baye, Daniel. "The Lagrange-mesh method." Physics reports 565 (2015): 1-107.

  // potential
  const auto y = Yamaguchi();
  auto p = [&y](auto r, auto rp, auto d) {
    return y.eval_reduced(r,rp,d);
  };

  const auto ch  = Channel(0., 15, constants::n_mass_amu, 0, frac(1,2), 
                           constants::p_mass_amu, 1);

  // S-Wave, 0+
  const auto am   = Channel::AngularMomentum(2);
  SECTION("0.1MeV") {
    const auto e  = ch.set_erg_cms(0.1);
    const auto k  = e.k;

    REQUIRE( Solver(ch, e, am, p).kmatrix() == y.analytic_swave_kmatrix(k) );
  }
  
  SECTION("10MeV") {
    const auto e  = ch.set_erg_cms(10.0);
    const auto k  = e.k;

    REQUIRE( Solver(ch, e, am, p).kmatrix() == y.analytic_swave_kmatrix(k) );
  }
}
