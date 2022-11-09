#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "solver/channel.hpp"

using namespace omplib;
using Catch::Approx;

TEST_CASE("Build channel at varying energy") {

  using namespace omplib;
  using namespace omplib::constants;
  
  constexpr int A = 139;
  constexpr int Z = 54;

  // constants
  constexpr real projectile_mass = n_mass_amu; 
  constexpr real target_mass     = A; 
  
  constexpr real threshold = 0.; // MeV
  constexpr real ch_radius = 12; // fm
  
  /// s-wave
  constexpr int l   = 0;
  constexpr int S2  = 2; // S2 == (2*(1/2) +1) == dimension of spin 1/2 SU(2) representation
  // on 0+ ground state
  constexpr int J2  = 1; // J2 == (2*J +1) == dimension of SU(2) representation
  
  const auto chn  = Channel(threshold, ch_radius, 
                            projectile_mass, 0, target_mass, Z);
  
  const auto am   = Channel::AngularMomentum(S2, l, J2);

  SECTION("neutron 1MeV") {
    const auto e    = chn.set_erg_cms(1.0);
    const auto asym = chn.set_angular_momentum(am, e.k);
    
    REQUIRE(e.erg_lab           == Approx(1.00725658212791));
    REQUIRE(e.k                 == Approx(0.2189449404));
    REQUIRE(e.h2ma              == Approx(1.7374714662639255));
    REQUIRE(e.reduced_mass      == Approx(933.7788172169));
    REQUIRE(e.sommerfield_param == Approx(0.));
    
    REQUIRE(asym.wvfxn_in.real()         == Approx(2.24661470877110));
    REQUIRE(asym.wvfxn_in.imag()         == Approx(-3.97661722658389));
    REQUIRE(asym.wvfxn_deriv_in.real()   == Approx(-0.870660221853020));
    REQUIRE(asym.wvfxn_deriv_in.imag()   == Approx(-0.491884923618166));
    REQUIRE(asym.wvfxn_out.real()        == Approx(2.24661470877110));
    REQUIRE(asym.wvfxn_out.imag()        == Approx(3.97661722658389));
    REQUIRE(asym.wvfxn_deriv_out.real()  == Approx(-0.870660221853020));
    REQUIRE(asym.wvfxn_deriv_out.imag()  == Approx(0.491884923618166));
  }
  
  SECTION("neutron 100MeV") {
    const auto e    = chn.set_erg_cms(100.0);
    const auto asym = chn.set_angular_momentum(am, e.k);
    
    REQUIRE(e.erg_lab           == Approx(100.725658212791));
    REQUIRE(e.k                 == Approx(2.24505752889899));
    REQUIRE(e.h2ma              == Approx(1.57382520561950));
    REQUIRE(e.reduced_mass      == Approx(1030.87383743656));
    REQUIRE(e.sommerfield_param == Approx(0.));
    
    REQUIRE(asym.wvfxn_in.real()         == Approx(0.432955865628073));
    REQUIRE(asym.wvfxn_in.imag()         == Approx(-0.104645899583040));
    
    REQUIRE(asym.wvfxn_deriv_in.real()   == Approx(-0.234936064727312));
    REQUIRE(asym.wvfxn_deriv_in.imag()   == Approx(-0.972010825809283));
    
    REQUIRE(asym.wvfxn_out.real()        == Approx(0.432955865628073));
    REQUIRE(asym.wvfxn_out.imag()        == Approx(0.104645899583040));
    
    REQUIRE(asym.wvfxn_deriv_out.real()  == Approx(-0.234936064727312));
    REQUIRE(asym.wvfxn_deriv_out.imag()  == Approx(0.972010825809283));
  }
  
  /*
  SECTION("proton") {
    const auto chp   = Channel(threshold, erg_cms, ch_radius, 
                              constants::p_mass_amu, target_mass, Z, 1., 0, J2, pi);
    
    REQUIRE(chp.sommerfield_param == Approx(15.1412752890966));
  }
  */
}
