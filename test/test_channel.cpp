#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "solver/channel.hpp"

using namespace omplib;
using Catch::Approx;

TEST_CASE("Build simple channel") {

  using namespace omplib;
  using namespace omplib::constants;
  
  constexpr int A = 139;
  constexpr int Z = 54;

  // constants
  constexpr real projectile_mass = n_mass_amu; 
  constexpr real target_mass     = A; 
  
  constexpr real threshold = 0.; // MeV
  constexpr real erg_cms   = 1.; // MeV
  constexpr real ch_radius = 12; // fm
  
  // on 0+ ground state
  constexpr int J2  = 1; // J2 == (2*J +1) == dimension of SU(2) representation
  constexpr auto pi = Parity::even;

  SECTION("neutron 1MeV") {
    const auto chn   = Channel(threshold, erg_cms, ch_radius, 
                              projectile_mass, 0, target_mass, Z, 0, 0, 2);
    
    REQUIRE(chn.Tlab         == Approx(1.00725658212791));
    REQUIRE(chn.k            == Approx(0.2189449404));
    REQUIRE(chn.h2ma         == Approx(1.7374714662639255));
    REQUIRE(chn.reduced_mass      == Approx(933.7788172169));
    REQUIRE(chn.sommerfield_param == Approx(0.));
    
    REQUIRE(chn.asymptotic_wvfxn_in.real()         == Approx(2.24661470877110));
    REQUIRE(chn.asymptotic_wvfxn_in.imag()         == Approx(-3.97661722658389));
    
    REQUIRE(chn.asymptotic_wvfxn_deriv_in.real()   == Approx(-0.870660221853020));
    REQUIRE(chn.asymptotic_wvfxn_deriv_in.imag()   == Approx(-0.491884923618166));
    
    REQUIRE(chn.asymptotic_wvfxn_out.real()        == Approx(2.24661470877110));
    REQUIRE(chn.asymptotic_wvfxn_out.imag()        == Approx(3.97661722658389));
    
    REQUIRE(chn.asymptotic_wvfxn_deriv_out.real()  == Approx(-0.870660221853020));
    REQUIRE(chn.asymptotic_wvfxn_deriv_out.imag()  == Approx(0.491884923618166));
  }
  
  SECTION("neutron 100MeV") {
    const auto chn   = Channel(threshold, 100, ch_radius, 
                              projectile_mass, 0, target_mass, Z, 0, 0, 2);
    
    REQUIRE(chn.Tlab         == Approx(100.725658212791));
    REQUIRE(chn.k            == Approx(2.24505752889899));
    REQUIRE(chn.h2ma         == Approx(1.57382520561950));
    REQUIRE(chn.reduced_mass      == Approx(1030.87383743656));
    REQUIRE(chn.sommerfield_param == Approx(0.));
    
    REQUIRE(chn.asymptotic_wvfxn_in.real()         == Approx(0.432955865628073));
    REQUIRE(chn.asymptotic_wvfxn_in.imag()         == Approx(-0.104645899583040));
    
    REQUIRE(chn.asymptotic_wvfxn_deriv_in.real()   == Approx(-0.234936064727312));
    REQUIRE(chn.asymptotic_wvfxn_deriv_in.imag()   == Approx(-0.972010825809283));
    
    REQUIRE(chn.asymptotic_wvfxn_out.real()        == Approx(0.432955865628073));
    REQUIRE(chn.asymptotic_wvfxn_out.imag()        == Approx(0.104645899583040));
    
    REQUIRE(chn.asymptotic_wvfxn_deriv_out.real()  == Approx(-0.234936064727312));
    REQUIRE(chn.asymptotic_wvfxn_deriv_out.imag()  == Approx(0.972010825809283));
  }
  
  /*
  SECTION("proton") {
    const auto chp   = Channel(threshold, erg_cms, ch_radius, 
                              constants::p_mass_amu, target_mass, Z, 1., 0, J2, pi);
    
    REQUIRE(chp.sommerfield_param == Approx(15.1412752890966));
  }
  */
}
