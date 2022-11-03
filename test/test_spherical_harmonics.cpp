#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "util/spherical_harmonics.hpp"
#include "util/constants.hpp"

#include <iostream>
#include <cmath>

using namespace omplib;
using Catch::Approx;

TEST_CASE("SphericalHarmonics","[sphericalHarmonics]") {
  SECTION("Test associated legendre polynomials") {

    // P^0_0(x) = 1.0
    REQUIRE( evalLegendrePolynomial(0,0,1.0) == 1.0 );
    REQUIRE( evalLegendrePolynomial(0,0,0.2) == 1.0 );
    REQUIRE( evalLegendrePolynomial(0,0,-0.7) == 1.0 );

    // P^1_0(x) = x
    REQUIRE( evalLegendrePolynomial(1,0,1.0) == 1.0 );
    REQUIRE( evalLegendrePolynomial(1,0,0.8) == 0.8 );
    REQUIRE( evalLegendrePolynomial(1,0,-0.6) == -0.6 );

    // P^1_1(x) = -(1-x^2)^(1/2)
    REQUIRE( evalLegendrePolynomial(1,1,1.0) == 0.0 );
    REQUIRE( evalLegendrePolynomial(1,1,0.2) == -sqrt(1. - 0.2*0.2) );
    REQUIRE( evalLegendrePolynomial(1,1,-0.4) == -sqrt(1. - 0.4*0.4) );

    // P^(-1)_1(x) = 0.5*(1-x^2)^(1/2)
    REQUIRE( evalLegendrePolynomial(1,-1,1.0) == 0.0 );
    REQUIRE( evalLegendrePolynomial(1,-1,0.2) == 0.5*sqrt(1. - 0.2*0.2) );
    REQUIRE( evalLegendrePolynomial(1,-1,-0.4) == 0.5*sqrt(1. - 0.4*0.4) );

    // P^(2)_2(x) = 3*(1-x^2)
    REQUIRE( evalLegendrePolynomial(2,2,1.0) == 0.0 );
    REQUIRE( evalLegendrePolynomial(2,2,0.2) == 3.*(1. - 0.2*0.2) );
    REQUIRE( evalLegendrePolynomial(2,2,-0.4) == 3.*(1. - 0.4*0.4) );

    // P^(-2)_2(x) = 1/8*(1-x^2)
    REQUIRE( evalLegendrePolynomial(2,-2,1.0) == 0.0 );
    REQUIRE( evalLegendrePolynomial(2,-2,0.2) == 1./8.*(1. - 0.2*0.2) );
    REQUIRE( evalLegendrePolynomial(2,-2,-0.4) == 1./8.*(1. - 0.4*0.4) );
  }

  SECTION("Test spherical harmonic function") {
    // P^0_0(x) = 1.0, let omega = pi/4 -> cos(m*omega) = 1
    REQUIRE( realSphericalHarmonic(0,0,1.0,0.) == 1.0 );
    REQUIRE( realSphericalHarmonic(0,0,0.2,0.) == 1.0 );
    REQUIRE( realSphericalHarmonic(0,0,-0.7,0.) == 1.0 );

    // P^1_0(x) = x, C_{l,m} = 3, let omega = pi/4 -> cos(m*omega) = 1
    REQUIRE( realSphericalHarmonic(1,0,1.0,0.) == sqrt(3) );
    REQUIRE( realSphericalHarmonic(1,0,0.8,0.) == sqrt(3)*0.8 );
    REQUIRE( realSphericalHarmonic(1,0,-0.6,0.) == sqrt(3)*-0.6 );

    // P^1_1(x) = -(1-x^2)^(1/2), C_{l,m} = 3/2, let omega = pi/4 -> cos(m*omega) = 1/sqrt(2)
    REQUIRE( realSphericalHarmonic(1,1,1.0,constants::pi/4) == 0.0 );
    REQUIRE( realSphericalHarmonic(1,1,0.2,constants::pi/4) == Approx( sqrt(1.5) * -sqrt(1. - 0.2*0.2) * sqrt(2)*0.5 ));
    REQUIRE( realSphericalHarmonic(1,1,-0.4,constants::pi/4) == Approx( sqrt(1.5) * -sqrt(1. - 0.4*0.4) * sqrt(2)*0.5 ));

    // P^(-1)_1(x) = 0.5*(1-x^2)^(1/2), C_{l,m} = 6, let omega = pi/4 -> cos(m*omega) = 1/sqrt(2)
    REQUIRE( realSphericalHarmonic(1,-1,1.0,constants::pi/4) == 0.0 );
    REQUIRE( realSphericalHarmonic(1,-1,0.2,constants::pi/4) == Approx( sqrt(6) * 0.5*sqrt(1. - 0.2*0.2) * sqrt(2)*0.5 ));
    REQUIRE( realSphericalHarmonic(1,-1,-0.4,constants::pi/4) == Approx( sqrt(6) * 0.5*sqrt(1. - 0.4*0.4) * sqrt(2)*0.5));

    // P^(2)_2(x) = 3*(1-x^2), C_{l,m} = 5/24, let omega = pi/6 -> cos(m*omega) = 1/2
    REQUIRE( realSphericalHarmonic(2,2,1.0,constants::pi/6) == 0.0 );
    REQUIRE( realSphericalHarmonic(2,2,0.2,constants::pi/6) == Approx( sqrt(5./24.) * 3.*(1. - 0.2*0.2) * 0.5 ));
    REQUIRE( realSphericalHarmonic(2,2,-0.4,constants::pi/6) == Approx( sqrt(5./24.) * 3.*(1. - 0.4*0.4) * 0.5 ));

    // P^(-2)_2(x) = 1/8*(1-x^2), C_{l,m} = 120, let omega = pi/6 -> cos(m*omega) = 1/2
    REQUIRE( realSphericalHarmonic(2,-2,1.0,constants::pi/6) == 0.0 );
    REQUIRE( realSphericalHarmonic(2,-2,0.2,constants::pi/6) == Approx( sqrt(120.) * 1./8.*(1. - 0.2*0.2) * 0.5 ));
    REQUIRE( realSphericalHarmonic(2,-2,-0.4,constants::pi/6) == Approx( sqrt(120.) * 1./8.*(1. - 0.4*0.4) * 0.5 ));
  }
}
