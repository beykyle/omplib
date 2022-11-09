#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "util/gauss_legendre_weights.hpp"

using namespace omplib;
using Catch::Approx;

TEST_CASE("regular GL weights") {
  
  SECTION("N=10") {
    constexpr auto gl = GaussLegendre<10>{};
    
    REQUIRE(gl.weights[0] == Approx(0.06667134));
    REQUIRE(gl.weights[3] == Approx(0.26926672));
    REQUIRE(gl.weights[8] == Approx(0.14945135));
    
    REQUIRE(gl.abscissa[0] == Approx(-0.97390653));
    REQUIRE(gl.abscissa[3] == Approx(-0.43339539));
    REQUIRE(gl.abscissa[8] == Approx(0.86506337));
  }
}

TEST_CASE("GL weights shift reduce") {
  
  SECTION("N=10") {
    constexpr auto gl = GaussLegendre<10,0,1>{};
    
    REQUIRE(gl.weights[0] == Approx(0.03333567));
    REQUIRE(gl.weights[3] == Approx(0.13463336));
    REQUIRE(gl.weights[8] == Approx(0.07472567));
    
    REQUIRE(gl.abscissa[0] == Approx(0.01304674));
    REQUIRE(gl.abscissa[3] == Approx(0.2833023));
    pEQUIRE(gl.abscissa[8] == Approx(0.93253168));
  }
}
