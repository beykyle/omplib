#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "util/factorial.hpp"

using namespace omplib;

TEST_CASE("Factorial","[factorial]") {

  SECTION("Test factorial") {
    REQUIRE( factorial(0) == 1 );
    REQUIRE( factorial(1) == 1 );
    REQUIRE( factorial(2) == 2 );
    REQUIRE( factorial(5) == 120 );
    REQUIRE( factorial(7) == 5040 );
    REQUIRE( factorial(18) == 6.402373705728e15 );
  }

  SECTION("Test double factorial") {
    REQUIRE( doubleFactorial(0) == 1 );
    REQUIRE( doubleFactorial(1) == 1 );
    REQUIRE( doubleFactorial(2) == 2 );
    REQUIRE( doubleFactorial(5) == 15 );
    REQUIRE( doubleFactorial(7) == 105 );
    REQUIRE( doubleFactorial(18) == 185794560 );
  }
}
