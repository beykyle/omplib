#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include "potential/params.hpp"

#include <fstream>


using namespace omplib;

TEST_CASE("Read json params from file") {
  
  // read default
  const auto path  = "KDGlobal.json";
  auto fstr  = std::ifstream(path);
  json pfile = json::parse(fstr);
  const auto kdn   = KoningDelaroche03<Proj::neutron>(pfile);

  // built in default
  const auto kdn_def = KoningDelaroche03<Proj::neutron>();

  REQUIRE( kdn.real_cent_r(66,156,189.23 ) == kdn_def.real_cent_r(66,156,189.23));
}
