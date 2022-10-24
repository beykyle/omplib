#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include "potential/params.hpp"

#include <fstream>


using namespace omplib;

TEST_CASE("Read KD json params from file") {
  
  SECTION("neutron") {
    // read default
    const auto path  = "KDGlobal.json";
    auto fstr  = std::ifstream(path);
    json pfile = json::parse(fstr);
    const auto kdn   = KoningDelaroche03<Proj::neutron>(pfile);

    // built in default
    const auto kdn_def = KoningDelaroche03<Proj::neutron>();

    REQUIRE( kdn.real_cent_r(66,156,189.23 ) == kdn_def.real_cent_r(66,156,189.23));
  }
  
  SECTION("neutron") {
    // read default
    const auto path  = "KDGlobal.json";
    auto fstr  = std::ifstream(path);
    json pfile = json::parse(fstr);
    const auto kdn   = KoningDelaroche03<Proj::proton>(pfile);

    // built in default
    const auto kdn_def = KoningDelaroche03<Proj::proton>();

    REQUIRE( kdn.real_cent_r(66,156,189.23 ) == kdn_def.real_cent_r(66,156,189.23));
  }
}

TEST_CASE("Read CH json params from file") {
  
  SECTION("neutron") {
    // read default
    const auto path  = "CH89_default.json";
    auto fstr  = std::ifstream(path);
    json pfile = json::parse(fstr);
    const auto chn   = ChapelHill89<Proj::neutron>(pfile);

    // built in default
    const auto chn_def = ChapelHill89<Proj::neutron>();

    REQUIRE( chn.real_cent_r(66,156,189.23 ) == chn_def.real_cent_r(66,156,189.23));
  }

  SECTION("proton") {
    // read default
    const auto path  = "CH89_default.json";
    auto fstr  = std::ifstream(path);
    json pfile = json::parse(fstr);
    const auto chn   = ChapelHill89<Proj::proton>(pfile);

    // built in default
    const auto chn_def = ChapelHill89<Proj::proton>();

    REQUIRE( chn.real_cent_r(66,156,189.23 ) == chn_def.real_cent_r(66,156,189.23));
  }
}
