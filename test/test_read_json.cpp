#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include "potential/params.hpp"

#include <fstream>


using namespace omplib;

TEST_CASE("Read KD json params from file") {
  
  SECTION("neutron") {
    // read default
    const auto path  = "KD_default.json";
    auto fstr  = std::ifstream(path);
    json pfile = json::parse(fstr);
    const auto kdn   = KD03Params<Proj::neutron>(pfile);

    // built in default
    const auto kdn_def = KD03Params<Proj::neutron>();

    REQUIRE( kdn.real_cent_r(66,156,189.23 ) == kdn_def.real_cent_r(66,156,189.23));
  }
  
  SECTION("neutron") {
    // read default
    const auto path  = "KD_default.json";
    auto fstr  = std::ifstream(path);
    json pfile = json::parse(fstr);
    const auto kdn   = KD03Params<Proj::proton>(pfile);

    // built in default
    const auto kdn_def = KD03Params<Proj::proton>();

    REQUIRE( kdn.real_cent_r(66,156,189.23 ) == kdn_def.real_cent_r(66,156,189.23));
  }
}

TEST_CASE("Read CH json params from file") {
  
  SECTION("neutron") {
    // read default
    const auto path  = "CH89_default.json";
    auto fstr  = std::ifstream(path);
    json pfile = json::parse(fstr);
    const auto chn   = CH89Params<Proj::neutron>(pfile);

    // built in default
    const auto chn_def = CH89Params<Proj::neutron>();

    REQUIRE( chn.real_cent_r(66,156,189.23 ) == chn_def.real_cent_r(66,156,189.23));
  }

  SECTION("proton") {
    // read default
    const auto path  = "CH89_default.json";
    auto fstr  = std::ifstream(path);
    json pfile = json::parse(fstr);
    const auto chn   = CH89Params<Proj::proton>(pfile);

    // built in default
    const auto chn_def = CH89Params<Proj::proton>();

    REQUIRE( chn.real_cent_r(66,156,189.23 ) == chn_def.real_cent_r(66,156,189.23));
  }
}

TEST_CASE("Read WLH json params from file") {
  
  SECTION("neutron") {
    // read default
    const auto path  = "WLH_mean.json";
    auto fstr  = std::ifstream(path);
    json pfile = json::parse(fstr);

    const auto wlh = WLH21Params<Proj::neutron>(pfile);

    const auto wlh_def = WLH21Params<Proj::neutron>(); 
    
    REQUIRE( wlh.real_cent_r(66,156,189.23 ) == wlh_def.real_cent_r(66,156,189.23));
  }
  
  SECTION("proton") {
    // read default
    const auto path  = "WLH_mean.json";
    auto fstr  = std::ifstream(path);
    json pfile = json::parse(fstr);

    const auto wlh = WLH21Params<Proj::proton>(pfile);

    const auto wlh_def = WLH21Params<Proj::proton>(); 
    
    REQUIRE( wlh.real_cent_r(66,156,189.23 ) == wlh_def.real_cent_r(66,156,189.23));
  }
}
