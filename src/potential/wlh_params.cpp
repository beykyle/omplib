#include "potential/wlh_params.hpp"

using namespace omplib;

WLH21<Proj::proton>::WLH21(json p) {
  v0 = p["WLHReal_V0_p"];
  v1 = p["WLHReal_V1_p"];
  v2 = p["WLHReal_V2_p"];
  v3 = p["WLHReal_V3_p"];
  v4 = p["WLHReal_V4_p"];
  v5 = p["WLHReal_V5_p"];
  v6 = p["WLHReal_V6_p"];
  r0 = p["WLHReal_r0_p"];
  r1 = p["WLHReal_r1_p"];
  r2 = p["WLHReal_r2_p"];
  r3 = p["WLHReal_r3_p"];
  a0 = p["WLHReal_a0_p"];
  a1 = p["WLHReal_a1_p"];
  a2 = p["WLHReal_a2_p"];
  a3 = p["WLHReal_a3_p"];
  a4 = p["WLHReal_a4_p"];
  w0 = p["WLHImagVolume_W0_p"];
  w1 = p["WLHImagVolume_W1_p"];
  w2 = p["WLHImagVolume_W2_p"];
  w3 = p["WLHImagVolume_W3_p"];
  w4 = p["WLHImagVolume_W4_p"];
  rw0 = p["WLHImagVolume_r0_p"];
  rw1 = p["WLHImagVolume_r1_p"];
  rw2 = p["WLHImagVolume_r2_p"];
  rw3 = p["WLHImagVolume_r3_p"];
  rw4 = p["WLHImagVolume_r4_p"];
  rw5 = p["WLHImagVolume_r5_p"];
  aw0 = p["WLHImagVolume_a0_p"];
  aw1 = p["WLHImagVolume_a1_p"];
  aw2 = p["WLHImagVolume_a2_p"];
  aw3 = p["WLHImagVolume_a3_p"];
  aw4 = p["WLHImagVolume_a4_p"];
  d0 = p["WLHImagSurface_W0_p"];
  d1 = p["WLHImagSurface_W1_p"];
  d2 = p["WLHImagSurface_W2_p"];
  d3 = p["WLHImagSurface_W3_p"];
  rs0 = p["WLHImagSurface_r0_p"];
  rs1 = p["WLHImagSurface_r1_p"];
  rs2 = p["WLHImagSurface_r2_p"];
  as0 = p["WLHImagSurface_a0_p"];
  vso_0 = p["WLHRealSpinOrbit_V0_p"];
  vso_1 = p["WLHRealSpinOrbit_V1_p"];
  rso_0 = p["WLHRealSpinOrbit_r0_p"];
  rso_1 = p["WLHRealSpinOrbit_r1_p"];
  aso_0 = p["WLHRealSpinOrbit_a0_p"];
  aso_1 = p["WLHRealSpinOrbit_a1_p"];
}

template<>
WLH21<Proj::neutron>::WLH21(json p):
  v0(p["WLHReal_V0_n"]),
  v1(p["WLHReal_V1_n"]),
  v2(p["WLHReal_V2_n"]),
  v3(p["WLHReal_V3_n"]),
  v4(p["WLHReal_V4_n"]),
  v5(p["WLHReal_V5_n"]),
  v6(p["WLHReal_V6_n"]),
  r0(p["WLHReal_r0_n"]),
  r1(p["WLHReal_r1_n"]),
  r2(p["WLHReal_r2_n"]),
  r3(p["WLHReal_r3_n"]),
  a0(p["WLHReal_a0_n"]),
  a1(p["WLHReal_a1_n"]),
  a2(p["WLHReal_a2_n"]),
  a3(p["WLHReal_a3_n"]),
  a4(p["WLHReal_a4_n"]),
  w0(p["WLHImagVolume_W0_n"]),
  w1(p["WLHImagVolume_W1_n"]),
  w2(p["WLHImagVolume_W2_n"]),
  w3(p["WLHImagVolume_W3_n"]),
  w4(p["WLHImagVolume_W4_n"]),
  rw0(p["WLHImagVolume_r0_n"]),
  rw1(p["WLHImagVolume_r1_n"]),
  rw2(p["WLHImagVolume_r2_n"]),
  rw3(p["WLHImagVolume_r3_n"]),
  rw4(p["WLHImagVolume_r4_n"]),
  rw5(p["WLHImagVolume_r5_n"]),
  aw0(p["WLHImagVolume_a0_n"]),
  aw1(p["WLHImagVolume_a1_n"]),
  aw2(p["WLHImagVolume_a2_n"]),
  aw3(p["WLHImagVolume_a3_n"]),
  aw4(p["WLHImagVolume_a4_n"]),
  d0(p["WLHImagSurface_W0_n"]),
  d1(p["WLHImagSurface_W1_n"]),
  d2(p["WLHImagSurface_W2_n"]),
  d3(p["WLHImagSurface_W3_n"]),
  rs0(p["WLHImagSurface_r0_n"]),
  rs1(p["WLHImagSurface_r1_n"]),
  rs2(p["WLHImagSurface_r2_n"]),
  as0(p["WLHImagSurface_a0_n"]),
  vso_0(p["WLHRealSpinOrbit_V0_n"]),
  vso_1(p["WLHRealSpinOrbit_V1_n"]),
  rso_0(p["WLHRealSpinOrbit_r0_n"]),
  rso_1(p["WLHRealSpinOrbit_r1_n"]),
  aso_0(p["WLHRealSpinOrbit_a0_n"]),
  aso_1(p["WLHRealSpinOrbit_a1_n"])
{}
