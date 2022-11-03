#include "potential/wlh_params.hpp"

using namespace omplib;

WLH21Params<Proj::proton>::WLH21Params(json p) {
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
WLH21Params<Proj::neutron>::WLH21Params(json p):
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

WLH21Params<Proj::proton>::WLH21Params()
  : WLH21Params<Proj::neutron>(), OMParams<Proj::proton>()
{
  v0     = 53.6540319643;
  v1     = 0.300042293;
  v2     = -0.0001192337;
  v3     = 2.1471e-06;
  v4     = 12.918410608;
  v5     = 0.1330858829;
  v6     = 0.0003328253;
  r0     = 1.3041479419;
  r1     = 0.2938537968;
  r2     = 0.0008039876;
  r3     = 6.8659e-06;
  a0     = 0.779602661;
  a1     = 0.0001182833;
  a2     = 5.4588e-06;
  a3     = 0.0445486378;
  a4     = 0.5571679352;
  w0     = 4.3701953072;
  w1     = 0.2267714208;
  w2     = 0.000338132;
  w3     = 11.4340928806;
  w4     = 0.077068475;
  rw0    = 0.6541271675;
  rw1    = 46.0846969827;
  rw2    = 0.6099691135;
  rw3    = 58.2530624354;
  rw4    = 0.6418349689;
  rw5    = 1.5491e-06;
  aw0    = 0.5616542778;
  aw1    = 0.2708915659;
  aw2    = 14.3691362537;
  aw3    = 0.4184493049;
  aw4    = 0.0024314838;
  d0     = 0.8982409214;   
  d1     = 0.0351390393;   
  d2     = 0.8135662888;   
  d3     = 0.1214099645;   
  rs0    = 0.9505407272;
  rs1    = 0.4391797103;
  rs2    = 0.003947669;
  as0    = 0.7649106373;
  vso_0  = 9.6195980947;
  vso_1  = 0.0085112416;
  rso_0  = 1.276250854;
  rso_1  = 0.8685769478;
  aso_0  = 0.8032951951;
  aso_1  = 0.0003508356;
}
