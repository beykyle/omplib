#include "potential/ch_params.hpp"

using namespace omplib;

template<>
CH89Params<Proj::neutron>::CH89Params()
  : OMParams<Proj::neutron>()
  ,  r_0(-0.225)
  ,  r_A(1.250)
  ,  a0(0.690)
  ,  rw_0(-0.42)
  ,  rw_A(1.33)
  ,  aw(0.69)
  ,  rso_0(-1.2)
  ,  rso_A(1.34)
  ,  aso(0.63)
  
  ,  v_0(52.90)
  ,  v_e(-0.299)
  ,  v_asym(13.10)
  ,  wv_0(7.8)
  ,  wve_0(35.0)
  ,  wv_ew(16.0)
  ,  ws_0(10.0)
  ,  ws_asym(18.0)
  ,  ws_e0(36.0)
  ,  ws_ew(37.0)

  ,  vso_0(5.9)
{};

CH89Params<Proj::proton>::CH89Params()
  : CH89Params<Proj::neutron>()
  , OMParams<Proj::proton>()
  , rc_0(0.12)
  , rc_A(1.24)
{};

CH89Params<Proj::proton>::CH89Params(json p)
  : CH89Params<Proj::neutron>(p)
  , OMParams<Proj::proton>()
  , rc_0(  p["CH89Coulomb_r_c_0"] )
  , rc_A(  p["CH89Coulomb_r_c"] )
{};


real omplib::CH89Params<Proj::proton>::Ec(int Z, int A, real erg) const {
  const real z = static_cast<real>(Z);

  return 6. * z * constants::e_sqr / (5 * real_coul_r(Z,A,erg));
}

real omplib::CH89Params<Proj::proton>::real_coul_r(
    int, int A, real) const {
  const real a = static_cast<real>(A);
  return rc_0 + rc_A * pow(a, 1./3.);
}
