#include "potential/ch_params.hpp"

using namespace omplib;

template<>
CH89Params<Proj::neutron>::CH89Params()
  : OMParams<Proj::neutron>()
  ,  v_0(52.90)
  ,  v_e(-0.299)
  ,  v_asym(13.10)
  ,  r_0(-0.225)
  ,  r_A(1.250)
  ,  a0(0.690)

  ,  wv_0(7.8)
  ,  wve_0(35.0)
  ,  wv_ew(16.0)
  
  ,  ws_0(10.0)
  ,  ws_e0(36.0)
  ,  ws_ew(37.0)
  ,  ws_asym(18.0)
  ,  rw_0(-0.42)
  ,  rw_A(1.33)
  ,  aw(0.69)

  ,  vso_0(5.9)
  ,  rso_0(-1.2)
  ,  rso_A(1.34)
  ,  aso(0.63)
{};

CH89Params<Proj::proton>::CH89Params()
  : OMParams<Proj::proton>()
  , CH89Params<Proj::neutron>()
  , rc_0(0.12)
  , rc_A(1.24)
{};

CH89Params<Proj::proton>::CH89Params(json p)
  : OMParams<Proj::proton>()
  , CH89Params<Proj::neutron>(p)
  , rc_0(  p["CH89Coulomb_r_c_0"] )
  , rc_A(  p["CH89Coulomb_r_c"] )
{};


double omplib::CH89Params<Proj::proton>::Ec(int Z, int A, double erg) const {
  const double z = static_cast<double>(Z);
  const double a = static_cast<double>(A);

  return 6. * z * constants::e_sqr / (5 * real_coul_r(Z,A,erg));
}

double omplib::CH89Params<Proj::proton>::real_coul_r(
    int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return rc_0 + rc_A * pow(a, 1./3.);
}
