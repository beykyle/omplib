#include "potential/ch_params.hpp"

using namespace omplib;
template<>
ChapelHill89<Proj::neutron>::ChapelHill89()
  : OMParams<Proj::neutron>()
{};

ChapelHill89<Proj::proton>::ChapelHill89()
  : OMParams<Proj::proton>()
  , ChapelHill89<Proj::neutron>()
{};

template<>
ChapelHill89<Proj::neutron>
ChapelHill89<Proj::neutron>::build_CHUQ(){};

ChapelHill89<Proj::proton>
ChapelHill89<Proj::proton>::build_CHUQ() {};

ChapelHill89<Proj::proton>::ChapelHill89(json p)
  : OMParams<Proj::proton>()
  , ChapelHill89<Proj::neutron>(p)
  , rc_0(  p["CH89Coulomb"]["r_c_0"] )
  , rc_A(  p["CH89Coulomb"]["r_c_A"] )
{};


double omplib::ChapelHill89<Proj::proton>::Ec(int Z, int A, double erg) const {
  const double z = static_cast<double>(Z);
  const double a = static_cast<double>(A);

  return 6. * z * e_sqr / (5 * real_coul_r(Z,A,erg));
}

double omplib::ChapelHill89<Proj::proton>::real_coul_V(
    int Z, int A, double erg) const {
  const double z = static_cast<double>(Z);
  const double a = static_cast<double>(A);
  return 6. * z * e_sqr / ( 5. * real_coul_r(z,a,erg) * pow(a,1./3.));
}

double omplib::ChapelHill89<Proj::proton>::real_coul_r(
    int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return rc_0 + rc_A * pow(a, -1./3.);
}
