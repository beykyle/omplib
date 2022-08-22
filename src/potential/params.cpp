#include "potential/params.hpp"
#include "util/constants.hpp"

// KD03
using namespace omplib;
using omplib::Projectile;

double KoningDelaroche03<Projectile::proton>::real_coul_V(
    int Z, int A, double erg) const {
  const double z = static_cast<double>(Z);
  const double a = static_cast<double>(A);
  return 6. * z * e * e / ( 5. * real_coul_r(z,a,erg) * pow(a,1./3.));
}

double KoningDelaroche03<Projectile::proton>::real_coul_r(
    int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return rc_0 + rc_A * pow(a, -1./3.) + rc_A2 * pow(a, -5./3.);
}
