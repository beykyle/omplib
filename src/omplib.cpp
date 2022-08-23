#include "potential/params.hpp"

#include <iostream>

using omplib::Proj;

int main(int argc, char** argv) {
  std::cout << "Running test problem!\n" << std::flush;
  
  auto kdn_def = omplib::KoningDelaroche03<Proj::neutron>();
  auto kdn_uq = omplib::KDUQ<Proj::neutron>();
  auto kdp_def = omplib::KoningDelaroche03<Proj::proton>();
  auto kdp_uq = omplib::KDUQ<Proj::proton>();

  static_assert(kdn_def.projectile == Proj::neutron);
  static_assert(kdp_def.projectile == Proj::proton);
  static_assert(kdn_uq.projectile == Proj::neutron);
  static_assert(kdp_uq.projectile == Proj::proton);


  return 0;
};
