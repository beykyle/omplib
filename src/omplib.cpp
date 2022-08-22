#include "potential/params.hpp"

#include <iostream>

using omplib::Projectile;

int main(int argc, char** argv) {
  std::cout << "Running test problem!\n" << std::flush;
  
  auto kdn = omplib::KoningDelaroche03<Projectile::neutron>{};
  auto kdp = omplib::KoningDelaroche03<Projectile::proton>{};

  static_assert(kdn.projectile == Projectile::neutron);
  static_assert(kdp.projectile == Projectile::proton);


    

  return 0;
};
