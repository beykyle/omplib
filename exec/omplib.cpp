#include "potential/kd_params.hpp"
#include "potential/params.hpp"

#include <iostream>

constexpr auto n = omplib::Proj::neutron;
constexpr auto p = omplib::Proj::neutron;

int main(int argc, char** argv) {
  std::cout << "Running test problem!\n" << std::flush;
  
  auto kdn_def  = omplib::KoningDelaroche03<n>();
  auto kdn_uq   = omplib::KoningDelaroche03<n>::build_KDUQ();
  auto kdp_def  = omplib::KoningDelaroche03<p>();
  auto kdp_uq   = omplib::KoningDelaroche03<p>::build_KDUQ();

  static_assert(kdn_def.projectile == n);
  static_assert(kdn_uq.projectile  == n);
  
  static_assert(kdp_def.projectile == p);
  static_assert(kdp_uq.projectile  == p);


  return 0;
};
