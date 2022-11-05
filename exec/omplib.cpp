#include "potential/ch_params.hpp"
#include "potential/kd_params.hpp"
#include "potential/params.hpp"
#include "potential/wlh_params.hpp"

#include "potential/potential.hpp"

#include "solver/rmatrix.hpp"

#include <iomanip>
#include <ios>
#include <iostream>
#include <fstream>

constexpr auto n = omplib::Proj::neutron;
constexpr auto p = omplib::Proj::neutron;

int main(int argc, char** argv) {
  
  auto kdn_uq   = omplib::KD03Params<n>::build_KDUQ();
  auto wlh_mean = omplib::WLH21Params<n>();

  // R-Matrix
  const auto p = OMP(96,142,wlh_mean);

  const auto ch     = Channel(0., 1., 12, 10., 96, 0, 0, 0, Parity::even);
  const auto solver = RmatrixSolverSingleChannel<10>(ch, p);
  const auto soln   = solver.solve();
  
  // print pot stuff on energy grid
  constexpr auto erg_min = 0.01;
  constexpr auto erg_max = 10.;
  constexpr auto range = erg_max - erg_min;
  constexpr auto e_grid_sz = 500;
  
  auto e_grid = std::array<double,e_grid_sz> {};
  for (int i = 0; i < e_grid_sz; ++i) {
    e_grid[i] = erg_min 
      + range * static_cast<double>(i) / static_cast<double>(e_grid_sz);
  }

  // mass 144 isotopes
  using Isotope = std::pair<int,int>;
  constexpr size_t niso = 6;
  
  constexpr auto isotopes =  std::array<Isotope,niso>{
    Isotope{58,86},Isotope{57,87},Isotope{56,88},
    Isotope{55,89},Isotope{54,90},Isotope{53,91}
  };

  std::cout << "Calculating surface potential depths for " << niso << " isotopes.\n" << std::flush;

  auto out_wlh = std::ofstream("./cmpl_vs_wlh.csv");
  auto out_kd = std::ofstream("./cmpl_vs_kd.csv");

  out_wlh << "E" << "\t";
  out_kd  << "E" << "\t";
  for (const auto& [Z,A] : isotopes) {
    out_wlh << Z << "_" << A << "\t";
    out_kd  << Z << "_" << A << "\t";
  }
  out_wlh << "\n";
  out_kd  << "\n";

  for (int j = 0; j < e_grid_sz; ++j) {
    
    const auto erg = e_grid[j];
    out_wlh << std::scientific << std::setprecision(5)
            << erg << "\t";
    out_kd  << std::scientific << std::setprecision(5)
            << erg << "\t";
    
    for (int i = 0; i <  niso; ++i) {
      const auto& [Z,A] = isotopes[i];
      
      out_wlh << std::scientific << std::setprecision(5)
              << wlh_mean.cmpl_surf_V(Z,A,erg) << "\t";
      out_kd  << std::scientific << std::setprecision(5)
              << kdn_uq.cmpl_surf_V(Z,A,erg) << "\t";
    }
    out_wlh << "\n";
    out_kd  << "\n";
  }

  return 0;
};
