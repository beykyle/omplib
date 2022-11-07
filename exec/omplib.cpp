#include "potential/ch_params.hpp"
#include "potential/kd_params.hpp"
#include "potential/params.hpp"
#include "potential/wlh_params.hpp"

#include "util/constants.hpp"
#include "util/types.hpp"

#include "potential/potential.hpp"
#include "solver/rmatrix.hpp"


#include <iomanip>
#include <ios>
#include <iostream>
#include <fstream>

// this file is a loose collection of examples of things you can do with this library

// nice shorthand
constexpr auto n = omplib::Proj::neutron;
constexpr auto p = omplib::Proj::proton;

// l-wave scattering on 0+ ground state
void rmatrix_neutron_scatter(int Z, int A, int l) {
  
  using namespace omplib;
  using namespace omplib::constants;
  
  // constants
  constexpr real n_spin          = 1./2.;
  constexpr real projectile_mass = n_mass_amu; 
  const real target_mass         = A; 
  
  constexpr real threshold = 0.; // MeV
  constexpr real erg_cms   = 1.; // MeV
  constexpr real ch_radius = 12; // fm
  
  // on 0+ ground state
  constexpr int J2  = 1; // J2 == (2*J +1) == dimension of SU(2) representation
  constexpr auto pi = Parity::even;

  // build potentials and define scattering channel
  auto wlh_mean   = WLH21Params<n>();
  auto pwlh = OMP(Z, A, (2*l+1), (2*n_spin + 1), (2*(l+n_spin) + 1), wlh_mean);
  const auto ch   = Channel(threshold, erg_cms, ch_radius, 
                            projectile_mass, target_mass, Z, 0, l, J2, pi);
  
  // build radial grid to print the wavefunction
  constexpr auto r_grid_sz = 200;
  std::vector<real> r_grid(r_grid_sz,0.);
  for (int i = 0; i < r_grid_sz; ++i) {
    r_grid[i] = ch.radius * static_cast<real>(i) / (static_cast<real>(r_grid_sz) +1);
  }
  
  // run the R-Matrix solver with N basis functions
  constexpr unsigned int N = 10;

  // spin down
  // J = L+1/2
  // j2 = 2*J+1
  auto solver = RMatrixSolverSingleChannel<N>(ch, 
      [&pwlh, &ch](real r, real rp) -> cmpl {
        if ( r != rp ) return 0.;
        return pwlh.eval(r, ch.energy);
      } 
      );
  const auto [Rp, Sp, Tp, Kp, wvfxnp] = solver.solve();

  if ( l > 0 ) {
    // spin down
    // J = L-1/2
    // j2 = 2*J+1
    pwlh.set_proj_tot_am(2*(l-n_spin) + 1);
    solver = RMatrixSolverSingleChannel<N>(ch, 
        [&pwlh, &ch](real r, real rp) -> cmpl {
          if ( r != rp ) return 0.;
          return pwlh.eval(r, ch.energy);
        } 
      );
    const auto [Rm, Sm, Tm, Km, wvfxnm] = solver.solve();
  }

  //TODO tot, el, rxn xs, analyzing powers, differential el xs
}

void print_potential_vals() {
  
  using namespace omplib;
  auto kdn_uq     = KD03Params<n>::build_KDUQ();
  auto wlh_mean   = WLH21Params<n>();
  
  // let's look at potnetial values for various isotopes on an energy grid
  constexpr auto erg_min   = 0.01;
  constexpr auto erg_max   = 10.;
  constexpr auto e_range   = erg_max - erg_min;
  constexpr auto e_grid_sz = 500;
  auto e_grid = std::array<real,e_grid_sz> {};
  for (int i = 0; i < e_grid_sz; ++i) {
    e_grid[i] = erg_min 
      + e_range * static_cast<real>(i) / static_cast<real>(e_grid_sz);
  }

  // mass 144 isotopes
  using Isotope = std::pair<int,int>;
  constexpr auto niso = 6;
  constexpr auto isotopes =  std::array<Isotope,niso>{
    Isotope{58,144},Isotope{57,144},Isotope{56,144},
    Isotope{55,144},Isotope{54,144},Isotope{53,144}
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
    
    for (int i = 0; i < niso; ++i) {
      const auto& [Z,A] = isotopes[i];
      
      out_wlh << std::scientific << std::setprecision(5)
              << wlh_mean.cmpl_surf_V(Z,A,erg) << "\t";
      out_kd  << std::scientific << std::setprecision(5)
              << kdn_uq.cmpl_surf_V(Z,A,erg) << "\t";
    }
    out_wlh << "\n";
    out_kd  << "\n";
  }
}

int main(int argc, char** argv) {
  rmatrix_neutron_scatter(54,139,0); // s-wave scatter on 139-Xe
  print_potential_vals();
  return 0;
};
