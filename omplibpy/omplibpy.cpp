#ifndef PYBINDINGS
#define PYBINDINGS

#include <iostream>

#include <pybind11/pybind11.h>

#include "solver/scatter.hpp"

using namespace omplib;

constexpr auto n = omplib::Proj::neutron;


std::pair<double,double> wlh_xs_n(double erg, int A, int Z) {
  const auto nuc = Isotope{Z,A,(real)A};
  const auto wlh = WLH21Params<n>();
  const auto pot = OMP(nuc, wlh); 

  const auto solver = NAScatter<n>(nuc); 

  auto p = [&pot](double r, double rp, const PData& d) -> cmpl { 
    return pot.eval_reduced(r, rp, d); 
  };

  const auto& [xs_tot, xs_rxn]  = solver.xs(erg, p);
  return {xs_tot, xs_rxn};
}




namespace py = pybind11;

PYBIND11_MODULE(omplibpy, m) {
  m.doc() = "omplibpy";
  m.def("wlh_xs_n" , &wlh_xs_n, "get WLH xs for isotope A,Z at given energy" );
}

#endif
