import omplibpy as omp

# 0.1 Mev neutron incident on 48-Ca, using WLH microscopic optical potential
erg_cms_Mev = 0.1
A = 20
Z = 48
xs_tot, xs_rxn = omp.wlh_xs_n(0.1, A, Z)
