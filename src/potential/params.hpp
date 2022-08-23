#ifndef PARAMS_HPP
#define PARAMS_HPP

#include "util/constants.hpp"

#include "nlohmann/json.hpp"
using nlohmann::json;

namespace omplib {

/// @brief incident particle type
enum class Proj : bool {
  proton,
  neutron,
};
 
/// @tparam  incident neutron or proton 
template<Proj p>
/// @brief  Pure abstract base class for OM potential parameters
struct OMParams {
  constexpr static Proj projectile = p;

  // Woods-Saxon term radii 
  virtual double real_volu_r(int Z, int A, double erg) const = 0;
  virtual double cmpl_volu_r(int Z, int A, double erg) const = 0;
  virtual double real_surf_r(int Z, int A, double erg) const = 0;
  virtual double cmpl_surf_r(int Z, int A, double erg) const = 0;
  virtual double real_spin_r(int Z, int A, double erg) const = 0;
  virtual double cmpl_spin_r(int Z, int A, double erg) const = 0;
  
  // Woods-Saxon term diffusivity
  virtual double real_volu_a(int Z, int A, double erg) const = 0;
  virtual double cmpl_volu_a(int Z, int A, double erg) const = 0;
  virtual double real_surf_a(int Z, int A, double erg) const = 0;
  virtual double cmpl_surf_a(int Z, int A, double erg) const = 0;
  virtual double real_spin_a(int Z, int A, double erg) const = 0;
  virtual double cmpl_spin_a(int Z, int A, double erg) const = 0;

  // Woods-Saxon term depth
  virtual double real_volu_V(int Z, int A, double erg) const = 0;
  virtual double cmpl_volu_V(int Z, int A, double erg) const = 0;
  virtual double real_surf_V(int Z, int A, double erg) const = 0;
  virtual double cmpl_surf_V(int Z, int A, double erg) const = 0;
  virtual double real_spin_V(int Z, int A, double erg) const = 0;
  virtual double cmpl_spin_V(int Z, int A, double erg) const = 0;
};

template<>
struct OMParams<Proj::proton> {
  virtual double real_coul_r(int Z, int A, double erg) const = 0;
  virtual double real_coul_V(int Z, int A, double erg) const = 0;
};


/// @brief Phenomenological global OM potential parameterization from
/// A. Koning and J. Delaroche, 
/// Nuclear Physics A 713, 231 (2003), ISSN 0375-9474, 
/// URL https://www.sciencedirect.com/science/article/pii/S0375947402013210.
template <Proj projectile>
class KoningDelaroche03 : public OMParams<projectile> {
protected:
  // fermi energy
  double e_fermi_0, e_fermi_A;
  
  // central terms shape 
  double rv_0, rv_A, av_0, av_A;
  
  // complex central term shape
  double rd_0, rd_A, ad_0, ad_A;

  // spin orbit terms shape
  double rso_0, rso_A, aso_0;

  // real central depth
  double v1_0, v1_asym, v1_A, v2_0, v2_A, v3_0, v3_A, v4_0; 

  // complex central depth
  double w1_0, w1_A, w2_0, w2_A;
  
  // complex surface depth
  double d1_0, d1_asym, d2_0, d2_A, d2_A2, d2_A3, d3_0;

  // real spin orbit depth
  double vso1_0, vso1_A, vso2_0;
  
  // cmplex spin orbit depth
  double wso1_0, wso2_0;

  // structure and energy factors
  double Ef(int A) const;
  double asym(int Z, int A) const;

public:

  double real_volu_r(int Z, int A, double erg) const override;
  double cmpl_volu_r(int Z, int A, double erg) const override;
  double real_surf_r(int Z, int A, double erg) const override;
  double cmpl_surf_r(int Z, int A, double erg) const override;
  double real_spin_r(int Z, int A, double erg) const override;
  double cmpl_spin_r(int Z, int A, double erg) const override;
  
  double real_volu_a(int Z, int A, double erg) const override;
  double cmpl_volu_a(int Z, int A, double erg) const override;
  double real_surf_a(int Z, int A, double erg) const override;
  double cmpl_surf_a(int Z, int A, double erg) const override;
  double real_spin_a(int Z, int A, double erg) const override;
  double cmpl_spin_a(int Z, int A, double erg) const override;

  double real_volu_V(int Z, int A, double erg) const override;
  double cmpl_volu_V(int Z, int A, double erg) const override;
  double real_surf_V(int Z, int A, double erg) const override;
  double cmpl_surf_V(int Z, int A, double erg) const override;
  double real_spin_V(int Z, int A, double erg) const override;
  double cmpl_spin_V(int Z, int A, double erg) const override;

  KoningDelaroche03<projectile>( const KoningDelaroche03<projectile>& rhs) = default;

  // @brief  Construct using params supplied in a json file
  KoningDelaroche03<projectile>(json param_file);

  // @brief Construct using the default KD03 params
  KoningDelaroche03<projectile>();
};


template <Proj p>
/// @brief constructs a KoningDelaroche03<p> with params refit w/ MCMC; from
/// Pruitt, C. D. et al, 
/// “Uncertainty-Quantified Phenomenological Optical Potentials 
/// for Single-Nucleon Scattering”, 
/// LLNL release number LLNL-JRNL-835671-DRAFT (to be published).
class KDUQ : public KoningDelaroche03<p> { 
public: 
    KDUQ(); 
};

template<>
class KoningDelaroche03<Proj::proton> 
  : public KoningDelaroche03<Proj::neutron> 
  , public OMParams<Proj::proton> {
protected: 
  double rc_0, rc_A, rc_A2;

public:
  constexpr static Proj projectile = Proj::proton;
  double real_coul_r(int Z, int A, double erg) const final;
  double real_coul_V(int Z, int A, double erg) const final;
  
  KoningDelaroche03( 
      const KoningDelaroche03<Proj::proton>& rhs) = default;
  KoningDelaroche03(json param_file);
  KoningDelaroche03();
};

// constructors
template<>
KoningDelaroche03<Proj::neutron>::KoningDelaroche03() 
  : OMParams<Proj::neutron>()
  , e_fermi_0(-11.2814), e_fermi_A(0.02646)
  , v1_0(5.93E1)  , v1_asym(2.10E1), v1_A(2.40E-2)
  , v2_0(7.228E-3), v2_A(1.48e-6)  
  , v3_0(1.994E-5), v3_A(2E-8)
  , v4_0(7.00E-9)
  , w1_0(12.195)  , w1_A(0.0167)    , w2_0(73.55)    , w2_A(0.0795)
  , d1_0(16.0)    , d1_asym(16.0)
  , d2_0(0.0180)  , d2_A(0.003802)  , d2_A2(8)       , d2_A3(156)
  , d3_0(1.15E1)
  , vso1_0(5.922) , vso1_A(0.0030)
  , vso2_0(0.0040)
  , wso1_0(-3.1)
  , wso2_0(160)
  , rv_0(1.3039E0), av_0(6.778E-1)  , av_A(1.487E-4)
  , rd_0(1.3424)  , rd_A(0.01585)   , ad_0(0.5446)   , ad_A(1.656E-4)
  , rso_0(1.1854) , rso_A(0.647)    , aso_0(0.59)
   {};

KoningDelaroche03<Proj::proton>::KoningDelaroche03() 
  : OMParams<Proj::proton>()
  , KoningDelaroche03<Proj::neutron>()
  , rc_0(1.2E0), rc_A(6.97E-1), rc_A2(1.3E1)
 {
   // params which are different for protons
   e_fermi_0 = -8.4075;
   e_fermi_A = 1.01378;
   v2_0 = 7.067E-3;
   v2_A = 4.23E-6;
   v3_0 = 1.729E-5;
   v3_A = 1.136E-8;
   w1_0 = 1.4667E1;
   w1_A = 9.629E-3;
   av_0 = 5.19E-1;
   av_A = 5.21E-4;
 };

template<>
KDUQ<Proj::neutron>::KDUQ() 
  : KoningDelaroche03<Proj::neutron>() {
  v1_0    = 5.86E1;  
  v1_asym = 1.34E1;  
  v1_A    = 2.61E-2; 
  v4_0    = -4.3E-9; 
  rv_0    = 1.27E0;  
  rv_A    = 3.61E-1; 
  av_0    = 6.89E-1; 
  av_A    = -0.42E-4;
  vso1_0  = 5.99E0;  
  vso1_A  = 1.95E-3; 
  vso2_0  = 4.75E-3; 
  rso_0   = 1.21E0;  
  rso_A   = 7.35E-1; 
  aso_0   = 6.00E-1; 
  wso1_0  = -3.79E0; 
  wso2_0  = 2.19E2;  
  w2_0    = 10.29E1; 
  w2_A    = 2.43E-2; 
  d1_0    = 1.67E1;  
  d1_asym = 1.11E1;  
  d2_0    = 2.34E-2; 
  d2_A    = 3.73E-3; 
  d2_A2   = 8.57E0;  
  d2_A3   = 2.51E2;  
  d3_0    = 1.38E1;  
  rd_0    = 1.35E0;  
  rd_A    = 1.75E-2; 
                      
  // different for incident protons nand neutrons
  w1_0    = 2.09E1;   
  w1_A    = 0.61E-2;  
  v2_0    = 6.35E-3;
  v2_A    = 1.82E-6;
  v3_0    = 1.08E-5;
  v3_A    = 1.45E-8;
  ad_0    = 5.43E-1;
  ad_A    = -2.14E-4;
};

template<>
KDUQ<Proj::proton>::KDUQ() 
  : KoningDelaroche03<Proj::proton>() {
  v1_0    = 5.86E1;
  v1_asym = 1.34E1;
  v1_A    = 2.61E-2;
  v4_0    = -4.3E-9;
  rv_0    = 1.27E0;
  rv_A    = 3.61E-1;
  av_0    = 6.89E-1;
  av_A    = -0.42E-4;
  vso1_0  = 5.99E0;
  vso1_A  = 1.95E-3;
  vso2_0  = 4.75E-3;
  rso_0   = 1.21E0;
  rso_A   = 7.35E-1;
  aso_0   = 6.00E-1;
  wso1_0  = -3.79E0;
  wso2_0  = 2.19E2;
  w2_0    = 10.29E1;
  w2_A    = 2.43E-2;
  d1_0    = 1.67E1;
  d1_asym = 1.11E1;
  d2_0    = 2.34E-2;
  d2_A    = 3.73E-3;
  d2_A2   = 8.57E0;
  d2_A3   = 2.51E2;
  d3_0    = 1.38E1;
  rd_0    = 1.35E0;
  rd_A    = 1.75E-2;
  rc_0    = 1.19E0;
  rc_A    = 6.72E-1;
  rc_A2   = 1.3E1;

  // different for incident protons nand neutrons
  v2_0    = 6.76E-3;
  v2_A    = 2.91E-6;
  v3_0    = 6.76E-3;
  v3_A    = 1.43E-8;
  w1_0    = 1.86E1;
  w1_A    = 32.5E-3;
  ad_0    = 5.08E-1;
  ad_A    = 14.10E-4;
};

// potential terms
template<Proj p>
double KoningDelaroche03<p>::Ef(int A) const {
  return e_fermi_0 + e_fermi_A * static_cast<double>(A);
};

template<Proj p>
double KoningDelaroche03<p>::asym(int Z, int A) const {
  const double a = static_cast<double>(A);
  const double z = static_cast<double>(Z);
  const double n = a - z;
  return (n - z)/a;
}

template<Proj p>
double KoningDelaroche03<p>::real_volu_r(
    int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return rv_0 - rv_A / pow(A , -1./3.);
};

template<Proj p>
double KoningDelaroche03<p>::cmpl_volu_r(
    int Z, int A, double erg) const {
  return real_volu_r(Z,A,erg); // complex and real volume terms share geometry
};

template<Proj p>
double KoningDelaroche03<p>::real_surf_r(
    int Z, int A, double erg) const {
  return 0; // no real surf term in KD03
};

template<Proj p>
double KoningDelaroche03<p>::cmpl_surf_r(
    int Z, int A, double erg) const {
  return rd_0 - rd_A * pow(static_cast<double>(A), -1./3.);;
};

template<Proj p>
double KoningDelaroche03<p>::real_spin_r(
    int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return rso_0 - rso_A * pow(a, -1./3.); 
};

template<Proj p>
double KoningDelaroche03<p>::cmpl_spin_r(
    int Z, int A, double erg) const {
  return real_spin_r(Z,A,erg);// complex and real SO terms share geometry
};

template<Proj p>
double KoningDelaroche03<p>::real_volu_a(
    int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return av_0 - av_A / a;
};

template<Proj p>
double KoningDelaroche03<p>::cmpl_volu_a(
    int Z, int A, double erg) const {
  return real_volu_a(Z,A,erg); // complex and real volume terms share geometry
};

template<Proj p>
double KoningDelaroche03<p>::real_surf_a(
    int Z, int A, double erg) const {
  return 0; // no real surface term in KD
};

template<Proj p>
double KoningDelaroche03<p>::cmpl_surf_a(
    int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  if constexpr (p == Proj::neutron)
    return ad_0 - ad_A * a;
  else if constexpr(p == Proj::proton)
    return ad_0 + ad_A * a;
};

template<Proj p>
double KoningDelaroche03<p>::real_spin_a(
    int Z, int A, double erg) const {
  return aso_0;
};

template<Proj p>
double KoningDelaroche03<p>::cmpl_spin_a(
    int Z, int A, double erg) const {
  return real_spin_a(Z,A,erg); // complex and real SO term share geometry
};

template<Proj p>
double KoningDelaroche03<p>::real_volu_V(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  const double alpha = asym(Z,A);
  const double dE    = erg - Ef(A); 

  if constexpr (p == Proj::neutron) {
    const double v1 = v1_0 - v1_A * a - v1_asym * alpha;
    const double v2 = v2_0 - v2_A * a;
    const double v3 = v3_0 - v3_A * a;
    const double v4 = v4_0;

    return v1 * (1 - v2 * dE + v3 * dE * dE - v4 * dE * dE * dE );
  }
  else if constexpr (p == Proj::proton) {
    const double v1 = v1_0 - v1_A * a + v1_asym * alpha;
    const double v2 = v2_0 + v2_A * a;
    const double v3 = v3_0 + v3_A * a;
    const double v4 = v4_0;
    const double vc = KoningDelaroche03<Proj::proton>::real_coul_V(Z,A,erg);

    return v1 * (1 - v2 * dE + v3 * dE * dE - v4 * dE * dE * dE )
      + vc * v1 * ( v2 - 2 * v3 * dE + 3 * v4 * dE * dE);
    }
}

template<Proj p>
double KoningDelaroche03<p>::cmpl_volu_V(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  const double alpha = asym(Z,A);
  const double dE    = erg - Ef(A); 
  
  const double w1 = w1_0 + w1_A * a;
  const double w2 = w2_0 + w2_A * a;

  return w1 * dE * dE / ( dE * dE + w2 * w2);
}

template<Proj p>
double KoningDelaroche03<p>::real_surf_V(int Z, int A, double erg) const {
  return 0; // no real surf term in KD03
}

template<Proj p>
double KoningDelaroche03<p>::cmpl_surf_V(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  const double alpha = asym(Z,A);
  const double dE    = erg - Ef(A); 

  const double d2 = d2_0 + d2_A / (1 + exp( (a - d2_A3) / d2_A2 ) );
  const double d3 = d3_0;

  auto WD = [dE](double d1, double d2, double d3){
    return d1 * dE * dE / (dE * dE + d3 * d3) * exp( - d2 * dE);
  };

  if constexpr (p == Proj::neutron) {
    const double d1 = d1_0 - d1_asym * alpha;
    return WD(d1,d2,d3);
  }
  else if constexpr (p == Proj::proton) {
    const double d1 = d1_0 + d1_asym * alpha;
    return WD(d1,d2,d3);
  }
}

template<Proj p>
double KoningDelaroche03<p>::real_spin_V(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  const double dE    = erg - Ef(A); 

  const double vso1 = vso1_0 + vso1_A * a;
  const double vso2 = vso2_0;

  return vso1 * exp( - vso2 * dE);
}

template<Proj p>
double KoningDelaroche03<p>::cmpl_spin_V(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  const double dE    = erg - Ef(A); 

  const double wso1 = wso1_0;
  const double wso2 = wso2_0;

  return wso1 * dE * dE / ( dE * dE + wso2 * wso2 );
}

double omplib::KoningDelaroche03<Proj::proton>::real_coul_V(
    int Z, int A, double erg) const {
  const double z = static_cast<double>(Z);
  const double a = static_cast<double>(A);
  return 6. * z * e_sqr / ( 5. * real_coul_r(z,a,erg) * pow(a,1./3.));
}

double omplib::KoningDelaroche03<Proj::proton>::real_coul_r(
    int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return rc_0 + rc_A * pow(a, -1./3.) + rc_A2 * pow(a, -5./3.);
}



/// @brief  Phenomenological global OM potential parameterization from:
/// R. Varner, W. Thompson, T. McAbee, E. Ludwig, and T. Clegg, 
/// Physics Reports 201, 57 (1991), ISSN 0370-1573, 
/// URL https://www.sciencedirect.com/science/article/pii/037015739190039O
template <Proj projectile>
class ChapelHill89 : public OMParams<projectile> {
public:
  
  double real_volu_r(int Z, int A, double erg) const final;
  double cmpl_volu_r(int Z, int A, double erg) const final;
  double real_surf_r(int Z, int A, double erg) const final;
  double cmpl_surf_r(int Z, int A, double erg) const final;
  double real_spin_r(int Z, int A, double erg) const final;
  double cmpl_spin_r(int Z, int A, double erg) const final;
  
  double real_volu_a(int Z, int A, double erg) const final;
  double cmpl_volu_a(int Z, int A, double erg) const final;
  double real_surf_a(int Z, int A, double erg) const final;
  double cmpl_surf_a(int Z, int A, double erg) const final;
  double real_spin_a(int Z, int A, double erg) const final;
  double cmpl_spin_a(int Z, int A, double erg) const final;

  double real_volu_V(int Z, int A, double erg) const final;
  double cmpl_volu_V(int Z, int A, double erg) const final;
  double real_surf_V(int Z, int A, double erg) const final;
  double cmpl_surf_V(int Z, int A, double erg) const final;
  double real_spin_V(int Z, int A, double erg) const final;

  ChapelHill89<projectile>(json param_file);
  ChapelHill89<projectile>();
};

template <Proj p>
/// @brief constructs a ChapelHill89<p> with params refit w/ MCMC; from
/// Pruitt, C. D. et al, 
/// “Uncertainty-Quantified Phenomenological Optical Potentials 
/// for Single-Nucleon Scattering”, 
/// LLNL release number LLNL-JRNL-835671-DRAFT (to be published).
class CHUQ : public ChapelHill89<p> { public: CHUQ(); };

template <>
class CHUQ<Proj::proton> { public: CHUQ(); };


/// @brief  Microscopic global OM potential parameterization using XEFT, from
/// T. R. Whitehead, Y. Lim, and J. W. Holt, 
/// Phys. Rev. Lett. 127, 182502 (2021), 
/// URL https://link.aps.org/doi/10.1103/PhysRevLett.127.182502.
template <Proj projectile>
class WLH21 : public OMParams<projectile> {
public:
  
  double real_volu_r(int Z, int A, double erg) const final;
  double cmpl_volu_r(int Z, int A, double erg) const final;
  double real_surf_r(int Z, int A, double erg) const final;
  double cmpl_surf_r(int Z, int A, double erg) const final;
  double real_spin_r(int Z, int A, double erg) const final;
  double cmpl_spin_r(int Z, int A, double erg) const final;
  
  double real_volu_a(int Z, int A, double erg) const final;
  double cmpl_volu_a(int Z, int A, double erg) const final;
  double real_surf_a(int Z, int A, double erg) const final;
  double cmpl_surf_a(int Z, int A, double erg) const final;
  double real_spin_a(int Z, int A, double erg) const final;
  double cmpl_spin_a(int Z, int A, double erg) const final;

  double real_volu_V(int Z, int A, double erg) const final;
  double cmpl_volu_V(int Z, int A, double erg) const final;
  double real_surf_V(int Z, int A, double erg) const final;
  double cmpl_surf_V(int Z, int A, double erg) const final;
  double real_spin_V(int Z, int A, double erg) const final;

  WLH21<projectile>(json param_file);
  WLH21<projectile>();
};
};


#endif
