#ifndef KD_PARAMS_HEADER
#define KD_PARAMS_HEADER

#include "potential/params_base.hpp"

namespace omplib {

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

  using OMParams<projectile>::asym;

public:

  double real_cent_r(int Z, int A, double erg) const override;
  double cmpl_cent_r(int Z, int A, double erg) const override;
  double cmpl_surf_r(int Z, int A, double erg) const override;
  double real_spin_r(int Z, int A, double erg) const override;
  double cmpl_spin_r(int Z, int A, double erg) const override;
  
  double real_cent_a(int Z, int A, double erg) const override;
  double cmpl_cent_a(int Z, int A, double erg) const override;
  double cmpl_surf_a(int Z, int A, double erg) const override;
  double real_spin_a(int Z, int A, double erg) const override;
  double cmpl_spin_a(int Z, int A, double erg) const override;

  double real_cent_V(int Z, int A, double erg) const override;
  double cmpl_cent_V(int Z, int A, double erg) const override;
  double cmpl_surf_V(int Z, int A, double erg) const override;
  double real_spin_V(int Z, int A, double erg) const override;
  double cmpl_spin_V(int Z, int A, double erg) const override;
  
  // KD03 does not have a real surface term
  double real_surf_a(int Z, int A, double erg) const override { return 0; }
  double real_surf_V(int Z, int A, double erg) const override { return 0; }
  double real_surf_r(int Z, int A, double erg) const override { return 0; }

  KoningDelaroche03( const KoningDelaroche03<projectile>& rhs) = default;

  // @brief Construct using the default KD03 params
  KoningDelaroche03();
  // @brief  Construct using params supplied in a json file
  KoningDelaroche03(json p): 
    // same for n's and p's 
    v1_0(   p["KDHartreeFock"]["V1_0"]     ) , v1_asym(  p["KDHartreeFock"]["V1_asymm"]   ),
    v1_A(   p["KDHartreeFock"]["V1_A"]     ) , v4_0(     p["KDHartreeFock"]["V4_0"]       ),
    rv_0(   p["KDHartreeFock"]["r_0"]      ) , rv_A(     p["KDHartreeFock"]["r_A"]        ),
    av_0(   p["KDHartreeFock"]["a_0"]      ) , av_A(     p["KDHartreeFock"]["a_A"]        ),
    w2_0(   p["KDImagVolume"]["W2_0"]      ) , w2_A(     p["KDImagVolume"]["W2_A"]        ),
    d1_0(   p["KDImagSurface"]["D1_0"]     ) , d1_asym(  p["KDImagSurface"]["D1_asymm"]   ),
    d2_0(   p["KDImagSurface"]["D2_0"]     ) , d2_A(     p["KDImagSurface"]["D2_A"]       ),
    d2_A2(  p["KDImagSurface"]["D2_A2"]    ) , d2_A3(    p["KDImagSurface"]["D2_A3"]      ),
    d3_0(   p["KDImagSurface"]["D3_0"]     ) , rd_0(     p["KDImagSurface"]["r_0"]        ),
    rd_A(   p["KDImagSurface"]["r_A"]      ) , vso1_0(   p["KDRealSpinOrbit"]["V1_0"]     ),
    vso1_A( p["KDRealSpinOrbit"]["V1_A"]   ) , rso_0(    p["KDRealSpinOrbit"]["r_0"]      ),
    vso2_0( p["KDRealSpinOrbit"]["V2_0"]   ) , aso_0(    p["KDRealSpinOrbit"]["a_0"]      ), 
    rso_A(  p["KDRealSpinOrbit"]["r_A"]    ) , wso2_0(   p["KDImagSpinOrbit"]["W2_0"]     ), 
    wso1_0( p["KDImagSpinOrbit"]["W1_0"]   ) 

    // different for neutrons and protons
    {
      if constexpr (projectile == Proj::neutron) {
        e_fermi_0 = -11.2814;
        e_fermi_A =  0.02646;
        ad_0 =  p["KDImagSurface"]["a_0_n" ];
        ad_A =  p["KDImagSurface"]["a_A_n" ]; 
        v2_0 =  p["KDHartreeFock"]["V2_0_n"];
        v2_A =  p["KDHartreeFock"]["V2_A_n"]; 
        v3_0 =  p["KDHartreeFock"]["V3_0_n"];
        v3_A =  p["KDHartreeFock"]["V3_A_n"]; 
        w1_0 =  p["KDImagVolume" ]["W1_0_n"]; 
        w1_A =  p["KDImagVolume" ]["W1_A_n"];
      }
      else if constexpr (projectile == Proj::proton) {
        e_fermi_0 = -8.4075;
        e_fermi_A = 0.01378;
        ad_0 =  p["KDImagSurface"]["a_0_p" ];
        ad_A =  p["KDImagSurface"]["a_A_p" ]; 
        v2_0 =  p["KDHartreeFock"]["V2_0_p"];
        v2_A =  p["KDHartreeFock"]["V2_A_p"]; 
        v3_0 =  p["KDHartreeFock"]["V3_0_p"];
        v3_A =  p["KDHartreeFock"]["V3_A_p"]; 
        w1_0 =  p["KDImagVolume" ]["W1_0_p"]; 
        w1_A =  p["KDImagVolume" ]["W1_A_p"];
      }
    }
    
  
  /// @brief constructs a KoningDelaroche03<p> with params refit w/ MCMC; from
  /// @brief constructs a KoningDelaroche03\<p\> with params refit w/ MCMC; from
  /// Pruitt, C. D. et al, 
  /// “Uncertainty-Quantified Phenomenological Optical Potentials 
  /// for Single-Nucleon Scattering”, 
  static KoningDelaroche03<projectile> build_KDUQ();
};

template<>
class KoningDelaroche03<Proj::proton> : 
  public KoningDelaroche03<Proj::neutron> , OMParams<Proj::proton> {
protected: 
  double rc_0, rc_A, rc_A2;

public:
  constexpr static Proj projectile = Proj::proton;
  double real_coul_r(int Z, int A, double erg) const final;
  
  KoningDelaroche03( 
      const KoningDelaroche03<Proj::proton>& rhs) = default;
  KoningDelaroche03();
  KoningDelaroche03(json p);
  static KoningDelaroche03<Proj::proton> build_KDUQ();
};


// potential terms
template<Proj proj>
double KoningDelaroche03<proj>::Ef(int A) const {
  return e_fermi_0 + e_fermi_A * static_cast<double>(A);
};

template<Proj proj>
double KoningDelaroche03<proj>::real_cent_r(
    int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return rv_0 * pow(a,1./3.)- rv_A ;
};

template<Proj proj>
double KoningDelaroche03<proj>::cmpl_cent_r(
    int Z, int A, double erg) const {
  return real_cent_r(Z,A,erg); // complex and real centme terms share geometry
};

template<Proj proj>
double KoningDelaroche03<proj>::cmpl_surf_r(
    int Z, int A, double erg) const {
  const auto a3 = pow(static_cast<double>(A), 1./3.);
  return rd_0 *a3- rd_A * a3*a3;
};

template<Proj proj>
double KoningDelaroche03<proj>::real_spin_r(
    int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return rso_0 * pow(a, 1./3.) - rso_A ; 
};

template<Proj proj>
double KoningDelaroche03<proj>::cmpl_spin_r(
    int Z, int A, double erg) const {
  return real_spin_r(Z,A,erg);// complex and real SO terms share geometry
};

template<Proj proj>
double KoningDelaroche03<proj>::real_cent_a(
    int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  return av_0 - av_A * a;
};

template<Proj proj>
double KoningDelaroche03<proj>::cmpl_cent_a(
    int Z, int A, double erg) const {
  return real_cent_a(Z,A,erg); // complex and real centme terms share geometry
};

template<Proj proj>
double KoningDelaroche03<proj>::cmpl_surf_a(
    int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  if constexpr (proj == Proj::neutron)
    return ad_0 - ad_A * a;
  else if constexpr(proj == Proj::proton)
    return ad_0 + ad_A * a;
};

template<Proj proj>
double KoningDelaroche03<proj>::real_spin_a(
    int Z, int A, double erg) const {
  return aso_0;
};

template<Proj proj>
double KoningDelaroche03<proj>::cmpl_spin_a(
    int Z, int A, double erg) const {
  return real_spin_a(Z,A,erg); // complex and real SO term share geometry
};

template<Proj proj>
double KoningDelaroche03<proj>::real_cent_V(int Z, int A, double erg) const {
  const double z = static_cast<double>(Z);
  const double a = static_cast<double>(A);
  const double alpha = asym(Z,A);
  const double dE    = erg - Ef(A); 
  const double v1 = v1_0 - v1_A * a + v1_asym * alpha;
  const double v4 = v4_0;

  if constexpr (proj == Proj::neutron) {
    const double v2 = v2_0 - v2_A * a;
    const double v3 = v3_0 - v3_A * a;

    return v1 * (1. - v2 * dE + v3 * dE * dE - v4 * dE * dE * dE );
  }
  else if constexpr (proj == Proj::proton) {
    const double v2 = v2_0 + v2_A * a;
    const double v3 = v3_0 + v3_A * a;
    const double rc = KoningDelaroche03<Proj::proton>::real_coul_r(Z,A,erg);
    const double vc = 6 * z * e_sqr / (5 * rc * pow(a,1./3.));

    return v1 * (1. - v2 * dE + v3 * dE * dE - v4 * dE * dE * dE )
      + vc * v1 * ( v2 - 2. * v3 * dE + 3. * v4 * dE * dE);
  }
}

template<Proj proj>
double KoningDelaroche03<proj>::cmpl_cent_V(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  const double dE    = erg - Ef(A); 
  
  const double w1 = w1_0 + w1_A * a;
  const double w2 = w2_0 + w2_A * a;

  return w1 * dE * dE / ( dE * dE + w2 * w2);
}

template<Proj proj>
double KoningDelaroche03<proj>::cmpl_surf_V(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  const double alpha = asym(Z,A);
  const double dE    = erg - Ef(A); 

  const double d1 = d1_0 + d1_asym * alpha;
  const double d2 = d2_0 + d2_A / (1 + exp( (a - d2_A3) / d2_A2 ) );
  const double d3 = d3_0;

  return d1 * dE * dE / (dE * dE + d3 * d3) * exp( - d2 * dE);
}

template<Proj proj>
double KoningDelaroche03<proj>::real_spin_V(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  const double dE    = erg - Ef(A); 

  const double vso1 = vso1_0 + vso1_A * a;
  const double vso2 = vso2_0;

  return vso1 * exp( - vso2 * dE);
}

template<Proj proj>
double KoningDelaroche03<proj>::cmpl_spin_V(int Z, int A, double erg) const {
  const double a = static_cast<double>(A);
  const double dE    = erg - Ef(A); 

  const double wso1 = wso1_0;
  const double wso2 = wso2_0;

  return wso1 * dE * dE / ( dE * dE + wso2 * wso2 );
}

}
#endif
