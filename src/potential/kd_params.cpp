#include "potential/kd_params.hpp"

using namespace omplib;

template<>
KD03Params<Proj::neutron>::KD03Params() 
  : OMParams<Proj::neutron>(),
  
  e_fermi_0(-11.2814), e_fermi_A(0.02646),

  rv_0(1.3039)       , rv_A(0.4054)      , av_0(6.778e-1) ,
  av_A(1.487e-4)     , rd_0(1.3424)      , rd_A(0.01585)  , 
  ad_0(0.5446)       , ad_A(1.656e-4)    , rso_0(1.1854)  ,
  rso_A(0.647)       , aso_0(0.59)       ,
  
  v1_0(5.93e1)       , v1_asym(2.10e1)   , v1_A(0.024)    ,
  v2_0(7.228e-3)     , v2_A(1.48e-6)     , v3_0(1.994e-5) ,
  v3_A(2.0e-8)       , v4_0(7.00e-9)     , w1_0(12.195)   ,  
  w1_A(0.0167)       , w2_0(73.55)       , w2_A(0.0795)   ,
  d1_0(16.0)         , d1_asym(16.0)     , d2_0(0.0180)   ,
  d2_A(0.003802)     , d2_A2(8)          , d2_A3(156)     ,       
  d3_0(1.15e1)       , vso1_0(5.922)     , vso1_A(0.0030) , 
  vso2_0(0.0040)     , wso1_0(-3.1)      , wso2_0(160)    
{};

KD03Params<Proj::proton>::KD03Params(json p) 
  : KD03Params<Proj::neutron>(p)
  , OMParams<Proj::proton>()
  , rc_0(  p["KDCoulomb_r_C_0"] )
  , rc_A(  p["KDCoulomb_r_C_A"] )
  , rc_A2( p["KDCoulomb_r_C_A2"])
{};

KD03Params<Proj::proton>::KD03Params() 
  : KD03Params<Proj::neutron>()
  , OMParams<Proj::proton>()
  , rc_0(1.2E0), rc_A(6.97E-1), rc_A2(1.3E1)
{
  // params which are different for protons
  e_fermi_0 = -8.4075;
  e_fermi_A = 0.01378;
  v2_0 = 7.067E-3;
  v2_A = 4.23E-6;
  v3_0 = 1.729E-5;
  v3_A = 1.136E-8;
  w1_0 = 1.4667E1;
  w1_A = 9.629E-3;
  av_0 = 5.19E-1;
  av_A = 5.21E-4;
};

real omplib::KD03Params<Proj::proton>::real_coul_r(
    int, int A, real) const {
  const real a = static_cast<real>(A);
  return rc_0 + rc_A * pow(a, -1./3.) + rc_A2 * pow(a, -5./3.);
}

template<>
KD03Params<Proj::neutron> 
KD03Params<Proj::neutron>::build_KDUQ() {
  KD03Params<Proj::neutron> p{};
  p.v1_0    = 5.86E1;  
  p.v1_asym = 1.34E1;  
  p.v1_A    = 2.61E-2; 
  p.v4_0    = -4.3E-9; 
  p.rv_0    = 1.27E0;  
  p.rv_A    = 3.61E-1; 
  p.av_0    = 6.89E-1; 
  p.av_A    = -0.42E-4;
  p.vso1_0  = 5.99E0;  
  p.vso1_A  = 1.95E-3; 
  p.vso2_0  = 4.75E-3; 
  p.rso_0   = 1.21E0;  
  p.rso_A   = 7.35E-1; 
  p.aso_0   = 6.00E-1; 
  p.wso1_0  = -3.79E0; 
  p.wso2_0  = 2.19E2;  
  p.w2_0    = 10.29E1; 
  p.w2_A    = 2.43E-2; 
  p.d1_0    = 1.67E1;  
  p.d1_asym = 1.11E1;  
  p.d2_0    = 2.34E-2; 
  p.d2_A    = 3.73E-3; 
  p.d2_A2   = 8.57E0;  
  p.d2_A3   = 2.51E2;  
  p.d3_0    = 1.38E1;  
  p.rd_0    = 1.35E0;  
  p.rd_A    = 1.75E-2; 
                      
  // different for incident protons and neutrons
  p.w1_0    = 2.09E1;   
  p.w1_A    = 0.61E-2;  
  p.v2_0    = 6.35E-3;
  p.v2_A    = 1.82E-6;
  p.v3_0    = 1.08E-5;
  p.v3_A    = 1.45E-8;
  p.ad_0    = 5.43E-1;
  p.ad_A    = -2.14E-4;
  
  p.e_fermi_0 = -11.2815;
  p.e_fermi_A = 0.02646;

  return p;
};

KD03Params<Proj::proton> 
KD03Params<Proj::proton>::build_KDUQ() {
  
  KD03Params<Proj::proton> p;
  p.v1_0    = 5.86E1;
  p.v1_asym = 1.34E1;
  p.v1_A    = 2.61E-2;
  p.v4_0    = -4.3E-9;
  p.rv_0    = 1.27E0;
  p.rv_A    = 3.61E-1;
  p.av_0    = 6.89E-1;
  p.av_A    = -0.42E-4;
  p.vso1_0  = 5.99E0;
  p.vso1_A  = 1.95E-3;
  p.vso2_0  = 4.75E-3;
  p.rso_0   = 1.21E0;
  p.rso_A   = 7.35E-1;
  p.aso_0   = 6.00E-1;
  p.wso1_0  = -3.79E0;
  p.wso2_0  = 2.19E2;
  p.w2_0    = 10.29E1;
  p.w2_A    = 2.43E-2;
  p.d1_0    = 1.67E1;
  p.d1_asym = 1.11E1;
  p.d2_0    = 2.34E-2;
  p.d2_A    = 3.73E-3;
  p.d2_A2   = 8.57E0;
  p.d2_A3   = 2.51E2;
  p.d3_0    = 1.38E1;
  p.rd_0    = 1.35E0;
  p.rd_A    = 1.75E-2;

  // different for incident protons and neutrons
  p.v2_0    = 6.76E-3;
  p.v2_A    = 2.91E-6;
  p.v3_0    = 6.76E-3;
  p.v3_A    = 1.43E-8;
  p.w1_0    = 1.86E1;
  p.w1_A    = 32.5E-3;
  p.ad_0    = 5.08E-1;
  p.ad_A    = 14.10E-4;
  
  p.rc_0  = 1.19E0;
  p.rc_A  = 6.72E-1;
  p.rc_A2 = 1.3E1;

  p.e_fermi_0 = -8.4075;
  p.e_fermi_A = 0.01378;
  
  return p;
};

