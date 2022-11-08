#ifndef GAUS_LEGENDRE_WEIGHTS
#define GAUS_LEGENDRE_WEIGHTS

#include <array>

#include "util/types.hpp"

namespace omplib {

/// @brief Number of basis functions
template<unsigned int N>
/// @brief  weights and abscissacissa of basis polynomials for use in 
/// Gauss-Legendre quadrature integration for x on shifted domain [0,1]
/// https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature
struct GaussLegendre {
  std::array<real,N> weights;
  std::array<real,N> abscissa;

  constexpr GaussLegendre() {
    static_assert(N >= 2);
    static_assert(N <= 20);

    // hard coded values for typical domain 
    // x on [-1,1]
    if constexpr (N == 2){
      abscissa = {
        -0.5773502691896257645092,
        0.5773502691896257645092
      };
      weights = {1,1};
    }
    if constexpr (N == 3){
      abscissa = {
        -0.7745966692414833770359,
        0,
        0.7745966692414833770359
      };
      weights = {
        0.5555555555555555555556, 
        0.8888888888888888888889, 
        0.555555555555555555556
      };
    }
    if constexpr (N == 4){
      abscissa = {
        -0.861136311594052575224
        -0.3399810435848562648027,
        0.3399810435848562648027,
        0.861136311594052575224,
      };
      weights = {
        0.3478548451374538573731,
        0.6521451548625461426269,
        0.6521451548625461426269,
        0.3478548451374538573731,
      };
    }
    if constexpr (N == 5){
      abscissa = {
        -0.9061798459386639927976,
        -0.5384693101056830910363,
        0,
        0.5384693101056830910363,
        0.9061798459386639927976
      };
      weights = {
        0.2369268850561890875143,
        0.4786286704993664680413,
        0.5688888888888888888889,
        0.4786286704993664680413,
        0.2369268850561890875143
      };
    }
    if constexpr (N == 6){
      abscissa = {
        -0.9324695142031520278123,
        -0.661209386466264513661,
        -0.2386191860831969086305,
        0.238619186083196908631,
        0.661209386466264513661,
        0.9324695142031520278123
      };
      weights = {
        0.1713244923791703450403,
        0.3607615730481386075698,
        0.4679139345726910473899,
        0.46791393457269104739,
        0.3607615730481386075698,
        0.1713244923791703450403
      };
    }
    if constexpr (N == 7){
      abscissa = {
        -0.9491079123427585245262,
        -0.7415311855993944398639,
        -0.4058451513773971669066,
        0,
        0.4058451513773971669066,
        0.7415311855993944398639,
        0.9491079123427585245262
      };
      weights = {
        0.1294849661688696932706,
        0.2797053914892766679015,
        0.38183005050511894495,
        0.417959183673469387755,
        0.38183005050511894495,
        0.279705391489276667901,
        0.129484966168869693271
      };
    }
    if constexpr (N == 8){
      abscissa = {
        -0.9602898564975362316836,
        -0.7966664774136267395916,
        -0.5255324099163289858177,
        -0.1834346424956498049395,
        0.1834346424956498049395,
        0.5255324099163289858177,
        0.7966664774136267395916,
        0.9602898564975362316836
      };
      weights = {
        0.1012285362903762591525,
        0.2223810344533744705444,
        0.313706645877887287338,
        0.3626837833783619829652,
        0.3626837833783619829652,
        0.313706645877887287338,
        0.222381034453374470544,
        0.1012285362903762591525
      };
    }
    if constexpr (N == 9){
      abscissa = {
        -0.9681602395076260898356,
        -0.8360311073266357942994,
        -0.6133714327005903973087,
        -0.3242534234038089290385,
        0,
        0.3242534234038089290385,
        0.6133714327005903973087,
        0.8360311073266357942994,
        0.9681602395076260898356
      };
      weights = {
        0.0812743883615744119719,
        0.1806481606948574040585,
        0.2606106964029354623187,
        0.312347077040002840069,
        0.330239355001259763165,
        0.312347077040002840069,
        0.260610696402935462319,
        0.1806481606948574040585,
        0.081274388361574411972
      };
    }
    if constexpr (N == 10){
      abscissa = {
        -0.973906528517171720078,
        -0.8650633666889845107321,
        -0.6794095682990244062343,
        -0.4333953941292471907993,
        -0.1488743389816312108848,
        0.1488743389816312108848,
        0.4333953941292471907993,
        0.6794095682990244062343,
        0.8650633666889845107321,
        0.973906528517171720078
      };
      weights = {
        0.0666713443086881375936,
        0.149451349150580593146,
        0.219086362515982043996,
        0.2692667193099963550912,
        0.2955242247147528701739,
        0.295524224714752870174,
        0.269266719309996355091,
        0.2190863625159820439955,
        0.1494513491505805931458,
        0.0666713443086881375936
      };
    }
    if constexpr (N == 11){
      abscissa = {
        -0.9782286581460569928039,
        -0.8870625997680952990752,
        -0.7301520055740493240934,
        -0.5190961292068118159257,
        -0.2695431559523449723315,
        0,
        0.269543155952344972332,
        0.5190961292068118159257,
        0.7301520055740493240934,
        0.887062599768095299075,
        0.9782286581460569928039
      };
      weights = {
        0.0556685671161736664828,
        0.1255803694649046246347,
        0.1862902109277342514261,
        0.2331937645919904799185,
        0.2628045445102466621807,
        0.2729250867779006307145,
        0.262804544510246662181,
        0.2331937645919904799185,
        0.1862902109277342514261,
        0.1255803694649046246347,
        0.055668567116173666483
      };
    }
    if constexpr (N == 12){
      abscissa = {
        -0.9815606342467192506906,
        -0.9041172563704748566785,
        -0.769902674194304687037,
        -0.5873179542866174472967,
        -0.3678314989981801937527,
        -0.1252334085114689154724,
        0.1252334085114689154724,
        0.3678314989981801937527,
        0.5873179542866174472967,
        0.7699026741943046870369,
        0.9041172563704748566785,
        0.9815606342467192506906
      };
      weights = {
        0.0471753363865118271946,
        0.1069393259953184309603,
        0.1600783285433462263347,
        0.2031674267230659217491,
        0.233492536538354808761,
        0.2491470458134027850006,
        0.2491470458134027850006,
        0.233492536538354808761,
        0.203167426723065921749,
        0.160078328543346226335,
        0.1069393259953184309603,
        0.0471753363865118271946
      };
    }
    if constexpr (N == 13){
      abscissa = {
        -0.9841830547185881494728,
        -0.9175983992229779652066,
        -0.8015780907333099127942,
        -0.642349339440340220644,
        -0.4484927510364468528779,
        -0.2304583159551347940655,
        0,
        0.2304583159551347940655,
        0.448492751036446852878,
        0.642349339440340220644,
        0.8015780907333099127942,
        0.9175983992229779652066,
        0.9841830547185881494728
      };
      weights = {
        0.04048400476531587952,
        0.0921214998377284479144,
        0.1388735102197872384636,
        0.1781459807619457382801,
        0.2078160475368885023125,
        0.2262831802628972384121,
        0.2325515532308739101946,
        0.2262831802628972384121,
        0.2078160475368885023125,
        0.17814598076194573828,
        0.138873510219787238464,
        0.0921214998377284479144,
        0.04048400476531587952
      };
    }
    if constexpr (N == 14){
      abscissa = {
        -0.9862838086968123388416,
        -0.9284348836635735173364,
        -0.82720131506976499319,
        -0.687292904811685470148,
        -0.515248636358154091965,
        -0.3191123689278897604357,
        -0.1080549487073436620662,
        0.1080549487073436620662,
        0.3191123689278897604357,
        0.515248636358154091965,
        0.687292904811685470148,
        0.82720131506976499319,
        0.9284348836635735173364,
        0.9862838086968123388416
      };
      weights = {
        0.0351194603317518630318,
        0.0801580871597602098056,
        0.1215185706879031846894,
        0.1572031671581935345696,
        0.185538397477937813742,
        0.2051984637212956039659,
        0.2152638534631577901959,
        0.215263853463157790196,
        0.205198463721295603966,
        0.185538397477937813742,
        0.15720316715819353457,
        0.121518570687903184689,
        0.0801580871597602098056,
        0.0351194603317518630318
      };
    }
    if constexpr (N == 15){
      abscissa = {
        -0.9879925180204854284896,
        -0.9372733924007059043078,
        -0.8482065834104272162007,
        -0.7244177313601700474162,
        -0.5709721726085388475372,
        -0.394151347077563369897,
        -0.2011940939974345223006,
        0,
        0.2011940939974345223006,
        0.394151347077563369897,
        0.5709721726085388475372,
        0.7244177313601700474162,
        0.8482065834104272162007,
        0.9372733924007059043078,
        0.9879925180204854284896
      };
      weights = {
        0.03075324199611726835463,
        0.0703660474881081247093,
        0.107159220467171935012,
        0.1395706779261543144478,
        0.1662692058169939335532,
        0.186161000015562211027,
        0.198431485327111576456,
        0.2025782419255612728806,
        0.198431485327111576456,
        0.1861610000155622110268,
        0.1662692058169939335532,
        0.1395706779261543144478,
        0.1071592204671719350119,
        0.070366047488108124709,
        0.0307532419961172683546
      };
    }
    if constexpr (N == 16){
      abscissa = {
        -0.9894009349916499325962,
        -0.944575023073232576078,
        -0.8656312023878317438805,
        -0.7554044083550030338951,
        -0.6178762444026437484467,
        -0.4580167776572273863424,
        -0.2816035507792589132305,
        -0.0950125098376374401853,
        0.0950125098376374401853,
        0.28160355077925891323,
        0.458016777657227386342,
        0.617876244402643748447,
        0.7554044083550030338951,
        0.8656312023878317438805,
        0.944575023073232576078,
        0.9894009349916499325962
      };
      weights = {
        0.02715245941175409485178,
        0.0622535239386478928628,
        0.0951585116824927848099,
        0.1246289712555338720525,
        0.1495959888165767320815,
        0.169156519395002538189,
        0.1826034150449235888668,
        0.1894506104550684962854,
        0.189450610455068496285,
        0.182603415044923588867,
        0.1691565193950025381893,
        0.149595988816576732082,
        0.124628971255533872052,
        0.0951585116824927848099,
        0.0622535239386478928628,
        0.0271524594117540948518
      };
    }
    if constexpr (N == 17){
      abscissa = {
         -0.9905754753144173356754,
        -0.9506755217687677612227,
        -0.880239153726985902123,
        -0.7815140038968014069252,
        -0.6576711592166907658503,
        -0.5126905370864769678863,
        -0.3512317634538763152972,
        -0.1784841814958478558507,
        0,
        0.1784841814958478558507,
        0.3512317634538763152972,
        0.5126905370864769678863,
        0.6576711592166907658503,
        0.7815140038968014069252,
        0.880239153726985902123,
        0.9506755217687677612227,
        0.9905754753144173356754     
      };
      weights = {
        0.0241483028685479319601,
        0.0554595293739872011294,
        0.0850361483171791808835,
        0.111883847193403971095,
        0.1351363684685254732863,
        0.1540457610768102880814,
        0.16800410215645004451,
        0.1765627053669926463253,
        0.1794464703562065254583,
        0.1765627053669926463253,
        0.16800410215645004451,
        0.1540457610768102880814,
        0.1351363684685254732863,
        0.111883847193403971095,
        0.0850361483171791808835,
        0.055459529373987201129,
        0.0241483028685479319601
      };
    }
    if constexpr (N == 18){
      abscissa = {
        -0.99156516842093094673,
        -0.9558239495713977551812,
        -0.8926024664975557392061,
        -0.8037049589725231156824,
        -0.6916870430603532078749,
        -0.559770831073947534608,
        -0.4117511614628426460359,
        -0.251886225691505509589,
        -0.0847750130417353012423,
        0.0847750130417353012423,
        0.251886225691505509589,
        0.4117511614628426460359,
        0.5597708310739475346079,
        0.6916870430603532078749,
        0.803704958972523115682,
        0.8926024664975557392061,
        0.9558239495713977551812,
        0.99156516842093094673
      };
      weights = {
        0.0216160135264833103133,
        0.0497145488949697964533,
        0.0764257302548890565291,
        0.100942044106287165563,
        0.122555206711478460185,
        0.1406429146706506512047,
        0.1546846751262652449254,
        0.1642764837458327229861,
        0.1691423829631435918407,
        0.169142382963143591841,
        0.164276483745832722986,
        0.1546846751262652449254,
        0.140642914670650651205,
        0.1225552067114784601845,
        0.100942044106287165563,
        0.0764257302548890565291,
        0.0497145488949697964533,
        0.0216160135264833103133
      };
    }
    if constexpr (N == 19){
      abscissa = {
        -0.992406843843584403189,
        -0.9602081521348300308528,
        -0.9031559036148179016427,
        -0.8227146565371428249789,
        -0.7209661773352293786171,
        -0.6005453046616810234696,
        -0.4645707413759609457173,
        -0.3165640999636298319901,
        -0.1603586456402253758681,
        0,
        0.1603586456402253758681,
        0.3165640999636298319901,
        0.4645707413759609457173,
        0.6005453046616810234696,
        0.7209661773352293786171,
        0.8227146565371428249789,
        0.9031559036148179016427,
        0.9602081521348300308528,
        0.992406843843584403189
      };
      weights = {
        0.0194617882297264770363,
        0.0448142267656996003328,
        0.0690445427376412265807,
        0.0914900216224499994645,
        0.111566645547333994716,
        0.1287539625393362276755,
        0.1426067021736066117758,
        0.152766042065859666779,
        0.15896884339395434765,
        0.1610544498487836959792,
        0.15896884339395434765,
        0.1527660420658596667789,
        0.1426067021736066117758,
        0.1287539625393362276755,
        0.111566645547333994716,
        0.091490021622449999464,
        0.0690445427376412265807,
        0.0448142267656996003328,
        0.0194617882297264770363
      };
    }
    if constexpr (N == 20){
      abscissa = {
        -0.9931285991850949247861,
        -0.9639719272779137912677,
        -0.9122344282513259058678,
        -0.8391169718222188233945,
        -0.7463319064601507926143,
        -0.6360536807265150254528,
        -0.5108670019508270980044,
        -0.3737060887154195606725,
        -0.2277858511416450780805,
        -0.07652652113349733375464,
        0.0765265211334973337546,
        0.2277858511416450780805,
        0.3737060887154195606725,
        0.5108670019508270980044,
        0.6360536807265150254528,
        0.7463319064601507926143,
        0.8391169718222188233945,
        0.9122344282513259058678,
        0.9639719272779137912677,
        0.9931285991850949247861
      };
      weights = {
        0.0176140071391521183119,
        0.04060142980038694133104,
        0.0626720483341090635695,
        0.0832767415767047487248,
        0.1019301198172404350368,
        0.1181945319615184173124,
        0.1316886384491766268985,
        0.1420961093183820513293,
        0.1491729864726037467878,
        0.1527533871307258506981,
        0.152753387130725850698,
        0.149172986472603746788,
        0.142096109318382051329,
        0.1316886384491766268985,
        0.1181945319615184173124,
        0.101930119817240435037,
        0.083276741576704748725,
        0.0626720483341090635695,
        0.040601429800386941331,
        0.0176140071391521183119
      };
    }
    
    // shift to domain x on [0,1]
    for(unsigned int i =0; i < N; ++i ) {
      abscissa[i] = 0.5*(abscissa[i] + 1);
      weights[i] *= 0.5;;
    }
  }
  

};

}

#endif
