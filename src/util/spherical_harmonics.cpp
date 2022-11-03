#include "util/spherical_harmonics.hpp"
#include "util/factorial.hpp"

#include <cassert>
#include <cmath>

using namespace omplib;

const double
omplib::evalLegendrePolynomial(const uint8_t l, const int m, const double x) {
  assert((std::abs(m) <= l));
  // Compute Pmm(x) = (-1)^m(2m - 1)!!(1 - x^2)^(m/2), where !! is the double
  // factorial.
  double pmm = 1.0, minus_m_factor = 1.0;

  if (m < 0) {
    minus_m_factor = (1 - ((m & 1) << 1)) *
                     static_cast<double>(factorial(l + m)) /
                     static_cast<double>(factorial(l - m));
  }

  // P00 is defined as 1.0, do don't evaluate Pmm unless we know m > 0
  if (abs(m) > 0) {
    pmm = static_cast<double>((1 - ((m & 1) << 1))) *
          doubleFactorial(2 * abs(m) - 1) * pow(1 - x * x, abs(m) / 2.0);
  }

  if (l == abs(m)) {
    // Pml is the same as Pmm so there's no lifting to higher bands needed
    return pmm * minus_m_factor;
  }

  // Compute Pmm+1(x) = x(2m + 1)Pmm(x)
  double pmm1 = x * (2 * m + 1) * pmm;
  if (l == abs(m) + 1) {
    // Pml is the same as Pmm+1 so we are done as well
    return pmm1 * minus_m_factor;
  }

  // Use the last two computed bands to lift up to the next band until l is
  // reached, using the recurrence relationship:
  // Pml(x) = (x(2l - 1)Pml-1 - (l + m - 1)Pml-2) / (l - m)
  for (int n = abs(m) + 2; n <= l; n++) {
    double pmn =
        (x * (2 * n - 1) * pmm1 - (n + abs(m) - 1) * pmm) / (n - abs(m));
    pmm = pmm1;
    pmm1 = pmn;
  }

  // Pmm1 at the end of the above loop is equal to Pml
  return pmm1 * minus_m_factor;
}

const double omplib::realSphericalHarmonic(
    const uint8_t l, const int m, const double mu, const double omega) {
  assert((std::abs(m) <= l));

  return std::sqrt(
             (2 * l + 1) * static_cast<double>(factorial(l - m)) /
             static_cast<double>(factorial(l + m))) *
         evalLegendrePolynomial(l, m, mu) * std::cos(m * omega);
}
