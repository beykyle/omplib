// Copyright 2015 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// https://github.com/google/spherical-harmonics
//
// These functions have been modified from the source repository. Only the
// evalLegendrePolynomial() function has been retained from the original work.

#ifndef SPHERICAL_HARMONICS_HEADER
#define SPHERICAL_HARMONICS_HEADER

#include <cstdint>

namespace omplib {
/// @brief Evaluates the associated legendre polynomial functions
///
/// Evaluate the associated Legendre polynomial of degree, l, and order, m, at
/// coordinate, x. The inputs must satisfy:
/// 1. l >= 0
/// 2. -l <= m <= l
/// 3. -1 <= x <= 1
/// See http://en.wikipedia.org/wiki/Associated_Legendre_polynomials
///
/// Instead of computing Pml(x) directly, Pmm(x) is computed. Pmm can be
/// lifted to Pmm+1 recursively until Pml is found
const double
evalLegendrePolynomial(const uint8_t l, const int m, const double x);

/// @brief evaluates the (normalized) real spherical harmonics functions
const double realSphericalHarmonic(
    const uint8_t l, const int m, const double mu, const double omega);

} // namespace util

#endif
