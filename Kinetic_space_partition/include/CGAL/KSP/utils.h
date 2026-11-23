// Copyright (c) 2020 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot, Dmitry Anisimov

#ifndef CGAL_KSP_UTILS_H
#define CGAL_KSP_UTILS_H

#include <CGAL/license/Kinetic_space_partition.h>

// STL includes.
#include <set>
#include <algorithm>
#include <cmath>
#include <array>
#include <string>
#include <sstream>
#include <functional>
#include <fstream>
#include <vector>
#include <deque>
#include <queue>
#include <map>
#include <iostream>

// CGAL includes.
#include <CGAL/Iterator_range.h>
#include <CGAL/number_utils.h>
#include <CGAL/assertions.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Boost includes.
#include <boost/iterator/function_output_iterator.hpp>

namespace CGAL {
namespace KSP {
namespace internal {

#ifdef DOXYGEN_RUNNING
#else

// Diagonalization of 3x3 covariance matrices for the plane fitting of the input polygons if Eigen is not available.
// The cyclic Jacobi method iterates until the off-diagonal entries vanish in floating-point precision. Thus, the result
// agrees with the one of Eigen up to rounding errors. The deprecated CGAL::Diagonalize_traits, the fallback of
// CGAL::Default_diagonalize_traits if Eigen is not available, stops at a fixed tolerance and its results deviate
// considerably from the ones of Eigen.
// The model of `DiagonalizeTraits` is restricted to the dimension 3.
template<typename FT>
struct Jacobi_diagonalize_traits_3 {
  using Vector = std::array<FT, 3>;
  using Matrix = std::array<FT, 9>;
  using Covariance_matrix = std::array<FT, 6>;

  // Eigenvalues are sorted by increasing order. The eigenvector of eigenvalues[k] is stored in eigenvectors[3 * k + i], i = 0, 1, 2.
  static bool diagonalize_selfadjoint_covariance_matrix(const Covariance_matrix& cov, Vector& eigenvalues, Matrix& eigenvectors) {
    // The covariance matrix stores the upper triangle row-wise: (0,0), (0,1), (0,2), (1,1), (1,2), (2,2).
    double a[3][3] = {
      { CGAL::to_double(cov[0]), CGAL::to_double(cov[1]), CGAL::to_double(cov[2]) },
      { CGAL::to_double(cov[1]), CGAL::to_double(cov[3]), CGAL::to_double(cov[4]) },
      { CGAL::to_double(cov[2]), CGAL::to_double(cov[4]), CGAL::to_double(cov[5]) } };
    double v[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

    for (int sweep = 0; sweep < 64; sweep++) {
      if (a[0][1] == 0 && a[0][2] == 0 && a[1][2] == 0)
        break;

      for (int p = 0; p < 2; p++)
        for (int q = p + 1; q < 3; q++) {
          const double apq = a[p][q];
          if (apq == 0)
            continue;

          // The entry is negligible compared to both diagonal entries.
          const double g = 100 * std::abs(apq);
          if (sweep > 3 && std::abs(a[p][p]) + g == std::abs(a[p][p]) && std::abs(a[q][q]) + g == std::abs(a[q][q])) {
            a[p][q] = a[q][p] = 0;
            continue;
          }

          // Rotation in the (p, q) plane annihilating a[p][q].
          const double h = a[q][q] - a[p][p];
          double t;
          if (std::abs(h) + g == std::abs(h))
            t = apq / h;
          else {
            const double theta = 0.5 * h / apq;
            t = 1 / (std::abs(theta) + std::sqrt(1 + theta * theta));
            if (theta < 0)
              t = -t;
          }
          const double c = 1 / std::sqrt(1 + t * t);
          const double s = t * c;

          for (int k = 0; k < 3; k++) {
            const double akp = a[k][p], akq = a[k][q];
            a[k][p] = c * akp - s * akq;
            a[k][q] = s * akp + c * akq;
          }
          for (int k = 0; k < 3; k++) {
            const double apk = a[p][k], aqk = a[q][k];
            a[p][k] = c * apk - s * aqk;
            a[q][k] = s * apk + c * aqk;
          }
          a[p][q] = a[q][p] = 0;

          for (int k = 0; k < 3; k++) {
            const double vkp = v[k][p], vkq = v[k][q];
            v[k][p] = c * vkp - s * vkq;
            v[k][q] = s * vkp + c * vkq;
          }
        }
    }

    std::array<int, 3> order = { 0, 1, 2 };
    std::sort(order.begin(), order.end(), [&](int i, int j) { return a[i][i] < a[j][j] || (a[i][i] == a[j][j] && i < j); });
    for (int k = 0; k < 3; k++) {
      eigenvalues[k] = FT(a[order[k]][order[k]]);
      for (int i = 0; i < 3; i++)
        eigenvectors[3 * k + i] = FT(v[i][order[k]]);
    }
    return true;
  }

  static bool diagonalize_selfadjoint_covariance_matrix(const Covariance_matrix& cov, Vector& eigenvalues) {
    Matrix eigenvectors;
    return diagonalize_selfadjoint_covariance_matrix(cov, eigenvalues, eigenvectors);
  }

  static bool extract_largest_eigenvector_of_covariance_matrix(const Covariance_matrix& cov, Vector& normal) {
    Vector eigenvalues;
    Matrix eigenvectors;
    diagonalize_selfadjoint_covariance_matrix(cov, eigenvalues, eigenvectors);
    for (int i = 0; i < 3; i++)
      normal[i] = eigenvectors[6 + i];
    return true;
  }
};

// Collects the output of one line and writes it, including the line break, into the stream
// with a single output operation when going out of scope, e.g.:
//   Log_line("[P3] ") << "number of events: " << n;
// A single output operation is not interleaved with output of other threads by the standard
// stream implementations. Thus, lines written in parallel stay intact without additional locking.
class Log_line {
public:
  explicit Log_line(const std::string& prefix = std::string(), std::ostream& os = std::cout) : m_os(os) {
    m_buffer.precision(os.precision());
    m_buffer << prefix;
  }

  Log_line(const Log_line&) = delete;
  Log_line& operator=(const Log_line&) = delete;

  ~Log_line() {
    m_buffer << '\n';
    m_os << m_buffer.str();
    m_os.flush();
  }

  template<typename T>
  Log_line& operator<<(const T& t) {
    m_buffer << t;
    return *this;
  }

  // Manipulators, e.g., std::setprecision.
  Log_line& operator<<(std::ostream& (*manipulator)(std::ostream&)) {
    m_buffer << manipulator;
    return *this;
  }

private:
  std::ostringstream m_buffer;
  std::ostream& m_os;
};

// Convert point to string.
template<typename Point_d>
const std::string to_string(const Point_d& p) {
  std::ostringstream oss;
  oss.precision(20);
  oss << p;
  return oss.str();
}

// Distance between two points.
template<typename Point_d>
decltype(auto) distance(const Point_d& p, const Point_d& q) {
  using Traits = typename Kernel_traits<Point_d>::Kernel;
  using FT = typename Traits::FT;
  const FT sq_dist = CGAL::squared_distance(p, q);
  return static_cast<FT>(CGAL::approximate_sqrt(sq_dist));
}

// Normalize vector.
template<typename Vector_d>
inline const Vector_d normalize(const Vector_d& v) {
  using Traits = typename Kernel_traits<Vector_d>::Kernel;
  using FT = typename Traits::FT;
  const FT dot_product = CGAL::abs(v * v);
  //CGAL_assertion(dot_product != FT(0));
  return v / static_cast<FT>(CGAL::approximate_sqrt(dot_product));
}

// Intersections. Used only in the 2D version.
// For the 3D version, see conversions.h!
template<typename Type1, typename Type2, typename ResultType>
inline bool intersection(
  const Type1& t1, const Type2& t2, ResultType& result) {

  const auto inter = intersection(t1, t2);
  if (!inter) return false;

  if (CGAL::assign(result, inter))
    return true;

  return false;
}

template<typename ResultType, typename Type1, typename Type2>
inline const ResultType intersection(const Type1& t1, const Type2& t2) {

  ResultType out;
  CGAL_assertion_code(const bool is_intersection_found =) intersection(t1, t2, out);
  CGAL_assertion(is_intersection_found);
  return out;
}

// Angles.

// Converts radians to degrees.
template<typename FT>
const FT degrees_2(const FT angle_rad) {
  return angle_rad * FT(180) / static_cast<FT>(CGAL_PI);
}

// Computes an angle in degrees between two directions.
template<typename Direction_2>
const typename Kernel_traits<Direction_2>::Kernel::FT
compute_angle_2(const Direction_2& dir1, const Direction_2& dir2) {

  using Traits = typename Kernel_traits<Direction_2>::Kernel;
  using FT = typename Traits::FT;

  const auto v1 = dir2.to_vector();
  const auto v2 = -dir1.to_vector();

  const FT det = CGAL::determinant(v1, v2);
  const FT dot = CGAL::scalar_product(v1, v2);
  const FT angle_rad = static_cast<FT>(
    std::atan2(CGAL::to_double(det), CGAL::to_double(dot)));
  const FT angle_deg = degrees_2(angle_rad);
  return angle_deg;
}

// Converts an angle in degrees from the range [-180, 180]
// into the mod 90 angle.
template<typename FT>
const FT convert_angle_2(const FT angle_2) {

  FT angle = angle_2;
  if (angle > FT(90)) angle = FT(180) - angle;
  else if (angle < -FT(90)) angle = FT(180) + angle;
  return angle;
}

// Computes a positive angle in degrees that
// is always in the range [0, 90].
template<typename Direction_2>
const typename Kernel_traits<Direction_2>::Kernel::FT
angle_2(const Direction_2& dir1, const Direction_2& dir2) {
  const auto angle_2 = compute_angle_2(dir1, dir2);
  return CGAL::abs(convert_angle_2(angle_2));
}

#endif

} // namespace internal
} // namespace KSP
} // namespace CGAL

#endif // CGAL_KSP_UTILS_H
