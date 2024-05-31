// Copyright (c) 2023 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : William Wen, Sven Oesau, Martin Skrodzki

#ifndef CGAL_CONSTRAINT_BASED_SMOOTH_POINT_SET_H
#define CGAL_CONSTRAINT_BASED_SMOOTH_POINT_SET_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/Point_set_processing_3/internal/Neighbor_query.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Min_sphere_of_points_d_traits_3.h>
#include <CGAL/Min_sphere_of_spheres_d.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/iterator/zip_iterator.hpp>
#include <CGAL/Default_diagonalize_traits.h>

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL

namespace internal {

double dot(const std::array<double, 3>& a, const std::array<double, 3>& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double norm(const std::array<double, 3>& a) {
  return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
}

std::array<double, 6> dyadic(const std::array<double, 3>& a) {
  return { a[0] * a[0], a[0] * a[1], a[0] * a[2], a[1] * a[1], a[1] * a[2], a[2] * a[2] };
}

std::array<double, 3> operator+(const std::array<double, 3>& a, const std::array<double, 3>& b) {
  return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
}

std::array<double, 6> operator+(const std::array<double, 6>& a, const std::array<double, 6>& b) {
  return { a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3], a[4] + b[4], a[5] + b[5] };
}

std::array<double, 3> mul(double a, const std::array<double, 3>& v) {
  return { a * v[0], a * v[1], a * v[2] };
}

std::array<double, 3> sub(const std::array<double, 3>& a, const std::array<double, 3>& b) {
  return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
}

std::array<double, 3> mul(const std::array<double, 6>& m, const std::array<double, 3>& v) {
  return { m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
           m[1] * v[0] + m[3] * v[1] + m[4] * v[2],
           m[2] * v[0] + m[4] * v[1] + m[5] * v[2] };
}

std::array<double, 3> mul_inv(const std::array<double, 6>& m, const std::array<double, 3>& v) {
  std::array<double, 6> m_inv;
  m_inv[0] = m[5] * m[3] - m[4] * m[4];
  m_inv[1] = m[2] * m[4] - m[5] * m[1];
  m_inv[2] = m[1] * m[4] - m[2] * m[3];

  m_inv[3] = m[5] * m[0] - m[2] * m[2];
  m_inv[4] = m[1] * m[2] - m[0] * m[4];

  m_inv[5] = m[0] * m[3] - m[1] * m[1];

  // D = (m11 * a11) + (m12 * a12) + (m13 * a13)
  double det = m[0] * m_inv[0] + m[1] * m_inv[1] + m[2] * m_inv[2];
  if (det != 0) {
    m_inv[0] /= det;
    m_inv[1] /= det;
    m_inv[2] /= det;
    m_inv[3] /= det;
    m_inv[4] /= det;
    m_inv[5] /= det;

    return mul(m_inv, v);
  }
  else {
    return v;
  }
}

// Set eigenvalues of a matrix to 0 or 1 based on whether they are
// above/below a threshold
template <typename Kernel, typename Vector, typename Matrix>
void optimize_matrix_eigenvalues(
  std::vector<std::pair<Vector, Matrix>>& eigens,
  typename Kernel::FT eigenvalue_threshold
)
{
  typedef typename Kernel::FT FT;

  // get max eigenvalue
  FT max_eigenvalue = 0;
  FT avg_eigenvalue = 0;
  for (auto& eigen : eigens)
  {
    FT curr_max = eigen.first[2]; // largest eigen value
    avg_eigenvalue += curr_max;
    if (curr_max > max_eigenvalue)
    {
      max_eigenvalue = curr_max;
    }
  }
  avg_eigenvalue /= eigens.size();

  std::cout << "max eigenvalue: " << max_eigenvalue << "\n";
  std::cout << "avg eigenvalue: " << avg_eigenvalue << " / " << avg_eigenvalue / max_eigenvalue << "\n";

  // do binary optimization
  // scale so that eigenvalues are in range [0, 1]
  for (auto& eigen : eigens)
    for (int i = 0; i < 3; ++i)
      eigen.first[i] = (eigen.first[i] / max_eigenvalue) >= eigenvalue_threshold ? 1 : 0;
}

// Calculates the normal voting tensor and decomposes it into it's eigenvalues and eigenvectors
template <typename Kernel, typename PointRange,
  typename PointMap, typename VectorMap, typename DiagonalizeTraits>
std::pair<typename DiagonalizeTraits::Vector, typename DiagonalizeTraits::Matrix> calculate_nvt_eigenvalues(
  const typename PointRange::iterator::value_type& vt,
  PointMap point_map,
  VectorMap normal_map,
  const std::vector<typename PointRange::iterator>& neighbor_pwns,
  typename Kernel::FT normal_threshold,
  typename Kernel::FT eigenvalue_threshold,
  DiagonalizeTraits = DiagonalizeTraits())
{
  // basic geometric types
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Point_3 Point;
  typedef typename DiagonalizeTraits::Covariance_matrix Symmetric_matrix;
  typedef typename DiagonalizeTraits::Matrix Matrix;
  typedef typename DiagonalizeTraits::Vector Vector;

  FT weight = 0;
  Symmetric_matrix nvt = { 0, 0, 0, 0, 0, 0 };

  const Point& p = get(point_map, vt);
  const Vector_3& n = get(normal_map, vt);

  for (typename PointRange::iterator it : neighbor_pwns)
  {
    const Vector_3& nn = get(normal_map, *it);

    Vector vnn{ nn.x(), nn.y(), nn.z() };

    FT angle_diff = CGAL::approximate_angle(n, nn);
    if (angle_diff <= normal_threshold)
    {
      weight += 1;

      nvt[0] += vnn[0] * vnn[0];
      nvt[1] += vnn[0] * vnn[1];
      nvt[2] += vnn[0] * vnn[2];
      nvt[3] += vnn[1] * vnn[1];
      nvt[4] += vnn[1] * vnn[2];
      nvt[5] += vnn[2] * vnn[2];
    }
  }

  if (weight != 0)
  {
    nvt[0] /= weight;
    nvt[1] /= weight;
    nvt[2] /= weight;
    nvt[3] /= weight;
    nvt[4] /= weight;
    nvt[5] /= weight;
  }

  Vector eigenvalues;
  Matrix eigenvectors;
  DiagonalizeTraits::diagonalize_selfadjoint_covariance_matrix(nvt, eigenvalues, eigenvectors);

  return std::make_pair(eigenvalues, eigenvectors);
}

// Calculates the covariance matrix and decomposes it into it's eigenvalues and eigenvectors
template <typename Kernel, typename PointRange,
  typename PointMap, typename VectorMap, typename DiagonalizeTraits>
std::pair<typename DiagonalizeTraits::Vector, typename DiagonalizeTraits::Matrix> calculate_covm_eigenvalues(
  const typename PointRange::iterator::value_type& vt,
  PointMap point_map,
  VectorMap normal_map,
  const std::vector<typename PointRange::iterator>& neighbor_pwns,
  typename Kernel::FT normal_threshold,
  typename Kernel::FT eigenvalue_threshold,
  DiagonalizeTraits = DiagonalizeTraits())
{
  // basic geometric types
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Point_3 Point;
  typedef typename DiagonalizeTraits::Vector Vector;
  typedef typename DiagonalizeTraits::Covariance_matrix Symmetric_matrix;
  typedef typename DiagonalizeTraits::Matrix Matrix;

  FT weight = 0;
  Symmetric_matrix covm{ 0, 0, 0, 0, 0, 0 };
  Vector v_bar{ 0, 0, 0 };

  const Point& p = get(point_map, vt);
  const Vector_3& n = get(normal_map, vt);

  for (typename PointRange::iterator it : neighbor_pwns)
  {
    const Point& np = get(point_map, *it);
    const Vector_3& nn = get(normal_map, *it);

    Vector vnp{ np.x(), np.y(), np.z() };

    FT angle_diff = CGAL::approximate_angle(n, nn);
    if (angle_diff <= normal_threshold)
    {
      weight += 1;

      v_bar[0] += vnp[0];
      v_bar[1] += vnp[1];
      v_bar[2] += vnp[2];
    }
  }

  if (weight != 0)
  {
    v_bar[0] /= weight;
    v_bar[1] /= weight;
    v_bar[2] /= weight;
  }

  for (typename PointRange::iterator it : neighbor_pwns)
  {
    const Point& np = get(point_map, *it);
    const Vector_3& nn = get(normal_map, *it);

    Vector vnp{ np.x(), np.y(), np.z() };

    FT angle_diff = CGAL::approximate_angle(n, nn);
    if (angle_diff <= normal_threshold)
    {
      Vector temp_vec{ vnp[0] - v_bar[0], vnp[1] - v_bar[1], vnp[2] - v_bar[2] };
      covm[0] += temp_vec[0] * temp_vec[0];
      covm[1] += temp_vec[0] * temp_vec[1];
      covm[2] += temp_vec[0] * temp_vec[2];
      covm[3] += temp_vec[1] * temp_vec[1];
      covm[4] += temp_vec[1] * temp_vec[2];
      covm[5] += temp_vec[2] * temp_vec[2];
    }
  }

  if (weight != 0)
  {
    covm[0] /= weight;
    covm[1] /= weight;
    covm[2] /= weight;
    covm[3] /= weight;
    covm[4] /= weight;
    covm[5] /= weight;
  }

  Vector eigenvalues;
  Matrix eigenvectors;
  DiagonalizeTraits::diagonalize_selfadjoint_covariance_matrix(covm, eigenvalues, eigenvectors);

  return std::make_pair(eigenvalues, eigenvectors);
}

// Calculates new normal positions based off the normal voting tensor
template <typename Kernel, typename PointRange,
  typename PointMap, typename VectorMap,
  typename DiagonalizeTraits>
typename Kernel::Vector_3 nvt_normal_denoising(
  const typename PointRange::iterator::value_type& vt,
  PointMap point_map,
  VectorMap normal_map,
  const std::pair<typename DiagonalizeTraits::Vector, typename DiagonalizeTraits::Matrix> nvt_eigens,
  typename Kernel::FT damping_factor,
  DiagonalizeTraits = DiagonalizeTraits())
{
  // basic geometric types
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Point_3 Point;
  typedef typename DiagonalizeTraits::Covariance_matrix Symmetric_matrix;
  typedef typename DiagonalizeTraits::Vector Vector;

  Symmetric_matrix optimized_nvt{ 0, 0, 0, 0, 0, 0 };

  for (int i = 0; i < 3; ++i)
  {
    const int row_idx = 3 * i;
    optimized_nvt[0] += nvt_eigens.first[i] * nvt_eigens.second[row_idx] * nvt_eigens.second[row_idx];
    optimized_nvt[1] += nvt_eigens.first[i] * nvt_eigens.second[row_idx] * nvt_eigens.second[row_idx + 1];
    optimized_nvt[2] += nvt_eigens.first[i] * nvt_eigens.second[row_idx] * nvt_eigens.second[row_idx + 2];
    optimized_nvt[3] += nvt_eigens.first[i] * nvt_eigens.second[row_idx + 1] * nvt_eigens.second[row_idx + 1];
    optimized_nvt[4] += nvt_eigens.first[i] * nvt_eigens.second[row_idx + 1] * nvt_eigens.second[row_idx + 2];
    optimized_nvt[5] += nvt_eigens.first[i] * nvt_eigens.second[row_idx + 2] * nvt_eigens.second[row_idx + 2];
  }

  const Vector_3& n = get(normal_map, vt);
  Vector vn{ n.x(), n.y(), n.z() };

  Vector new_normal = mul(optimized_nvt, vn);
  new_normal[0] += damping_factor * vn[0];
  new_normal[1] += damping_factor * vn[1];
  new_normal[2] += damping_factor * vn[2];

  FT f = 1.0 / approximate_sqrt(dot(new_normal, new_normal));

  return Vector_3{ new_normal[0] * f, new_normal[1] * f, new_normal[2] * f };
}

enum point_type_t { corner = 3, edge = 1, flat = 2 };  // covm

// Detects type of point by using the covariance matrix eigenvalues
template <typename Kernel, typename Vector, typename Matrix>
point_type_t feature_detection(
  std::pair<Vector, Matrix> covm_eigens
) {
  int dominant_eigenvalue_count = 0;

  for (size_t i = 0; i < 3; ++i)
    if (covm_eigens.first[i] == 1)
      dominant_eigenvalue_count += 1;

  if (dominant_eigenvalue_count == 0)
    dominant_eigenvalue_count = 3; //corner

  return static_cast<point_type_t>(dominant_eigenvalue_count);
}

// Calculates new point position based on the point type and the covariance matrix
template <typename Kernel, typename PointRange,
  typename DiagonalizeTraits,
  typename PointMap, typename VectorMap>
typename Kernel::Point_3 calculate_new_point_position(
  const typename PointRange::iterator::value_type& vt,
  PointMap point_map,
  VectorMap normal_map,
  const std::vector<typename PointRange::iterator>& neighbor_pwns,
  const point_type_t point_type,
  std::pair<typename DiagonalizeTraits::Vector, typename DiagonalizeTraits::Matrix>& eigens,
  typename Kernel::FT update_threshold,
  typename Kernel::FT delta,
  DiagonalizeTraits = DiagonalizeTraits())
{
  // basic geometric types
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Point_3 Point;
  typedef typename DiagonalizeTraits::Vector Vector;
  typedef typename DiagonalizeTraits::Covariance_matrix Symmetric_matrix;

  // FT delta = 2.5;
  FT alpha = 0.1;

  Vector t{ 0, 0, 0 };

  const Point& p = get(point_map, vt);
  const Vector_3& n = get(normal_map, vt);

  Vector vp{ p.x(), p.y(), p.z() };
  Vector vn{ n.x(), n.y(), n.z() };

  if (point_type == point_type_t::corner)
  {
    Symmetric_matrix m_temp_symmetric{ 0, 0, 0, 0, 0 };
    Vector v_temp{ 0, 0, 0 };

    for (typename PointRange::iterator it : neighbor_pwns)
    {
      const Point& np = get(point_map, *it);
      const Vector_3& nn = get(normal_map, *it);

      Vector vnp{ np.x(), np.y(), np.z() };
      Vector vnn{ nn.x(), nn.y(), nn.z() };

      Symmetric_matrix dyadic_vnn = dyadic(vnn);
      m_temp_symmetric = m_temp_symmetric + dyadic_vnn;
      v_temp = v_temp + mul(dyadic_vnn, vnp);
    }

    t = mul_inv(m_temp_symmetric, v_temp);
  }
  else if (point_type == point_type_t::edge)
  {
    Vector dominant_eigenvec;

    for (int i = 0; i < 3; ++i)
      if (eigens.first[i] == 1)
      {
        dominant_eigenvec[0] = eigens.second[3 * i];
        dominant_eigenvec[1] = eigens.second[3 * i + 1];
        dominant_eigenvec[2] = eigens.second[3 * i + 2];
      }

    Symmetric_matrix m_temp{ 0, 0, 0, 0, 0, 0 };
    Vector v_temp{ 0, 0, 0 };

    for (typename PointRange::iterator it : neighbor_pwns)
    {
      const Point& np = get(point_map, *it);
      const Vector_3& nn = get(normal_map, *it);

      Vector vnp{ np.x(), np.y(), np.z() };
      Vector vnn{ nn.x(), nn.y(), nn.z() };

      Vector v_pi = sub(vnp, mul(dot(dominant_eigenvec, sub(vnp, vp)), dominant_eigenvec));
      Vector n_pi = sub(vnn, mul(dot(dominant_eigenvec, vnn), dominant_eigenvec));

      const std::array<double, 6> dyadic_n_pi = dyadic(n_pi);
      const std::array<double, 6> dyadic_eigenvec = dyadic(dominant_eigenvec);

      m_temp = m_temp + dyadic_n_pi + dyadic_eigenvec;

      v_temp = v_temp + mul(dyadic_n_pi, vnp) + mul(dyadic_eigenvec, vp);
    }

    t = mul_inv(m_temp, v_temp);
  }
  else if (point_type == point_type_t::flat)
  {
    FT cum_W = 0;
    Vector v_temp{ 0, 0, 0 };

    for (typename PointRange::iterator it : neighbor_pwns)
    {
      const Point& np = get(point_map, *it);
      const Vector_3& nn = get(normal_map, *it);

      Vector vnp{ np.x(), np.y(), np.z() };
      Vector vnn{ nn.x(), nn.y(), nn.z() };

      FT curr_W = std::exp(-(16 * norm(sub(vn, vnn)) / (delta * delta))) *
        std::exp(-(4 * norm(sub(vp, vnp)) / (delta * delta)));

      cum_W += curr_W;

      v_temp = v_temp + mul(curr_W * dot(vnn, sub(vnp, vp)), vn);
    }

    v_temp = mul((alpha / cum_W), v_temp);

    t = vp + v_temp;
  }

  double no = norm(sub(vp, t));

  if (no < update_threshold * update_threshold)
    return Point(t[0], t[1], t[2]);

  return Point(vp[0], vp[1], vp[2]);
}

} /* namespace internal */

/// \endcond



// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

template <typename ConcurrencyTag,
  typename PointRange,
  typename NamedParameters = parameters::Default_named_parameters
>
void
constraint_based_smooth_point_set(
  PointRange& points,
  const NamedParameters& np = parameters::default_values()
)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // basic geometric types
  typedef typename PointRange::iterator iterator;
  typedef std::vector<iterator> iterators;
  typedef typename iterator::value_type value_type;
  typedef Point_set_processing_3_np_helper<PointRange, NamedParameters> NP_helper;
  typedef typename NP_helper::Point_map PointMap;
  typedef typename NP_helper::Normal_map NormalMap;
  typedef typename NP_helper::Geom_traits Kernel;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;

  typedef typename GetDiagonalizeTraits<NamedParameters, double, 3>::type DiagonalizeTraits;
  typedef typename DiagonalizeTraits::Matrix Matrix;
  typedef typename DiagonalizeTraits::Vector Vector;

  CGAL_assertion_msg(NP_helper::has_normal_map(points, np), "Error: no normal map");

  typedef typename Kernel::FT FT;

  FT neighbor_radius = 0;
  FT normal_threshold = 0.9 * (180 / CGAL_PI);
  FT damping_factor = 3;
  FT eigenvalue_threshold_nvt = 0.7;
  FT eigenvalue_threshold_covm = .2;
  FT update_threshold = 0;

  bool do_point_smoothing = true;

  CGAL_precondition(points.begin() != points.end());

  // types for neighbors search structure
  typedef Point_set_processing_3::internal::Neighbor_query<Kernel, PointRange&, PointMap> Neighbor_query;

  PointMap point_map = NP_helper::get_point_map(points, np);
  NormalMap normal_map = NP_helper::get_normal_map(points, np);

  std::size_t nb_points = points.size();

  // initiate a KD-tree search for points
  Neighbor_query neighbor_query(points, point_map);

  // automatic parameter calculation
  if (neighbor_radius == 0)
  {
    FT avg_dist = CGAL::compute_average_spacing<ConcurrencyTag>(
      points, 16,
      CGAL::parameters::point_map(point_map));

    neighbor_radius = 2 * avg_dist;

    std::cout << "neighbor radius: " << neighbor_radius << "\n";
  }

  if (update_threshold == 0)
  {
    update_threshold = 2 * neighbor_radius;

    std::cout << "update threshold: " << update_threshold << "\n";
  }

  // compute all neighbors
  std::vector<iterators> pwns_neighbors;
  pwns_neighbors.resize(nb_points);

  typedef boost::zip_iterator<boost::tuple<iterator, typename std::vector<iterators>::iterator> > Zip_iterator;

  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range(boost::make_zip_iterator(boost::make_tuple(points.begin(), pwns_neighbors.begin())),
      boost::make_zip_iterator(boost::make_tuple(points.end(), pwns_neighbors.end()))),
      [&](const typename Zip_iterator::reference& t)
      {
        neighbor_query.get_iterators(get(point_map, get<0>(t)), 0, neighbor_radius,
          std::back_inserter(get<1>(t)));
        return true;
      });

  // compute the optimized normal voting tensors
  typedef boost::zip_iterator
    <boost::tuple<iterator,
    typename std::vector<iterators>::iterator,
    typename std::vector<std::pair<Vector, Matrix>>::iterator> > Zip_iterator_2;

  std::vector<std::pair<Vector, Matrix>> nvt_optimized_eigens(nb_points);

  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range(boost::make_zip_iterator(boost::make_tuple
    (points.begin(), pwns_neighbors.begin(), nvt_optimized_eigens.begin())),
      boost::make_zip_iterator(boost::make_tuple
      (points.end(), pwns_neighbors.end(), nvt_optimized_eigens.end()))),
      [&](const typename Zip_iterator_2::reference& t)
      {
        get<2>(t) = internal::calculate_nvt_eigenvalues<Kernel, PointRange>
          (get<0>(t),
            point_map, normal_map,
            get<1>(t),
            normal_threshold,
            eigenvalue_threshold_nvt,
            DiagonalizeTraits());
        return true;
      });

  // optimize eigenvalues
  internal::optimize_matrix_eigenvalues<Kernel>(nvt_optimized_eigens, eigenvalue_threshold_nvt);

  // compute the new normal for each point
  std::vector<typename Kernel::Vector_3> new_normals(nb_points);

  typedef boost::zip_iterator
    <boost::tuple<iterator,
    typename std::vector<std::pair<Vector, Matrix>>::iterator,
    typename std::vector<typename Kernel::Vector_3>::iterator> > Zip_iterator_3;

  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range(boost::make_zip_iterator(boost::make_tuple
    (points.begin(), nvt_optimized_eigens.begin(), new_normals.begin())),
      boost::make_zip_iterator(boost::make_tuple
      (points.end(), nvt_optimized_eigens.end(), new_normals.end()))),
      [&](const typename Zip_iterator_3::reference& t)
      {
        get<2>(t) = internal::nvt_normal_denoising<Kernel, PointRange>
          (get<0>(t),
            point_map, normal_map,
            get<1>(t),
            damping_factor,
            DiagonalizeTraits());
        return true;
      });

  // update the normals
  typedef boost::zip_iterator
    <boost::tuple<iterator,
    typename std::vector<typename Kernel::Vector_3>::iterator> > Zip_iterator_4;

  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range(boost::make_zip_iterator(boost::make_tuple(points.begin(), new_normals.begin())),
      boost::make_zip_iterator(boost::make_tuple(points.end(), new_normals.end()))),
      [&](const typename Zip_iterator_4::reference& t)
      {
        put(normal_map, get<0>(t), get<1>(t));
        return true;
      });

  // early return if point smoothing is turned off
  if (!do_point_smoothing)
    return;

  // compute the covariance matrix for each point
  std::vector<std::pair<Vector, Matrix>> covm_optimized_eigens(nb_points);

  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range(boost::make_zip_iterator(boost::make_tuple
    (points.begin(), pwns_neighbors.begin(), covm_optimized_eigens.begin())),
      boost::make_zip_iterator(boost::make_tuple
      (points.end(), pwns_neighbors.end(), covm_optimized_eigens.end()))),
      [&](const typename Zip_iterator_2::reference& t)
      {
        get<2>(t) = internal::calculate_covm_eigenvalues<Kernel, PointRange>
          (get<0>(t),
            point_map, normal_map,
            get<1>(t),
            normal_threshold,
            eigenvalue_threshold_covm,
            DiagonalizeTraits());
        return true;
      });

  // optimize eigenvalues
  internal::optimize_matrix_eigenvalues<Kernel>(covm_optimized_eigens, eigenvalue_threshold_covm);

  // classify each point based on the cov matrix
  std::vector<internal::point_type_t> point_classifications;
  point_classifications.reserve(nb_points);

  // for (auto it = nvt_optimized_eigens.begin() ; it < nvt_optimized_eigens.end() ; ++it)
  for (auto it = covm_optimized_eigens.begin(); it < covm_optimized_eigens.end(); ++it)
    point_classifications.push_back(internal::feature_detection<Kernel>(*it));

  // color the points FOR DEBUG ONLY
  for (size_t i = 0; i < point_classifications.size(); ++i)
    switch (point_classifications[i])
    {
    case internal::point_type_t::corner:
      get<2>(points[i]) = CGAL::make_array(static_cast<unsigned char>(255), static_cast<unsigned char>(0), static_cast<unsigned char>(0));
      break;
    case internal::point_type_t::edge:
      get<2>(points[i]) = CGAL::make_array(static_cast<unsigned char>(0), static_cast<unsigned char>(255), static_cast<unsigned char>(0));
      break;
    case internal::point_type_t::flat:
      get<2>(points[i]) = CGAL::make_array(static_cast<unsigned char>(0), static_cast<unsigned char>(0), static_cast<unsigned char>(255));
      break;
    }

  // calculate diameter of point set
  typedef CGAL::Min_sphere_of_points_d_traits_3<Kernel, FT> Traits;
  typedef CGAL::Min_sphere_of_spheres_d<Traits> Min_sphere;
  Min_sphere ms;
  for (auto it = points.begin(); it < points.end(); ++it)
    ms.insert(get<0>(*it));

  FT delta = ms.radius();

  std::cout << "radius of set: " << delta << "\n";

  // compute the new point
  std::vector<Point_3> new_points(nb_points);

  typedef boost::zip_iterator
    <boost::tuple<iterator,
    typename std::vector<iterators>::iterator,
    typename std::vector<internal::point_type_t>::iterator,
    typename std::vector<std::pair<Vector, Matrix>>::iterator,
    typename std::vector<Point_3>::iterator> > Zip_iterator_5;

  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range(boost::make_zip_iterator(boost::make_tuple
    (points.begin(), pwns_neighbors.begin(), point_classifications.begin(), covm_optimized_eigens.begin(), new_points.begin())),
      boost::make_zip_iterator(boost::make_tuple
      (points.end(), pwns_neighbors.end(), point_classifications.end(), covm_optimized_eigens.end(), new_points.end()))),
      [&](const typename Zip_iterator_5::reference& t)
      {
        get<4>(t) = internal::calculate_new_point_position<Kernel, PointRange>
          (get<0>(t),
            point_map, normal_map,
            get<1>(t),
            get<2>(t),
            get<3>(t),
            update_threshold,
            delta,
            DiagonalizeTraits());
        return true;
      });

  // update the new point
  typedef boost::zip_iterator
    <boost::tuple<iterator,
    typename std::vector<typename Kernel::Point_3>::iterator> > Zip_iterator_6;
  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range(boost::make_zip_iterator(boost::make_tuple(points.begin(), new_points.begin())),
      boost::make_zip_iterator(boost::make_tuple(points.end(), new_points.end()))),
      [&](const typename Zip_iterator_6::reference& t)
      {
        put(point_map, get<0>(t), get<1>(t));
        return true;
      });

}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_CONSTRAINT_BASED_SMOOTH_POINT_SET_H
