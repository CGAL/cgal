// Copyright (c) 2013-06  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : 

#ifndef CGAL_CONSTRAINT_BASED_SMOOTH_POINT_SET_H
#define CGAL_CONSTRAINT_BASED_SMOOTH_POINT_SET_H

#include <CGAL/license/Point_set_processing_3.h>

#include <cmath>
#include <CGAL/disable_warnings.h>

#include <CGAL/Point_set_processing_3/internal/Neighbor_query.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Min_sphere_of_points_d_traits_3.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
// #include <CGAL/Point_set_processing_3/internal/Callback_wrapper.h>
#include <CGAL/for_each.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/iterator/zip_iterator.hpp>

#include <CGAL/Eigen_vector.h>
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Kernel/global_functions.h>

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL

namespace internal {

template <typename Kernel>
void optimize_matrix_eigenvalues(
  std::vector<std::pair<Eigen::Vector3d, Eigen::Matrix3d>>& eigens,
  typename Kernel::FT eigenvalue_threshold
)
{
  typedef typename Kernel::FT FT;

  // get max eigenvalue
  FT max_eigenvalue = 0;
  FT avg_eigenvalue = 0;
  for(auto& eigen : eigens)
  {
    FT curr_max = eigen.first.maxCoeff();
    avg_eigenvalue += curr_max;
    if(curr_max > max_eigenvalue)
    {
      max_eigenvalue = curr_max;
    }
  }
  avg_eigenvalue /= eigens.size();

  std::cout << "max eigenvalue: " << max_eigenvalue << "\n";
  std::cout << "avg eigenvalue: " << avg_eigenvalue << " / " << avg_eigenvalue/max_eigenvalue << "\n";

  // do binary optimization
  // scale so that eigenvalues are in range [0, 1]
  for(auto& eigen : eigens)
  {
    for (int i=0; i<3; ++i)
    {
      // std::cout << eigens.first[i] / max_eigenvalue << "\n";
      eigen.first[i] = eigen.first[i] / max_eigenvalue >= eigenvalue_threshold ? 1 : 0;
      // eigen.first[i] = eigen.first[i] >= eigenvalue_threshold ? 1 : 0;
    }
  }
  // std::cout << "\n";

}

// Finds the optimized eigenvalues (and eigenvectors) of the normal voting
// tensors at each point
template <typename Kernel, typename PointRange,
        typename PointMap, typename VectorMap>
std::pair<Eigen::Vector3d, Eigen::Matrix3d> calculate_optimized_nvt_eigenvalues(
  const typename PointRange::iterator::value_type& vt,
  PointMap point_map,
  VectorMap normal_map,
  const std::vector<typename PointRange::iterator>& neighbor_pwns,
  typename Kernel::FT normal_threshold,
  typename Kernel::FT eigenvalue_threshold)
{
  // basic geometric types
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Point_3 Point;

  FT weight = 0;
  Eigen::Matrix3d nvt = Eigen::Matrix3d::Zero();

  const Point& p = get(point_map, vt);
  const Vector& n = get(normal_map, vt);

  for (typename PointRange::iterator it : neighbor_pwns)
  {
    const Vector& nn = get(normal_map, *it);

    Eigen::Vector3d vnn(nn.x(), nn.y(), nn.z());

    FT angle_diff = CGAL::approximate_angle(n, nn);
    if (angle_diff <= normal_threshold)
    {
      weight += 1;

      nvt += vnn * vnn.transpose();
    }
  }

  if (weight != 0)  // idk if exact equality is okay
  {
    nvt /= weight;
  }

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(nvt);
  std::pair<Eigen::Vector3d, Eigen::Matrix3d> eigens(solver.eigenvalues(), solver.eigenvectors());

  return eigens;
}

template <typename Kernel, typename PointRange,
        typename PointMap, typename VectorMap>
std::pair<Eigen::Vector3d, Eigen::Matrix3d> calculate_optimized_covm_eigenvalues(
  const typename PointRange::iterator::value_type& vt,
  PointMap point_map,
  VectorMap normal_map,
  const std::vector<typename PointRange::iterator>& neighbor_pwns,
  typename Kernel::FT normal_threshold,
  typename Kernel::FT eigenvalue_threshold)
{
  // basic geometric types
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Point_3 Point;

  FT weight = 0;
  Eigen::Matrix3d covm = Eigen::Matrix3d::Zero();
  Eigen::Vector3d v_bar = Eigen::Vector3d::Zero();

  const Point& p = get(point_map, vt);
  const Vector& n = get(normal_map, vt);

  for (typename PointRange::iterator it : neighbor_pwns)
  {
    const Point& np = get(point_map, *it);
    const Vector& nn = get(normal_map, *it);

    Eigen::Vector3d vnp(np.x(), np.y(), np.z());

    FT angle_diff = CGAL::approximate_angle(n, nn);
    if (angle_diff <= normal_threshold)
    {
      weight += 1;

      v_bar += vnp;
    }
  }

  if (weight != 0)  // idk if exact equality is okay
  {
    v_bar /= weight;
  }

  for (typename PointRange::iterator it : neighbor_pwns)
  {
    const Point& np = get(point_map, *it);
    const Vector& nn = get(normal_map, *it);

    Eigen::Vector3d vnp(np.x(), np.y(), np.z());

    FT angle_diff = CGAL::approximate_angle(n, nn);
    if (angle_diff <= normal_threshold)
    {
      auto temp_vec = vnp  - v_bar;
      covm += temp_vec * temp_vec.transpose();
    }
  }

  if (weight != 0)  // idk if exact equality is okay
  {
    covm /= weight;
  }

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covm);
  std::pair<Eigen::Vector3d, Eigen::Matrix3d> eigens(solver.eigenvalues(), solver.eigenvectors());

  return eigens;
}

template <typename Kernel, typename PointRange,
          typename PointMap, typename VectorMap>
typename Kernel::Vector_3 nvt_normal_denoising(
  const typename PointRange::iterator::value_type& vt,
  PointMap point_map,
  VectorMap normal_map,
  const std::pair<Eigen::Vector3d, Eigen::Matrix3d> nvt_eigens,
  typename Kernel::FT damping_factor)
{
  // basic geometric types
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Point_3 Point;

  Eigen::Matrix3d optimized_nvt = Eigen::Matrix3d::Zero();

  for (int i=0; i<3; ++i)
  {
    optimized_nvt += nvt_eigens.first[i] * nvt_eigens.second.col(i) * nvt_eigens.second.col(i).transpose();
  }

  const Vector& n = get(normal_map, vt);
  Eigen::Vector3d vn(n.x(), n.y(), n.z());

  Eigen::Vector3d new_normal = damping_factor * vn + optimized_nvt * vn;
  new_normal.normalize();

  return Vector{new_normal[0], new_normal[1], new_normal[2]};
}

enum point_type_t {corner = 3, edge = 1, flat = 2};  // covm
// enum point_type_t {corner = 3, edge = 2, flat = 1};  // nvt

template <typename Kernel>
point_type_t feature_detection(
  std::pair<Eigen::Vector3d, Eigen::Matrix3d> covm_eigens
) {
  int dominant_eigenvalue_count = 0;

  for (size_t i=0 ; i<3 ; ++i)
  {
    if (covm_eigens.first[i] == 1)
    {
      dominant_eigenvalue_count += 1;
    }
  }

  if (dominant_eigenvalue_count == 0)
  {
    dominant_eigenvalue_count = 3; //corner  
  }

  return static_cast<point_type_t>(dominant_eigenvalue_count);
}

template <typename Kernel, typename PointRange,
          typename PointMap, typename VectorMap>
typename Kernel::Point_3 calculate_new_point(
  const typename PointRange::iterator::value_type& vt,
  PointMap point_map,
  VectorMap normal_map,
  const std::vector<typename PointRange::iterator>& neighbor_pwns,
  const point_type_t point_type,
  std::pair<Eigen::Vector3d, Eigen::Matrix3d>& eigens,
  typename Kernel::FT update_threshold,
  typename Kernel::FT delta)
{ 
  // basic geometric types
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Point_3 Point;

  // FT delta = 2.5;
  FT alpha = 0.1;

  Eigen::Vector3d t = Eigen::Vector3d::Zero();

  const Point& p = get(point_map, vt);
  const Vector& n = get(normal_map, vt);

  Eigen::Vector3d vp(p.x(), p.y(), p.z());
  Eigen::Vector3d vn(n.x(), n.y(), n.z());

  if (point_type == point_type_t::corner)
  {
    Eigen::Matrix3d m_temp = Eigen::Matrix3d::Zero();
    Eigen::Vector3d v_temp = Eigen::Vector3d::Zero();

    for (typename PointRange::iterator it : neighbor_pwns)
    {
      const Point& np = get(point_map, *it);
      const Vector& nn = get(normal_map, *it);

      Eigen::Vector3d vnp(np.x(), np.y(), np.z());
      Eigen::Vector3d vnn(nn.x(), nn.y(), nn.z());

      m_temp += vnn * vnn.transpose();
      v_temp += (vnn * vnn.transpose()) * vnp;
    }

    t = m_temp.inverse() * v_temp;
  }
  else if (point_type == point_type_t::edge)
  {
    Eigen::Vector3d dominant_eigenvec;

    for (int i=0; i<3; ++i)
    {
      if(eigens.first[i] == 1)
      {
        dominant_eigenvec = eigens.second.col(i);
      }
    }

    Eigen::Matrix3d m_temp = Eigen::Matrix3d::Zero();
    Eigen::Vector3d v_temp = Eigen::Vector3d::Zero();

    for (typename PointRange::iterator it : neighbor_pwns)
    {
      const Point& np = get(point_map, *it);
      const Vector& nn = get(normal_map, *it);

      Eigen::Vector3d vnp(np.x(), np.y(), np.z());
      Eigen::Vector3d vnn(nn.x(), nn.y(), nn.z());

      Eigen::Vector3d v_pi = vnp - (dominant_eigenvec.dot(vnp - vp)) * dominant_eigenvec;
      Eigen::Vector3d n_pi = vnn - (dominant_eigenvec.dot(vnn)) * dominant_eigenvec;

      m_temp += (n_pi * n_pi.transpose()) + (dominant_eigenvec * dominant_eigenvec.transpose());

      v_temp += (n_pi * n_pi.transpose()) * vnp + (dominant_eigenvec * dominant_eigenvec.transpose()) * vp;
    }

    t = m_temp.inverse() * v_temp;
  }
  else if (point_type == point_type_t::flat)
  {
    FT cum_W = 0;
    Eigen::Vector3d v_temp = Eigen::Vector3d::Zero();

    for (typename PointRange::iterator it : neighbor_pwns)
    {
      const Point& np = get(point_map, *it);
      const Vector& nn = get(normal_map, *it);

      Eigen::Vector3d vnp(np.x(), np.y(), np.z());
      Eigen::Vector3d vnn(nn.x(), nn.y(), nn.z());

      FT curr_W = std::exp(- (16*(vn - vnn).squaredNorm()) / (delta * delta)) *
        std::exp(- (4*(vp - vnp).squaredNorm()) / (delta * delta));

      cum_W += curr_W;

      v_temp += curr_W * vnn.dot(vnp - vp) * vn;
    }

    v_temp = (alpha / cum_W) * v_temp;

    t = vp + v_temp;
  }

  if ((vp - t).squaredNorm() < update_threshold * update_threshold)
  {
    return Point(t[0], t[1], t[2]);
  }

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

  CGAL_assertion_msg(NP_helper::has_normal_map(points, np), "Error: no normal map");

  typedef typename Kernel::FT FT;

  FT neighbor_radius = 0;
  FT normal_threshold = 0.9 * (180/M_PI);
  FT damping_factor = 3;
  FT eigenvalue_threshold_nvt = .3;
  FT eigenvalue_threshold_covm = .3;
  FT update_threshold = 0;

  bool do_point_smoothing = true;

  CGAL_precondition(points.begin() != points.end());

  // types for neighbors search structure
  typedef Point_set_processing_3::internal::Neighbor_query<Kernel, PointRange&, PointMap> Neighbor_query;

  PointMap point_map = NP_helper::get_point_map(points, np);
  NormalMap normal_map = NP_helper::get_normal_map(points, np);

  std::size_t nb_points = points.size();

  // initiate a KD-tree search for points
  Neighbor_query neighbor_query (points, point_map);
  
  // automatic parameter calculation
  if (neighbor_radius == 0)
  {
    FT avg_dist = CGAL::compute_average_spacing<ConcurrencyTag>(
                         points, 6,
                         CGAL::parameters::point_map(point_map));

    neighbor_radius = 2 * avg_dist;
    update_threshold = 2 * neighbor_radius;

    std::cout << "neighbor radius: " << neighbor_radius << "\n";
    std::cout << "update threshold: " << update_threshold << "\n";
  }

  // compute all neighbors
  std::vector<iterators> pwns_neighbors;
  pwns_neighbors.resize(nb_points);

  typedef boost::zip_iterator<boost::tuple<iterator, typename std::vector<iterators>::iterator> > Zip_iterator;

  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range (boost::make_zip_iterator (boost::make_tuple (points.begin(), pwns_neighbors.begin())),
                      boost::make_zip_iterator (boost::make_tuple (points.end(), pwns_neighbors.end()))),
    [&](const typename Zip_iterator::reference& t)
    {
      neighbor_query.get_iterators (get(point_map, get<0>(t)), 0, neighbor_radius,
                                          std::back_inserter (get<1>(t)));
      return true;
    });

  // compute the optimized normal voting tensors
  typedef boost::zip_iterator
    <boost::tuple<iterator,
                  typename std::vector<iterators>::iterator,
                  typename std::vector<std::pair<Eigen::Vector3d, Eigen::Matrix3d>>::iterator> > Zip_iterator_2;

  std::vector<std::pair<Eigen::Vector3d, Eigen::Matrix3d>> nvt_optimized_eigens(nb_points);

  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range (boost::make_zip_iterator (boost::make_tuple
                                                (points.begin(), pwns_neighbors.begin(), nvt_optimized_eigens.begin())),
                      boost::make_zip_iterator (boost::make_tuple
                                                (points.end(), pwns_neighbors.end(), nvt_optimized_eigens.end()))),
    [&](const typename Zip_iterator_2::reference& t)
    {
      get<2>(t) = internal::calculate_optimized_nvt_eigenvalues<Kernel, PointRange>
          (get<0>(t),
          point_map, normal_map,
          get<1>(t),
          normal_threshold,
          eigenvalue_threshold_nvt);
      return true;
    });

  // optimize eigenvalues
  internal::optimize_matrix_eigenvalues<Kernel>(nvt_optimized_eigens, eigenvalue_threshold_nvt);

  // compute the new normal for each point
  std::vector<typename Kernel::Vector_3> new_normals(nb_points);

  typedef boost::zip_iterator
    <boost::tuple<iterator,
                  typename std::vector<std::pair<Eigen::Vector3d, Eigen::Matrix3d>>::iterator,
                  typename std::vector<typename Kernel::Vector_3>::iterator> > Zip_iterator_3;
  
  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range (boost::make_zip_iterator (boost::make_tuple
                                                (points.begin(), nvt_optimized_eigens.begin(), new_normals.begin())),
                      boost::make_zip_iterator (boost::make_tuple
                                                (points.end(), nvt_optimized_eigens.end(), new_normals.end()))),
    [&](const typename Zip_iterator_3::reference& t)
    {
      get<2>(t) = internal::nvt_normal_denoising<Kernel, PointRange>
          (get<0>(t),
          point_map, normal_map,
          get<1>(t),
          damping_factor);
      return true;
    });

  // update the normals
  typedef boost::zip_iterator
    <boost::tuple<iterator,
                  typename std::vector<typename Kernel::Vector_3>::iterator> > Zip_iterator_4;

  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range (boost::make_zip_iterator (boost::make_tuple (points.begin(), new_normals.begin())),
                      boost::make_zip_iterator (boost::make_tuple (points.end(), new_normals.end()))),
    [&](const typename Zip_iterator_4::reference& t)
    {
      put (normal_map, get<0>(t), get<1>(t));
      return true;
    });

  // early return if point smoothing is turned off
  if(!do_point_smoothing)
  {
    return;
  }

  // compute the covariance matrix for each point
  std::vector<std::pair<Eigen::Vector3d, Eigen::Matrix3d>> covm_optimized_eigens(nb_points);

  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range (boost::make_zip_iterator (boost::make_tuple
                                                (points.begin(), pwns_neighbors.begin(), covm_optimized_eigens.begin())),
                      boost::make_zip_iterator (boost::make_tuple
                                                (points.end(), pwns_neighbors.end(), covm_optimized_eigens.end()))),
    [&](const typename Zip_iterator_2::reference& t)
    {
      get<2>(t) = internal::calculate_optimized_covm_eigenvalues<Kernel, PointRange>
          (get<0>(t),
           point_map, normal_map,
           get<1>(t),
           normal_threshold,
           eigenvalue_threshold_covm);
      return true;
    });

  // optimize eigenvalues
  internal::optimize_matrix_eigenvalues<Kernel>(covm_optimized_eigens, eigenvalue_threshold_covm);

  // classify each point based on the cov matrix
  std::vector<internal::point_type_t> point_classifications;
  point_classifications.reserve(nb_points);

  // for (auto it = nvt_optimized_eigens.begin() ; it < nvt_optimized_eigens.end() ; ++it)
  for (auto it = covm_optimized_eigens.begin() ; it < covm_optimized_eigens.end() ; ++it)
  {
    point_classifications.push_back(internal::feature_detection<Kernel>(*it));
  }

  // colour the points FOR DEBUG ONLY
  for (size_t i=0; i<point_classifications.size(); ++i)
  {
    switch(point_classifications[i])
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
  }

  // calculate diameter of point set
  typedef CGAL::Min_sphere_of_points_d_traits_3<Kernel,FT> Traits;
  typedef CGAL::Min_sphere_of_spheres_d<Traits> Min_sphere;
  Min_sphere ms;
  for (auto it = points.begin() ; it < points.end() ; ++it)
  {
    ms.insert(get<0>(*it));
  }
  FT delta = ms.radius();

  std::cout << "radius of set: " << delta << "\n";

  // compute the new point
  std::vector<Point_3> new_points(nb_points);

  typedef boost::zip_iterator
    <boost::tuple<iterator,
                  typename std::vector<iterators>::iterator,
                  typename std::vector<internal::point_type_t>::iterator,
                  typename std::vector<std::pair<Eigen::Vector3d, Eigen::Matrix3d>>::iterator,
                  typename std::vector<Point_3>::iterator> > Zip_iterator_5;

  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range (boost::make_zip_iterator (boost::make_tuple
                                                (points.begin(), pwns_neighbors.begin(), point_classifications.begin(), covm_optimized_eigens.begin(), new_points.begin())),
                      boost::make_zip_iterator (boost::make_tuple
                                                (points.end(), pwns_neighbors.end(), point_classifications.end(), covm_optimized_eigens.end(), new_points.end()))),
    [&](const typename Zip_iterator_5::reference& t)
    {
      get<4>(t) = internal::calculate_new_point<Kernel, PointRange>
          (get<0>(t),
           point_map, normal_map,
           get<1>(t),
           get<2>(t),
           get<3>(t),
           update_threshold,
           delta);
      return true;
    });

  // update the new point
  typedef boost::zip_iterator
    <boost::tuple<iterator,
                  typename std::vector<typename Kernel::Point_3>::iterator> > Zip_iterator_6;
  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range (boost::make_zip_iterator (boost::make_tuple (points.begin(), new_points.begin())),
                       boost::make_zip_iterator (boost::make_tuple (points.end(), new_points.end()))),
     [&](const typename Zip_iterator_6::reference& t)
     {
       put (point_map, get<0>(t), get<1>(t));
       return true;
     });

}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_CONSTRAINT_BASED_SMOOTH_POINT_SET_H
