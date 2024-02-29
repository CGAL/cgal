// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Rao Fu

#ifndef CGAL_ESTIMATE_LFS_H
#define CGAL_ESTIMATE_LFS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/IO/trace.h>
#include <CGAL/centroid.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Point_set_processing_3/internal/Neighbor_query.h>
#include <CGAL/Point_set_processing_3/internal/Callback_wrapper.h>
#include <CGAL/for_each.h>
#include <CGAL/property_map.h>
#include <CGAL/Index_property_map.h>
#include <CGAL/assertions.h>
#include <CGAL/Memory_sizer.h>
#include <functional>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include <vector>
#include <map>
#include <unordered_map>
#include <queue>
#include <mutex>

namespace CGAL {

// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL
namespace internal {

template <typename Geom_traits, typename PointRange, typename PointMap>
typename Geom_traits::Sphere_3
calculate_loose_bsphere(const PointRange& points, const PointMap point_map,
  const typename Geom_traits::FT loose_ratio=0.1)
{
  typedef typename Geom_traits::FT  FT;
  typedef typename Geom_traits::Point_3  Point;
  typedef typename Geom_traits::Sphere_3  Sphere;

  CGAL_precondition(points.begin() != points.end());

  Point center = CGAL::centroid(CGAL::make_transform_iterator_from_property_map
                                       (points.begin(), point_map),
                  CGAL::make_transform_iterator_from_property_map
                                       (points.end(), point_map));

  FT max_distance = (std::numeric_limits<double>::min)();
  for (const auto& pt : points)
  {
    const Point& point = get(point_map, pt);
    FT distance = CGAL::squared_distance(center, point);
    if (distance > max_distance) {
            max_distance = distance;
    }
  }

  return Sphere(center, max_distance*(1.0 + loose_ratio));
}

template <typename NeighborQuery>
typename NeighborQuery::FT
classical_point_dist_func(const typename NeighborQuery::Point_3& query, ///< point
                          const NeighborQuery& neighbor_query)
{
  // basic geometric types
  typedef typename NeighborQuery::Kernel Kernel;
  typedef typename NeighborQuery::FT  FT;
  typedef typename Kernel::Point_3  Point;
  std::vector<Point> points;
  // find the nearest neighbor
  neighbor_query.get_points(query, 1, FT(0), std::back_inserter(points));

  const FT dsq = CGAL::squared_distance(query, points[0]);

  return std::sqrt(dsq);
}

template <typename NeighborQuery>
typename NeighborQuery::FT
av_sqrt_sum_sq_distance_to_knn(const typename NeighborQuery::Point_3& query, ///< point
                                 const NeighborQuery& neighbor_query,
                                 const int knn)
{
  // basic geometric types
  typedef typename NeighborQuery::Kernel Kernel;
  typedef typename NeighborQuery::FT  FT;
  typedef typename Kernel::Point_3  Point;
  std::vector<Point> points;
  // find the nearest neighbor
  neighbor_query.get_points(query, knn, FT(0), std::back_inserter(points));

  // get distances
  FT sum_sq_distances = FT(0.0);
  for (const auto& point : points)
  {
    sum_sq_distances += CGAL::squared_distance(query, point);
  }

  std::size_t sz = points.size();

  if (sz == 0)
        return 0.0;
    else
        return std::sqrt(sum_sq_distances / sz);
}

template <typename NeighborQuery, typename PointRange, typename PointMap>
typename NeighborQuery::FT
classical_point_dist_func_epsilon_band(const PointRange &points,
                                      const NeighborQuery& neighbor_query,
                                      const PointMap point_map,
                                      const int band_knn=12)
{
  // basic geometric types
  typedef typename NeighborQuery::Kernel Kernel;
  typedef typename NeighborQuery::FT  FT;
  typedef typename Kernel::Point_3  Point;

  CGAL_precondition(points.begin() != points.end());
  CGAL_precondition(band_knn > 0);

  FT l_min = std::numeric_limits<FT>::infinity();
  for (const auto& pt : points)
  {
    const Point &query = get(point_map, pt);

    FT l = av_sqrt_sum_sq_distance_to_knn(query, neighbor_query, band_knn);

    if (l == FT(0)) continue;

    if (l < l_min)
    {
      l_min = l;
    }
  }

  return l_min;
}

template <typename Point, typename Vector, typename Sphere>
Point calculate_ray_sphere_intersection(const Point &p, const Vector &n, const Sphere &sphere)
{

    Vector d = p - sphere.center();

    auto d_len = std::sqrt(d.squared_length());

    auto cos_theta = -1.0 * (d * n) / d_len / std::sqrt(n.squared_length());

    auto l = d_len * cos_theta +
           std::sqrt(sphere.squared_radius() - d.squared_length() * (1.0 - cos_theta * cos_theta));

    Point t = p + n * l;

    return t;
}

template <typename Point, typename Vector>
Point local2world(const Point &local, const Point &o,
                  const Vector &u_basis, const Vector &v_basis, const Vector &w_basis)
{
    Vector v = local.x() * u_basis + local.y() * v_basis + local.z() * w_basis;
    Point world = o + v;

    return world;
}


template <typename Geom_traits>
void random_dual_cone_search_segs(const typename Geom_traits::Point_3 &p,
                    const typename Geom_traits::Vector_3 &n,
                      const typename Geom_traits::Sphere_3 &bsphere,
                      std::vector<typename Geom_traits::Segment_3> &segs,
                      const typename Geom_traits::FT half_apex_angle, std::size_t N)
{
  typedef typename Geom_traits::FT  FT;
  typedef typename Geom_traits::Point_3  Point;
  typedef typename Geom_traits::Vector_3  Vector;
  typedef typename Geom_traits::Segment_3  Segment;

  FT half_apex_radian = half_apex_angle * CGAL_PI / 180.0;

  // add the normal directions
  const Point s = calculate_ray_sphere_intersection(p, n, bsphere);
  const Point t = calculate_ray_sphere_intersection(p, -1.0 * n, bsphere);
  segs.push_back(Segment(s, t));

  // add rays within the cone

  // define local coordinate system
  // normal should be normalized
  FT eps = 1e-5;
  assert(std::abs(n.squared_length() - 1.0) <= eps);
  const Vector w_basis = n;

  Vector v_vec(0.0, 0.0, 0.0);
  if (std::abs(n.x()) >= std::abs(n.y()))
      v_vec = Vector(n.z(), 0.0, -1.0 * n.x());
  else
      v_vec = Vector(0.0, n.z(), -1.0 * n.y());
  // normalize v_vec
  FT v_vec_norm = v_vec.squared_length();
  v_vec_norm = v_vec_norm > eps ? v_vec_norm : eps;
  const Vector v_basis = v_vec / std::sqrt(v_vec_norm);

  const Vector u_vec = CGAL::cross_product(v_basis, w_basis);
  // normalize u_vec
  FT u_vec_norm = u_vec.squared_length();
  u_vec_norm = u_vec_norm > eps ? u_vec_norm : eps;
  const Vector u_basis = u_vec / std::sqrt(u_vec_norm);

  FT w = std::cos(half_apex_radian);
  assert(w <= 1.0);

  FT R = std::sqrt(1.0 - w * w);

  for (int i = 0; i < N; i++)
  {
    FT rd1 = rand() / (double)RAND_MAX;
    FT rd2 = rand() / (double)RAND_MAX;

    FT random_r = R * std::sqrt(rd1);
    FT random_phi = 2.0 * CGAL_PI * rd2;

    FT u = random_r * std::cos(random_phi);
    FT v = random_r * std::sin(random_phi);

    // local to world
    const Point world = local2world(Point(u, v, w), p, u_basis, v_basis, w_basis);
    const Vector ray_direction = world - p;
    const Point s = calculate_ray_sphere_intersection(p, ray_direction, bsphere);
    const Point t = calculate_ray_sphere_intersection(p, -1.0 * ray_direction, bsphere);
    segs.push_back(Segment(s, t));
  }
}

template <typename NeighborQuery>
bool recursive_dichotomic_search_base(const typename NeighborQuery::Kernel::Segment_3 &e,
                                      const typename NeighborQuery::Kernel::Point_2 &s,
                                      const typename NeighborQuery::Kernel::Point_2 &t,
                                      const NeighborQuery& neighbor_query,
                                      const typename NeighborQuery::FT epsilon_band,
                                      std::vector<typename NeighborQuery::Kernel::Point_3>& intersections,
                                      typename NeighborQuery::FT lipschitz,
                                      typename NeighborQuery::FT eps)
{
  // basic geometric types
  typedef typename NeighborQuery::Kernel Kernel;
  typedef typename NeighborQuery::FT  FT;
  typedef typename Kernel::Point_3  Point;
  typedef typename Kernel::Point_2  Point_2;
  typedef typename Kernel::Vector_3  Vector;

  assert(t.x() + eps >= s.x() - eps);
  CGAL_precondition(lipschitz >= 0.0);

  Point ps = e.source();
  Point pt = e.target();
  Vector v = pt - ps;
  FT len = std::sqrt(v.squared_length());

  // should satisfy lipschitz continuity -- keep it?
  FT abs_dy = std::abs(t.y() - s.y());
  FT abs_dx = std::abs(t.x() - s.x());
  // it is better to abs_dy / abs_dx, when abs_dx is small
  if (abs_dy / abs_dx > lipschitz + eps)
  {
    return false;
  }

  FT lip_lbx = (s.x() + t.x()) * 0.5 - (t.y() - s.y()) * 0.5 / lipschitz;
  FT lip_lby = s.y() - lipschitz * (lip_lbx - s.x());

  assert(lip_lbx + eps >= s.x() - eps);
  assert(lip_lbx - eps <= t.x() + eps);

  FT lip_ubx = ((t.y() - s.y()) / lipschitz + t.x() + s.x()) * 0.5;
  FT lip_uby = s.y() + lipschitz * (lip_ubx - s.x());

  assert(lip_ubx + eps >= s.x() - eps);
  assert(lip_ubx - eps <= t.x() + eps);

  // stop
  if (lip_lby >= epsilon_band || lip_uby <= epsilon_band)
  {
    return true;
  }

  const bool s_sign = (s.y() > epsilon_band);
  const bool t_sign = (t.y() > epsilon_band);
  const FT s_slope = (s_sign) ? -lipschitz : lipschitz;
  const FT t_slope = (t_sign) ? lipschitz : -lipschitz;

  FT ms_x = s.x() + (epsilon_band - s.y()) / s_slope;
  const Point pms = ps + ms_x / len * v;
  FT ms_y = classical_point_dist_func(pms, neighbor_query);
  const Point_2 ms(ms_x, ms_y);
  assert(ms_x + eps >= s.x() - eps);

  FT mt_x = t.x() + (epsilon_band - t.y()) / t_slope;
  const Point pmt = ps + mt_x / len * v;
  FT mt_y = classical_point_dist_func(pmt, neighbor_query);
  const Point_2 mt(mt_x, mt_y);
  assert(mt_x - eps <= t.x() + eps);
  assert(mt_x + eps >= ms_x - eps);

  const FT mm_x = (ms_x + mt_x) / 2;
  const Point pmm = ps + mm_x / len * v;
  FT mm_y = classical_point_dist_func(pmm, neighbor_query);
  Point_2 mm(mm_x, mm_y);

  if (std::abs(mm_y - epsilon_band) <= eps && (mt_x - ms_x) <= epsilon_band / lipschitz)
  {
    if (s_sign == t_sign)
    {
            intersections.push_back(pmm);
            intersections.push_back(pmm);
    }
    else
    {
      intersections.push_back(pmm);
    }

      return true;
  }

  bool fl = recursive_dichotomic_search_base(e, ms, mm, neighbor_query, epsilon_band,
            intersections, lipschitz, eps);
  bool fr = recursive_dichotomic_search_base(e, mm, mt, neighbor_query, epsilon_band,
            intersections, lipschitz, eps);

  return (fl && fr);
}

template<class NeighborQuery>
bool recursive_dichotomic_search(const typename NeighborQuery::Kernel::Segment_3 &e,
                                const typename NeighborQuery::Kernel::Point_2 &s,
                                const typename NeighborQuery::Kernel::Point_2 &t,
                                const NeighborQuery &neighbor_query,
                                const typename NeighborQuery::FT epsilon_band,
                                std::vector<typename NeighborQuery::Kernel::Point_3>& intersections,
                                typename NeighborQuery::FT lipschitz,
                                typename NeighborQuery::FT eps=1e-5)
{
  // basic geometric types
  typedef typename NeighborQuery::Kernel Kernel;
  typedef typename Kernel::Point_3  Point;

  CGAL_precondition(intersections.empty());
  bool flag = recursive_dichotomic_search_base(e, s, t, neighbor_query, epsilon_band,
                    intersections, lipschitz, eps);

  if (intersections.size() % 2 != 0)
  {
    const bool s_sign = (s.y() > epsilon_band);
    const bool t_sign = (t.y() > epsilon_band);

    if ((!s_sign) && t_sign)
    {
      std::vector<Point> intersections_;
      intersections_.reserve(intersections.size());
      intersections_.push_back(e.source());
      for (std::size_t i = 1; i < intersections.size(); i = i + 2)
      {
        const Point p0 = intersections[i];
        const Point p1 = intersections[i + 1];
        const Point p = p1 + (p0 - p1) / 2;
        intersections_.push_back(p);
      }
      std::swap(intersections, intersections_);
    }

    if ((!t_sign) && s_sign)
    {
      std::vector<Point> intersections_;
      intersections_.reserve(intersections.size());
      for (std::size_t i = 0; i < intersections.size() - 1; i = i + 2)
      {
        const Point p0 = intersections[i];
        const Point p1 = intersections[i + 1];
        const Point p = p1 + (p0 - p1) / 2;
        intersections_.push_back(p);
      }
      intersections_.push_back(e.target());
      std::swap(intersections, intersections_);
    }

    if (s_sign && t_sign)
    {
      flag = false;
    }

    // should not happen
    if ((!s_sign) && (!t_sign))
    {
      flag = false;
    }
  }
  else
  {
    std::vector<Point> intersections_;
    for (std::size_t i = 0; i < intersections.size(); i = i + 2)
    {
      const Point p0 = intersections[i];
      const Point p1 = intersections[i + 1];
      const Point p = p1 + (p0 - p1) / 2;
      intersections_.push_back(p);
    }
      intersections = intersections_;
  }

  return flag;
}

template <typename NeighborQuery>
typename NeighborQuery::FT
estimate_shape_diameter(const typename NeighborQuery::Point_3& query, ///< point
           const typename NeighborQuery::Kernel::Vector_3& normal,
           const NeighborQuery& neighbor_query, ///< KD-tree
           const typename NeighborQuery::Kernel::Sphere_3& bsphere,
           const typename NeighborQuery::FT epsilon_band,
           const typename NeighborQuery::FT apex_angle,
           std::size_t N_rays,
           unsigned int reject_knn=6)
{
  // basic geometric types
  typedef typename NeighborQuery::Kernel Kernel;
  typedef typename NeighborQuery::FT  FT;
  typedef typename Kernel::Point_3  Point;
  typedef typename Kernel::Point_2  Point_2;
  typedef typename Kernel::Segment_3  Segment;

  std::vector<Point> points;
  neighbor_query.get_points (query, reject_knn, FT(0), std::back_inserter(points));
  FT project_radius = reject_knn == 0 ? 0.0 : std::sqrt(CGAL::squared_distance(query, points[reject_knn-1]));

  // dual cone search segments
  FT half_apex_angle = apex_angle / 2.0;
  std::vector<Segment> segs;
  random_dual_cone_search_segs<Kernel>(query, normal, bsphere, segs,
            half_apex_angle, N_rays);

  std::vector<FT> antipodal_point_dsqs;
  for (const auto seg : segs)
  {
    const Point ps = seg.source();
    const Point pt = seg.target();
    FT l = std::sqrt(seg.squared_length());

    FT ls = classical_point_dist_func(ps, neighbor_query);
    FT lt = classical_point_dist_func(pt, neighbor_query);
    FT lipschitz = 1.1; // slightly enlarge the lipschitz
    std::vector<Point> intersections;
    bool flag = recursive_dichotomic_search(seg, Point_2(0.0, ls), Point_2(l, lt), neighbor_query,
                                    epsilon_band, intersections, lipschitz);

    if (flag && intersections.size() >= 2)
    {
      std::vector<FT> dsqs;
      for (std::size_t i = 0; i < intersections.size(); i++)
      {
        const auto intersection = intersections[i];

        FT dsq = (query - intersection).squared_length();
        dsqs.push_back(dsq);
      }
      std::sort(dsqs.begin(), dsqs.end(), std::less<FT>());
      // dist should be larger than epsilon_band
      auto it = upper_bound(dsqs.begin(), dsqs.end(), project_radius*project_radius);

      if (it != dsqs.end())
      {
        FT dsq_min = *it;
        antipodal_point_dsqs.push_back(dsq_min);
      }
    }
  }

  FT dsq_upper_bound = 4.0 * bsphere.squared_radius();

  if(antipodal_point_dsqs.empty())
        antipodal_point_dsqs.push_back(dsq_upper_bound);

  // robust distance function here
  FT sum_squared_distances = std::accumulate(antipodal_point_dsqs.begin(), antipodal_point_dsqs.end(), 0.0);
  std::size_t sz = antipodal_point_dsqs.size();
  if (sz == 0)
      return 0.0;
  else
      return std::sqrt(sum_squared_distances / sz);
}

template <typename Monge_form>
typename Monge_form::FT
min_curvature_radius(Monge_form &monge_form)
{
  // basic geometric types
  typedef typename Monge_form::FT FT;

  FT max_principal_curvature = monge_form.principal_curvatures(0);
  FT min_principal_curvature = monge_form.principal_curvatures(1);

  FT r1 = 1.0 / std::abs(max_principal_curvature);
  FT r2 = 1.0 / std::abs(min_principal_curvature);
  FT r = (std::min)(r1, r2);

  return r;
}

template <typename NeighborQuery>
typename NeighborQuery::FT
estimate_local_feature_size(const typename NeighborQuery::Point_3& query, ///< point
                    typename NeighborQuery::Kernel::Vector_3& normal,
                    const NeighborQuery& neighbor_query, ///< KD-tree
                    const typename NeighborQuery::Kernel::Sphere_3& bsphere, ///< bounding sphere
                    unsigned int jet_k, ///< number of neighbors
                    typename NeighborQuery::FT neighbor_radius,
                    unsigned int degree_fitting,
                    unsigned int degree_monge,
                    typename NeighborQuery::FT epsilon_band,
                    typename NeighborQuery::FT apex_angle,
                    std::size_t N_rays
                    )
{
    // basic geometric types
    typedef typename NeighborQuery::Kernel Kernel;
    typedef typename Kernel::Point_3  Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::FT FT;

    // types for jet fitting
    typedef Monge_via_jet_fitting< Kernel > Monge_jet_fitting;
    typedef typename Monge_jet_fitting::Monge_form Monge_form;

    std::vector<Point> points;
    // query using as fallback minimum requires nb points for jet fitting (d+1)*(d+2)/2
    neighbor_query.get_points (query, jet_k, neighbor_radius, std::back_inserter(points),
                             (degree_fitting + 1) * (degree_fitting + 2) / 2);

    // estimate jet
    Monge_jet_fitting monge_fit;
    Monge_form monge_form = monge_fit(points.begin(), points.end(), degree_fitting, degree_monge);

    if (normal * normal == 0.0)
        normal = monge_form.normal_direction(); // already normalized in monge_form
    else
        monge_form.comply_wrt_given_normal(normal); // orient jet with valid given normal

    FT abs_curvature_r = min_curvature_radius(monge_form);

    // estimate shape diameter
    FT shape_diameter = estimate_shape_diameter(query, normal, neighbor_query, bsphere, epsilon_band, apex_angle, N_rays);
    FT half_shape_diameter = 0.5 * shape_diameter;

    FT lfs = (std::min)(abs_curvature_r, half_shape_diameter);

    return lfs;
}

template <typename ConcurrencyTag,
          typename PointRange,
          typename NamedParameters = parameters::Default_named_parameters>
void
construct_knn_graph(PointRange& points,
                    const int knn,
                    std::vector<std::vector<typename PointRange::iterator>>& knn_graph,
                    const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  // basic geometric types
  typedef Point_set_processing_3_np_helper<PointRange, NamedParameters> NP_helper;
  typedef typename NP_helper::Point_map PointMap;
  typedef typename NP_helper::Geom_traits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point;

  // types for K nearest neighbors search structure
  typedef typename PointRange::iterator iterator;
  typedef Point_set_processing_3::internal::Neighbor_query<Kernel, PointRange&, PointMap> Neighbor_query;
  // Property map typename PointRange::iterator -> index
  typedef Index_property_map<iterator> IndexMap;

  PointMap point_map = NP_helper::get_point_map(points, np);

  const std::function<bool(double)>& callback = choose_parameter(get_parameter(np, internal_np::callback),
                                                               std::function<bool(double)>());

  CGAL_precondition(points.begin() != points.end());

  std::size_t memory = CGAL::Memory_sizer().virtual_size();
  CGAL_TRACE_STREAM << (memory >> 20) << " Mb allocated\n";
  CGAL_TRACE_STREAM << "  Creates KD-tree\n";

  Neighbor_query neighbor_query(points, point_map);

  std::size_t nb_points = points.size();
  knn_graph.clear();
  knn_graph.resize(nb_points);

  Point_set_processing_3::internal::Callback_wrapper<ConcurrencyTag>
  callback_wrapper (callback, nb_points);
  typedef boost::zip_iterator<boost::tuple<iterator, typename std::vector<std::vector<iterator>>::iterator> > Zip_iterator;
  CGAL::for_each<ConcurrencyTag>
  (CGAL::make_range (boost::make_zip_iterator (boost::make_tuple (points.begin(), knn_graph.begin())),
                        boost::make_zip_iterator (boost::make_tuple (points.end(), knn_graph.end()))),
   [&](const typename Zip_iterator::reference& t)
   {
    if (callback_wrapper.interrupted())
      return false;

    const Point& point = get(point_map, boost::get<0>(t));
    neighbor_query.get_iterators(point, knn, FT(0), std::back_inserter(boost::get<1>(t)));

    ++ callback_wrapper.advancement();
    return true;
   });
}

} /* namespace internal */
/// \endcond

/**
   \ingroup PkgPointSetProcessing3Algorithms

   Estimate the local feature size (LFS) for the input 3D point cloud. The function only works for 3D point cloud.
   If the input 3D point cloud has no normals, the function will also estimate the normals using jet fitting.


   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range
   \param jet_k number of neighbors for jet fitting
   \param N_rays number of rays for dual cone search
   \param apex_angle angle for dual cone search
   \param LfsMap the map to store the LFS value
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \return The estimated local feature size is stored in `lfs_map`.
*/
template <typename ConcurrencyTag,
          typename PointRange,
          typename LfsMap,
          typename NamedParameters = parameters::Default_named_parameters>
void
estimate_local_feature_size(PointRange& points,
                            const unsigned int jet_k,
                            const std::size_t N_rays,
                            const typename Point_set_processing_3_np_helper<PointRange, NamedParameters>::Geom_traits::FT apex_angle,
                            LfsMap lfs_map,
                            const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  CGAL_TRACE_STREAM << "Calls estimate_local_feature_size()\n";

  // basic geometric types
  typedef Point_set_processing_3_np_helper<PointRange, NamedParameters> NP_helper;
  typedef typename NP_helper::Point_map PointMap;
  typedef typename NP_helper::Normal_map NormalMap;
  typedef typename NP_helper::Geom_traits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Sphere_3 Sphere;

  CGAL_assertion_msg(NP_helper::has_normal_map(points, np), "Error: no normal map");
  PointMap point_map = NP_helper::get_point_map(points, np);
  NormalMap normal_map = NP_helper::get_normal_map(points, np);
  unsigned int degree_fitting = choose_parameter(get_parameter(np, internal_np::degree_fitting), 2);
  unsigned int degree_monge = choose_parameter(get_parameter(np, internal_np::degree_monge), 2);
  FT neighbor_radius = choose_parameter(get_parameter(np, internal_np::neighbor_radius), FT(0.0));

  const std::function<bool(double)>& callback = choose_parameter(get_parameter(np, internal_np::callback),
                                                               std::function<bool(double)>());

  // types for K nearest neighbors search structure
  typedef Point_set_processing_3::internal::Neighbor_query<Kernel, PointRange&, PointMap> Neighbor_query;

  CGAL_precondition(points.begin() != points.end());

  // precondition: at least 2 nearest neighbors
  // @todo fix this k as it does not exist
  // CGAL_precondition(k >= 2);

  std::size_t memory = CGAL::Memory_sizer().virtual_size();
  CGAL_TRACE_STREAM << (memory >> 20) << " Mb allocated\n";
  CGAL_TRACE_STREAM << "  Creates KD-tree\n";

  Neighbor_query neighbor_query(points, point_map);

  // Input points types
  typedef typename PointRange::iterator iterator;
  typedef typename iterator::value_type value_type;

  std::size_t nb_points = points.size();

  // calculate a loose bounding shpere
  Sphere bsphere = CGAL::internal::calculate_loose_bsphere<Kernel>(points, point_map);
  // estimate an epsilon band for the classical distance function
  FT epsilon_band = CGAL::internal::classical_point_dist_func_epsilon_band(points, neighbor_query, point_map);

  Point_set_processing_3::internal::Callback_wrapper<ConcurrencyTag>
  callback_wrapper (callback, nb_points);

  CGAL::for_each<ConcurrencyTag>
  (points,
   [&](value_type& vt)
   {
    if (callback_wrapper.interrupted())
      return false;

    const Point& point = get(point_map, vt);
    Vector normal = get(normal_map, vt);

    bool need_jet_normal = false;

    FT squared_length = std::sqrt(normal.squared_length());
    if (squared_length == 0.0)
    {
      need_jet_normal = true;
    }
    else  if (std::abs(squared_length - 1.0) > 1e-5) // normalize to unit vector
    {
      normal = normal / std::sqrt(squared_length);
    }

    FT lfs = CGAL::internal::estimate_local_feature_size
          (point, normal, neighbor_query, bsphere,
            jet_k, neighbor_radius, degree_fitting, degree_monge,
            epsilon_band, apex_angle, N_rays);

    put(lfs_map, vt, lfs);

    if (need_jet_normal)
      put(normal_map, vt, normal);

    ++ callback_wrapper.advancement();

     return true;
   });
}

/**
   \ingroup PkgPointSetProcessing3Algorithms

   Smooth the local feature size using median filter.


   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range
   \param knn number of neighbors for median filter
   \param T number of iterations
   \param LfsMap the map for the LFS value
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \return The function will smooth the local feature size values directly in `lfs_map`.
*/
template <typename ConcurrencyTag,
          typename PointRange,
          typename LfsMap,
          typename NamedParameters = parameters::Default_named_parameters>
void
median_filter_smoothing_lfs(PointRange& points,
                        const unsigned int knn,
                        const unsigned int T,
                        LfsMap lfs_map,
                        const NamedParameters& np = parameters::default_values())
{
  // basic geometric types
  typedef Point_set_processing_3_np_helper<PointRange, NamedParameters> NP_helper;
  typedef typename NP_helper::Geom_traits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename PointRange::iterator iterator;
  typedef typename iterator::value_type value_type;

  CGAL_precondition(points.begin() != points.end());

  std::vector<std::vector<iterator>> knn_graph;
  CGAL::internal::construct_knn_graph<ConcurrencyTag>(points, knn, knn_graph, np);

  std::size_t nb_points = points.size();

  assert(nb_points == knn_graph.size());

  auto begin = boost::make_zip_iterator (boost::make_tuple (points.begin(), knn_graph.begin()));
  auto end = boost::make_zip_iterator (boost::make_tuple (points.end(), knn_graph.end()));

  for (int t = 0; t < T; t++)
  {
    for(auto it = begin; it != end; it++)
    {
      value_type& vt = boost::get<0>(*it);
      FT lfs = get(lfs_map, vt);
      const std::vector<iterator>& iterators = boost::get<1>(*it);

      std::vector<FT> knn_lfs_vals;
      for (const auto& iterator : iterators) knn_lfs_vals.push_back(get(lfs_map, *iterator));
      const auto median = knn_lfs_vals.begin() + knn_lfs_vals.size() / 2;
      std::nth_element(knn_lfs_vals.begin(), median, knn_lfs_vals.end());

      put(lfs_map, vt, *median);
    }
  }
}

/**
   \ingroup PkgPointSetProcessing3Algorithms

   Smooth the local feature size based on the lipschitz continuity.


   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range
   \param knn number of neighbors for the smoothing based on the lipschitz continuity, suggested set it to be 1.0
   \param LfsMap the map for the LFS value
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \return The function will smooth the local feature size values directly in `lfs_map`.
*/
template <typename ConcurrencyTag,
          typename PointRange,
          typename LfsMap,
          typename NamedParameters = parameters::Default_named_parameters>
void
lipschitz_continuity_smoothing_lfs(PointRange& points,
                        const unsigned int knn,
                        const typename Point_set_processing_3_np_helper<PointRange, NamedParameters>::Geom_traits::FT lipschitz,
                        LfsMap lfs_map,
                        const NamedParameters& np = parameters::default_values())
{
  // basic geometric types
  typedef Point_set_processing_3_np_helper<PointRange, NamedParameters> NP_helper;
  typedef typename NP_helper::Point_map PointMap;
  typedef typename NP_helper::Geom_traits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename PointRange::iterator iterator;
  typedef typename iterator::value_type value_type;

  CGAL_precondition(points.begin() != points.end());

  std::vector<std::vector<iterator>> knn_graph;
  CGAL::internal::construct_knn_graph<ConcurrencyTag>(points, knn, knn_graph, np);

  std::size_t nb_points = points.size();

  assert(nb_points == knn_graph.size());

  PointMap point_map = NP_helper::get_point_map(points, np);

  // Property map typename PointRange::iterator -> index
  typedef Index_property_map<iterator> IndexMap;
  IndexMap index_map(points.begin(), points.end());

  // starting from the minimum of the func_vals
  auto min_iterator = std::min_element(points.begin(), points.end(),
        [&](const value_type& vt1, const value_type& vt2) {
            return get(lfs_map, vt1) < get(lfs_map, vt2);
        }
  );
  std::size_t min_lfs_index = std::distance(std::begin(points), min_iterator);

  std::unordered_map<std::size_t, int> visited_map;
  visited_map[min_lfs_index] = 1;

  std::function<bool(std::size_t, std::size_t)> cmp =
        [&](std::size_t ind1, std::size_t ind2)
    { return get(lfs_map, *(points.begin() + ind1)) > get(lfs_map, *(points.begin() + ind2)); };
  std::priority_queue<std::size_t, std::vector<std::size_t>, decltype(cmp)> min_top_queue(cmp);

  for (iterator it = points.begin(); it != points.end(); it++) min_top_queue.push(get(index_map, it));

  while (!min_top_queue.empty())
  {
    std::size_t cur = min_top_queue.top();
    min_top_queue.pop();

    // check if the current node satisfying lipschitz continuity
    const std::vector<iterator>& iterators = *(knn_graph.begin() + cur);
    auto iterator_cur = points.begin() + cur;
    FT cur_lfs = get(lfs_map, *iterator_cur);
    for (const auto& iterator : iterators)
    {
      std::size_t nei = get(index_map, iterator);
      auto iterator_nei = points.begin() + nei;

      if (nei == cur) continue;

      if (visited_map.find(nei) != visited_map.end())
          continue;
      else
          visited_map[nei] = 1;

      FT dist = CGAL::squared_distance(get(point_map, *iterator_cur), get(point_map, *iterator_nei));
      dist = std::sqrt(dist);

      FT nei_lfs = get(lfs_map, *iterator_nei);

      if (nei_lfs > cur_lfs)
      {
        if (std::abs(nei_lfs - cur_lfs) > lipschitz * dist)
          nei_lfs = cur_lfs + dist;
          put(lfs_map, *iterator_nei, nei_lfs);
      }
      else
      {
        min_top_queue.push(cur);
        break;
      }
    }
  }
}

/**
   \ingroup PkgPointSetProcessing3Algorithms

   Smooth the local feature size using median filter.


   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range
   \param knn number of neighbors for laplacian smoothing
   \param T number of iterations
   \param lambda lambda value for laplacian smoothing
   \param LfsMap the map for the LFS value
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \return The function will smooth the local feature size values directly in `lfs_map`.
*/
template <typename ConcurrencyTag,
          typename PointRange,
          typename LfsMap,
          typename NamedParameters = parameters::Default_named_parameters>
void
laplacian_smoothing_lfs(PointRange& points,
                        const unsigned int knn,
                        const unsigned int T,
                        const typename Point_set_processing_3_np_helper<PointRange, NamedParameters>::Geom_traits::FT lambda,
                        LfsMap lfs_map,
                        const NamedParameters& np = parameters::default_values())
{
  // basic geometric types
  typedef Point_set_processing_3_np_helper<PointRange, NamedParameters> NP_helper;
  typedef typename NP_helper::Geom_traits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename PointRange::iterator iterator;
  typedef typename iterator::value_type value_type;

  CGAL_precondition(points.begin() != points.end());

  std::vector<std::vector<iterator>> knn_graph;
  CGAL::internal::construct_knn_graph<ConcurrencyTag>(points, knn, knn_graph, np);

  std::size_t nb_points = points.size();

  assert(nb_points == knn_graph.size());

  auto begin = boost::make_zip_iterator (boost::make_tuple (points.begin(), knn_graph.begin()));
  auto end = boost::make_zip_iterator (boost::make_tuple (points.end(), knn_graph.end()));

  LfsMap smoothed_lfs_map = lfs_map;
  for (int t = 0; t < T; t++)
  {
    for(auto it = begin; it != end; it++)
    {
      value_type& vt = boost::get<0>(*it);
      FT lfs = get(lfs_map, vt);
      const std::vector<iterator>& iterators = boost::get<1>(*it);

      FT avg_lfs = 0.0;
      for (const auto& iterator : iterators) avg_lfs += get(lfs_map, *iterator);
      avg_lfs /= iterators.size();

      FT smoothed_lfs = lfs + lambda * (avg_lfs - lfs);
      put(smoothed_lfs_map, vt, smoothed_lfs);
    }
  }

  std::swap(lfs_map, smoothed_lfs_map);
}

} //namespace CGAL

#endif // CGAL_ESTIMATE_LFS_H
