// Copyright (c) 2013-06  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Shihao Wu, Clement Jamin, Pierre Alliez

#ifndef CGAL_UPSAMPLE_POINT_SET_H
#define CGAL_UPSAMPLE_POINT_SET_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Point_set_processing_3/internal/Rich_grid.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/compute_average_spacing.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iterator>
#include <set>
#include <utility>


//#define  CGAL_PSP3_VERBOSE

namespace CGAL {

// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL

namespace upsample_internal{

/// For each query point, select a best "base point" in its neighborhoods.
/// Then, a new point will be interpolated between query point and "base point".
/// This is the key part of the upsample algorithm
///
/// \pre `radius > 0`
///
/// @tparam Kernel Geometric traits class.
///
/// @return local density length
template <typename Kernel>
typename Kernel::FT
base_point_selection(
  const rich_grid_internal::Rich_point<Kernel>& query, ///< 3D point to project
  const std::vector<rich_grid_internal::Rich_point<Kernel> >&
                    neighbor_points,///< neighbor sample points
  const typename Kernel::FT edge_sensitivity,///< edge sensitivity parameter
  unsigned int& output_base_index ///< base point index
  )
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;
  typedef typename rich_grid_internal::Rich_point<Kernel> Rich_point;

  if (neighbor_points.empty())
  {
#ifdef CGAL_PSP3_VERBOSE
    std::cout << "empty neighborhood" << std::endl;
#endif
    output_base_index = query.index;
    return 0.0;
  }

  FT best_dist2 = -10.0;
  const Rich_point& v = query;
  typename std::vector<Rich_point>::const_iterator iter = neighbor_points.begin();
  for (; iter != neighbor_points.end(); ++iter)
  {
    const Point& t = iter->pt;

    const Vector& vm = v.normal;
    const Vector& tm = iter->normal;

    Vector diff_v_t = t - v.pt;
    Point mid_point = v.pt + (diff_v_t * FT(0.5));

    FT dot_produce = std::pow((FT(2.0) - vm * tm), edge_sensitivity);

    Vector diff_t_mid = mid_point - t;
    FT project_t = diff_t_mid * tm;
    FT min_dist2 = diff_t_mid.squared_length() - project_t * project_t;

    typename std::vector<Rich_point>::const_iterator iter_in = neighbor_points.begin();
    for (; iter_in != neighbor_points.end(); ++iter_in)
    {
      Vector diff_s_mid = mid_point - iter_in->pt;
      FT project_s = diff_s_mid * iter_in->normal;

      FT proj_min2 = diff_s_mid.squared_length() - project_s * project_s;

      if (proj_min2 < min_dist2)
      {
        min_dist2 = proj_min2;
      }
    }
    min_dist2 *= dot_produce;

    if (min_dist2 > best_dist2)
    {
      best_dist2 = min_dist2;
      output_base_index = iter->index;
    }
  }

  return best_dist2;
}

/// For each new inserted point, we need to do the following job
/// 1, get neighbor information from the two "parent points"
/// 2, update position and determine normal by bilateral projection
/// 3, update neighbor information again
///
/// \pre `radius > 0`
///
/// @tparam Kernel Geometric traits class.
///
template <typename Kernel>
void
update_new_point(
  unsigned int new_point_index, ///< new inserted point
  unsigned int father_index, ///< father point index
  unsigned int mother_index, ///< mother point index
  std::vector<rich_grid_internal::Rich_point<Kernel> >& rich_point_set,
                                                           ///< all rich points
  const typename Kernel::FT radius,          ///< accept neighborhood radius
  const typename Kernel::FT sharpness_bandwidth  ///< control sharpness
)
{
  // basic geometric types
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;
  typedef typename rich_grid_internal::Rich_point<Kernel> Rich_point;

  CGAL_assertion_code( unsigned int size = static_cast<unsigned int>(rich_point_set.size()) );
  CGAL_point_set_processing_precondition(father_index < size);
  CGAL_point_set_processing_precondition(mother_index < size);

  // 1, get neighbor information from the two "parent points"
  Rich_point& new_v = rich_point_set[new_point_index];
  Rich_point& father_v = rich_point_set[father_index];
  Rich_point& mother_v = rich_point_set[mother_index];

  std::set<int> neighbor_indexes;
  std::vector<unsigned int>::iterator iter;
  for (iter = father_v.neighbors.begin();
       iter != father_v.neighbors.end();
       ++iter)
  {
    neighbor_indexes.insert(*iter);
  }

  for (iter = mother_v.neighbors.begin();
       iter != mother_v.neighbors.end();
       ++iter)
  {
    neighbor_indexes.insert(*iter);
  }

  neighbor_indexes.insert(father_v.index);
  neighbor_indexes.insert(mother_v.index);

  FT radius2 = radius * radius;

  new_v.neighbors.clear();
  std::set<int>::iterator set_iter;
  for (set_iter = neighbor_indexes.begin();
       set_iter != neighbor_indexes.end(); ++set_iter)
  {
    Rich_point& t = rich_point_set[*set_iter];
    FT dist2 =  CGAL::squared_distance(new_v.pt, t.pt);

    if (dist2 < radius2)
    {
      new_v.neighbors.push_back(t.index);
    }
  }

  // 2, update position and normal by bilateral projection
  const unsigned int candidate_num = 2; // we have two normal candidates:
                                        // we say father's is 0
                                        //        mother's is 1
  std::vector<Vector> normal_cadidate(candidate_num);
  normal_cadidate[0] = father_v.normal;
  normal_cadidate[1] = mother_v.normal;

  std::vector<FT> project_dist_sum(candidate_num, FT(0.0));
  std::vector<FT> weight_sum(candidate_num, FT(0.0));
  std::vector<Vector> normal_sum(candidate_num, NULL_VECTOR);

  FT radius16 = FT(-4.0) / radius2;

  for (unsigned int i = 0; i < new_v.neighbors.size(); ++i)
  {
    const Rich_point& t = rich_point_set[new_v.neighbors[i]];
    FT dist2 = CGAL::squared_distance(new_v.pt, t.pt);
    FT theta = std::exp(dist2 * radius16);

    for (unsigned int j = 0; j < candidate_num; j++)
    {
      FT psi = std::exp(-std::pow(1 - normal_cadidate[j] * t.normal, 2)
                       / sharpness_bandwidth);
      FT project_diff_t_v = (t.pt - new_v.pt) * t.normal;
      FT weight = psi * theta;

      project_dist_sum[j] += project_diff_t_v * weight;
      normal_sum[j] = normal_sum[j] + t.normal * weight;
      weight_sum[j] += weight;
    }
  }

  // select best candidate
  FT min_project_dist = (std::numeric_limits<FT>::max)();
  unsigned int best = 0;

  for (unsigned int i = 0; i < candidate_num; ++i)
  {
    FT absolute_dist = CGAL::abs(project_dist_sum[i] / weight_sum[i]);
    if (absolute_dist < min_project_dist)
    {
      min_project_dist = absolute_dist;
      best = i;
    }
  }

  // update position and normal
  Vector update_normal = normal_sum[best] / weight_sum[best];
  new_v.normal = update_normal / sqrt(update_normal.squared_length());

  FT project_dist = project_dist_sum[best] / weight_sum[best];
  new_v.pt = new_v.pt + new_v.normal * project_dist;


  // 3, update neighbor information again
  new_v.neighbors.clear();
  for (set_iter = neighbor_indexes.begin();
       set_iter != neighbor_indexes.end(); ++set_iter)
  {
    Rich_point& t = rich_point_set[*set_iter];
    FT dist2 =  CGAL::squared_distance(new_v.pt, t.pt);

    if (dist2 < radius2)
    {
      new_v.neighbors.push_back(t.index);
      t.neighbors.push_back(new_v.index);
    }
  }
}

} /* namespace upsample_internal */

/// \endcond

// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/**
   \ingroup PkgPointSetProcessing3Algorithms
   This method progressively upsamples the point set while
   approaching the edge singularities (detected by normal variation), which
   generates a denser point set from an input point set. This has applications
   in point-based rendering, hole filling, and sparse surface reconstruction.
   Normals of points are required as input. For more details, please refer to \cgalCite{ear-2013}.

   \tparam ConcurrencyTag enables sequential versus parallel versions
   of `compute_average_spacing()` (called internally). Possible
   values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.
   \tparam OutputIterator Type of the output iterator.
   The type of the objects put in it is
   `std::pair<geom_traits::Point_3, geom_traits::Vector_3>`.
   Note that the user may use a
   <A HREF="https://www.boost.org/libs/iterator/doc/function_output_iterator.html">function_output_iterator</A>
   to match specific needs.

   \param points input point range.
   \param output iterator where output points and normals are put.
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Vector_3`}
     \cgalParamNEnd

     \cgalParamNBegin{sharpness_angle}
       \cgalParamDescription{controls the sharpness of the result}
       \cgalParamType{floating scalar value}
       \cgalParamDefault{`30.00`}
       \cgalParamExtra{The larger the value is, the smoother the result will be.
                       The range of possible value is `[0, 90]`}
     \cgalParamNEnd

     \cgalParamNBegin{edge_sensitivity}
       \cgalParamDescription{controls the priority of points inserted along sharp features}
       \cgalParamType{floating scalar value}
       \cgalParamDefault{`1`}
       \cgalParamExtra{Larger values of edge-sensitivity give higher priority to inserting points
                       along sharp features. The range of possible values is `[0, 1]`.
                       See section \ref Point_set_processing_3Upsample_Parameter1 for an example}
     \cgalParamNEnd

     \cgalParamNBegin{number_of_output_points}
       \cgalParamDescription{the number of output points to generate}
       \cgalParamType{unsigned int}
       \cgalParamDefault{`1000`}
     \cgalParamNEnd

     \cgalParamNBegin{neighbor_radius}
       \cgalParamDescription{the spherical neighborhood radius}
       \cgalParamType{floating scalar value}
       \cgalParamDefault{`0` (no limit)}
       \cgalParamExtra{If provided, the neighborhood of a query point is computed with a fixed spherical
                       radius instead of a fixed number of neighbors. In that case, the parameter
                       `k` is used as a limit on the number of points returned by each spherical
                       query (to avoid overly large number of points in high density areas).}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd
*/
template <typename ConcurrencyTag,
          typename PointRange,
          typename OutputIterator,
          typename NamedParameters>
OutputIterator
edge_aware_upsample_point_set(
  const PointRange& points,
  OutputIterator output,
  const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // basic geometric types
  typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::type NormalMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;

  CGAL_static_assertion_msg(!(boost::is_same<NormalMap,
                              typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::NoMap>::value),
                            "Error: no normal map");

  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;
  typedef typename rich_grid_internal::Rich_point<Kernel> Rich_point;

  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));
  NormalMap normal_map = choose_parameter<NormalMap>(get_parameter(np, internal_np::normal_map));
  double sharpness_angle = choose_parameter(get_parameter(np, internal_np::sharpness_angle), 30.);
  double edge_sensitivity = choose_parameter(get_parameter(np, internal_np::edge_sensitivity), 1);
  double neighbor_radius = choose_parameter(get_parameter(np, internal_np::neighbor_radius), -1);
  std::size_t number_of_output_points = choose_parameter(get_parameter(np, internal_np::number_of_output_points), 1000);

  std::cerr << sharpness_angle << " " << edge_sensitivity << " " << neighbor_radius
            << " " << number_of_output_points << std::endl;
  // trick in case the output iterator add points to the input container
  typename PointRange::const_iterator begin = points.begin();
  typename PointRange::const_iterator end = points.end();

  // preconditions
  CGAL_point_set_processing_precondition(begin != end);
  CGAL_point_set_processing_precondition(sharpness_angle >= 0
                                       &&sharpness_angle <= 90);
  CGAL_point_set_processing_precondition(edge_sensitivity >= 0
                                       &&edge_sensitivity <= 1);

  edge_sensitivity *= 10;  // just project [0, 1] to [0, 10].

  std::size_t number_of_input = std::distance(begin, end);
  CGAL_point_set_processing_precondition(number_of_output_points > number_of_input);


  const unsigned int nb_neighbors = 6; // 1 ring
  FT average_spacing = CGAL::compute_average_spacing<ConcurrencyTag>(
    points, nb_neighbors, np);

  if (neighbor_radius < average_spacing)
  {
    neighbor_radius = average_spacing * 3.0f;
#ifdef CGAL_PSP3_VERBOSE
    std::cout << "neighbor radius: " << neighbor_radius << std::endl;
#endif
  }

  Real_timer task_timer;

  // copy rich point set
  std::vector<Rich_point> rich_point_set(number_of_input);
  CGAL::Bbox_3 bbox;

  typename PointRange::const_iterator it = begin; // point iterator
  for(unsigned int i = 0; it != end; ++it, ++i)
  {
    rich_point_set[i].pt = get(point_map, *it);
    rich_point_set[i].normal = get(normal_map, *it);

    rich_point_set[i].index = i;
    bbox += rich_point_set[i].pt.bbox();
    CGAL_point_set_processing_precondition(rich_point_set[i].normal.squared_length() > 1e-10);
  }

  // compute neighborhood
  rich_grid_internal::compute_ball_neighbors_one_self(rich_point_set,
                                                      bbox,
                                                      FT(neighbor_radius));

  //
  FT cos_sigma = static_cast<FT>(std::cos(CGAL::to_double(sharpness_angle) / 180.0 * CGAL_PI));
  FT sharpness_bandwidth = std::pow((CGAL::max)((FT)1e-8, (FT)1.0 - cos_sigma), 2);

  FT sum_density = 0.0;
  unsigned int count_density = 1;
  double max_iter_time = 20;
  FT current_radius = FT(neighbor_radius);
  FT density_pass_threshold = 0.0;

  for (unsigned int iter_time = 0; iter_time < max_iter_time; ++iter_time)
  {
  #ifdef CGAL_PSP3_VERBOSE
     std::cout << std::endl << "iter_time: " << iter_time + 1  << std::endl;
  #endif
    if (iter_time > 0)
    {
      current_radius *= 0.75;
      if (current_radius < density_pass_threshold * 3) //3 ring
      {
        current_radius = density_pass_threshold * 3;
      }
      rich_grid_internal::compute_ball_neighbors_one_self(rich_point_set,
                                                          bbox,
                                                          current_radius);
    }
 #ifdef CGAL_PSP3_VERBOSE
    std::cout << "current radius: " << current_radius << std::endl;
 #endif

    std::size_t current_size = rich_point_set.size();
    std::vector<bool> is_pass_threshold(current_size, false);

    if (iter_time == 0)
    {
      //estimate density threshold for the first time
      for (unsigned int i = 0; i < rich_point_set.size() * 0.05; ++i)
      {
        const Rich_point& v = rich_point_set[i];

        if (v.neighbors.empty())
          continue;

        // extract neighbor rich points by index
        std::vector<Rich_point> neighbor_rich_points(v.neighbors.size());
        for (unsigned int n = 0; n < v.neighbors.size(); n++)
        {
          neighbor_rich_points[n] = rich_point_set[v.neighbors[n]];
        }

        unsigned int base_index = 0;
        FT density2 = upsample_internal::
                              base_point_selection(v,
                                                   neighbor_rich_points,
                                                   FT(edge_sensitivity),
                                                   base_index);

        if (density2 < 0)
        {
          continue;
        }

        sum_density += density2;
        count_density++;
      }
    }

    density_pass_threshold = static_cast<FT>(sqrt(sum_density / count_density) * FT(0.65));
    sum_density = 0.;
    count_density = 1;

    FT density_pass_threshold2 = density_pass_threshold *
                                 density_pass_threshold;
 #ifdef CGAL_PSP3_VERBOSE
    std::cout << "pass_threshold:  " << density_pass_threshold << std::endl;
 #endif
    // insert new points until all the points' density pass the threshold
    unsigned int max_loop_time = 3;
    unsigned int loop = 0;
    while (true)
    {
   #ifdef CGAL_PSP3_VERBOSE
      std::cout << "loop_time: " << loop + 1 << std::endl;
   #endif
      unsigned int count_not_pass = 0;
      loop++;
      for (unsigned int i = 0; i < rich_point_set.size(); ++i)
      {
        if (is_pass_threshold[i])
        {
          continue;
        }

        const Rich_point& v = rich_point_set[i];

        if (v.neighbors.empty())
          continue;

        // extract neighbor rich points by index
        std::vector<Rich_point> neighbor_rich_points(v.neighbors.size());
        for (unsigned int n = 0; n < v.neighbors.size(); ++n)
        {
          neighbor_rich_points[n] = rich_point_set[v.neighbors[n]];
        }

        // select base point
        unsigned int base_index = 0;
        FT density2 = upsample_internal::
                              base_point_selection(v,
                                                   neighbor_rich_points,
                                                   FT(edge_sensitivity),
                                                   base_index);

        // test if it pass the density threshold
        if (density2 < density_pass_threshold2)
        {
          is_pass_threshold[i] = true;
          continue;
        }
        count_not_pass++;

        sum_density += density2;
        count_density++;

        // insert a new rich point
        unsigned int father_index = v.index;
        unsigned int mother_index = base_index;

        Rich_point new_v;
        Rich_point& base = rich_point_set[mother_index];
        Vector diff_v_base = base.pt - v.pt;
        new_v.pt = v.pt + (diff_v_base * FT(0.5));
        new_v.index = static_cast<unsigned int>(rich_point_set.size());

        unsigned int new_point_index = new_v.index;
        rich_point_set.push_back(new_v);
        is_pass_threshold.push_back(false);

        //update new rich point
        upsample_internal::update_new_point(new_point_index,
                                            father_index,
                                            mother_index,
                                            rich_point_set,
                                            current_radius,
                                            sharpness_bandwidth);

        if (rich_point_set.size() >= number_of_output_points)
        {
          break;
        }
      }
   #ifdef CGAL_PSP3_VERBOSE
      std::cout << "current size: " << rich_point_set.size() << std::endl;
   #endif
      if (count_not_pass == 0 ||
          loop >= max_loop_time ||
          rich_point_set.size() >= number_of_output_points)
      {
        break;
      }

    }

    if (rich_point_set.size() >= number_of_output_points)
    {
      break;
    }
  }

  for (std::size_t i = number_of_input; i < rich_point_set.size(); ++i)
  {
    const Rich_point& v = rich_point_set[i];
    Point point = v.pt;
    Vector normal = v.normal;
    *output++ = std::make_pair(point, normal);
  }

  return output;
}


/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename ConcurrencyTag,
          typename PointRange,
          typename OutputIterator>
OutputIterator
edge_aware_upsample_point_set(
  const PointRange& points,
  OutputIterator output)
{
  return edge_aware_upsample_point_set<ConcurrencyTag>
    (points, output, CGAL::Point_set_processing_3::parameters::all_default(points));
}
/// \endcond

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_UPSAMPLE_POINT_SET_H
