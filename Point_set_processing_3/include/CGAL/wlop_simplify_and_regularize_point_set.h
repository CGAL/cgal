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

#ifndef CGAL_wlop_simplify_and_regularize_point_set_H
#define CGAL_wlop_simplify_and_regularize_point_set_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/compute_average_spacing.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/algorithm.h>
#include <iterator>
#include <set>
#include <algorithm>
#include <cmath>
#include <ctime>

#include <CGAL/Point_set_processing_3/internal/Callback_wrapper.h>
#include <CGAL/for_each.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Bbox_3.h>

namespace CGAL {
// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL

namespace simplify_and_regularize_internal{

// Item in the Kd-tree: position (Point_3) + index
template <typename Kernel>
class Kd_tree_element : public Kernel::Point_3
{
public:
  unsigned int index;

  // basic geometric types
  typedef typename CGAL::Origin Origin;
  typedef typename Kernel::Point_3 Base;

  Kd_tree_element(const Origin& o = ORIGIN, unsigned int id=0)
    : Base(o), index(id)
  {}
  Kd_tree_element(const Base& p, unsigned int id=0)
    : Base(p), index(id)
  {}

};

// Helper class for the Kd-tree
template <typename Kernel>
class Kd_tree_gt : public Kernel
{
public:
  typedef Kd_tree_element<Kernel> Point_3;
};

template <typename Kernel>
class Kd_tree_traits : public CGAL::Search_traits_3<Kd_tree_gt<Kernel> >
{
public:
  typedef typename Kernel::Point_3 PointType;
};

/// Compute average and repulsion term, then
/// compute and update sample point locations
///
/// \pre `radius > 0`
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree Kd-tree.
///
/// @return average term vector
template <typename Kernel,
          typename Tree,
          typename RandomAccessIterator>
typename Kernel::Point_3
compute_update_sample_point(
  const typename Kernel::Point_3& query, ///< 3D point to project
  const Tree& original_kd_tree,          ///< original Kd-tree
  const Tree& sample_kd_tree,            ///< sample Kd-tree
  const typename Kernel::FT radius,      ///< neighborhood radius square
  const std::vector<typename Kernel::FT>& original_densities, ///<
  const std::vector<typename Kernel::FT>& sample_densities ///<
)
{
  CGAL_point_set_processing_precondition(radius > 0);
  bool is_original_densities_empty = original_densities.empty();
  bool is_sample_densities_empty = sample_densities.empty();

  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;

  //types for range search
  typedef simplify_and_regularize_internal::Kd_tree_element<Kernel> Kd_tree_point;
  typedef simplify_and_regularize_internal::Kd_tree_traits<Kernel> Traits;
  typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

  //range search for original neighborhood
  Fuzzy_sphere fs(query, radius, 0.0);
  std::vector<Kd_tree_point> neighbor_original_points;
  original_kd_tree.search(std::back_inserter(neighbor_original_points), fs);

  //Compute average term
  FT radius2 = radius * radius;
  Vector average = CGAL::NULL_VECTOR;
  FT weight = (FT)0.0, average_weight_sum = (FT)0.0;
  FT iradius16 = -(FT)4.0 / radius2;

  typename std::vector<Kd_tree_point>::iterator iter;
  iter = neighbor_original_points.begin();
  for (; iter != neighbor_original_points.end(); ++iter)
  {
    const Point& np = *iter;

    Kd_tree_point& kp = *iter;
    int original_index = kp.index;

    FT dist2 = CGAL::squared_distance(query, np);
    if (dist2 < 1e-10) continue;

    weight = exp(dist2 * iradius16);

    if (!is_original_densities_empty)
    {
      weight *= original_densities[original_index];
    }
    average_weight_sum += weight;
    average = average + (np - CGAL::ORIGIN) * weight;
  }

  if (neighbor_original_points.empty() || average_weight_sum < FT(1e-10))
  {
    average = query - CGAL::ORIGIN;
  }
  else
  {
    average = average / average_weight_sum;
  }
  neighbor_original_points.clear();


  //Compute repulsion term

  Fuzzy_sphere fs2(query, radius, 0.0);
  std::vector<Kd_tree_point> neighbor_sample_points;
  sample_kd_tree.search(std::back_inserter(neighbor_sample_points), fs2);

  weight = (FT)0.0;
  FT repulsion_weight_sum = (FT)0.0;
  Vector repulsion = CGAL::NULL_VECTOR;

  iter = neighbor_sample_points.begin();
  for(; iter != neighbor_sample_points.end(); ++iter)
  {
    const Point& np = *iter;

    Kd_tree_point& kp = *iter;
    int sample_index = kp.index;

    FT dist2 = CGAL::squared_distance(query, np);
    if (dist2 < 1e-10) continue;
    FT dist = std::sqrt(dist2);

    weight = std::exp(dist2 * iradius16) * std::pow(FT(1.0) / dist, 2); // L1

    if (!is_sample_densities_empty)
    {
      weight *= sample_densities[sample_index];
    }

    Vector diff = query - np;

    repulsion_weight_sum += weight;
    repulsion = repulsion + diff * weight;
  }

  if (neighbor_sample_points.size() < 3 || repulsion_weight_sum < FT(1e-10))
  {
    repulsion = CGAL::NULL_VECTOR;
  }
  else
  {
    repulsion = repulsion / repulsion_weight_sum;
  }
  neighbor_sample_points.clear();

  // Compute update sample point
  Point update_sample = CGAL::ORIGIN + average + FT(0.45) * repulsion;
  return update_sample;
}


/// Compute density weight for each original points,
/// according to their neighbor original points
///
/// \pre `k >= 2`, radius > 0
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree Kd-tree.
///
/// @return computed point
template <typename Kernel, typename Tree>
typename Kernel::FT
compute_density_weight_for_original_point(
  const typename Kernel::Point_3& query, ///< 3D point to project
  Tree& original_kd_tree,                       ///< Kd-tree
  const typename Kernel::FT radius       ///< neighbor radius square
)
{
  CGAL_point_set_processing_precondition(radius > 0);

  // basic geometric types
  typedef typename Kernel::Point_3                         Point;
  typedef typename Kernel::FT                              FT;

  //types for range search
  typedef simplify_and_regularize_internal::Kd_tree_element<Kernel> Kd_tree_point;
  typedef simplify_and_regularize_internal::Kd_tree_traits<Kernel> Traits;
  typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

  //range search for original neighborhood
  Fuzzy_sphere fs(query, radius, 0.0);
  std::vector<Kd_tree_point> neighbor_original_points;

  original_kd_tree.search(std::back_inserter(neighbor_original_points), fs);

  //Compute density weight
  FT radius2 = radius * radius;
  FT density_weight = (FT)1.0;
  FT iradius16 = -(FT)4.0 / radius2;

  typename std::vector<Kd_tree_point>::iterator iter;
  iter = neighbor_original_points.begin();

  for (; iter != neighbor_original_points.end(); iter++)
  {
    const Point& np = *iter;

    FT dist2 = CGAL::squared_distance(query, np);
    if (dist2 < 1e-8) continue;

    density_weight += std::exp(dist2 * iradius16);
  }

  // output
  return FT(1.0) / density_weight;
}


/// Compute density weight for sample point,
/// according to their neighbor sample points
///
/// \pre `k >= 2`, radius > 0
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree Kd-tree.
///
/// @return computed point
template <typename Kernel, typename Tree>
typename Kernel::FT
compute_density_weight_for_sample_point(
  const typename Kernel::Point_3& query, ///< 3D point to project
  Tree& sample_kd_tree,                ///< Kd-tree
  const typename Kernel::FT radius       ///< neighbor radius square
)
{
  // basic geometric types
  typedef typename Kernel::Point_3                          Point;
  typedef typename Kernel::FT                               FT;

  //types for range search
  typedef simplify_and_regularize_internal::Kd_tree_element<Kernel> Kd_tree_point;
  typedef simplify_and_regularize_internal::Kd_tree_traits<Kernel> Traits;
  typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

  //range search for original neighborhood
  Fuzzy_sphere fs(query, radius, 0.0);
  std::vector<Kd_tree_point> neighbor_sample_points;
  sample_kd_tree.search(std::back_inserter(neighbor_sample_points), fs);

  //Compute density weight
  FT radius2 = radius * radius;
  FT density_weight = (FT)1.0;
  FT iradius16 = -(FT)4.0 / radius2;

  typename std::vector<Kd_tree_point>::iterator iter;
  iter = neighbor_sample_points.begin();

  for (; iter != neighbor_sample_points.end(); iter++)
  {
    const Point& np = *iter;

    FT dist2 = CGAL::squared_distance(query, np);
    density_weight += std::exp(dist2 * iradius16);
  }

  return density_weight;
}

} // namespace simplify_and_regularize_internal

/// \endcond


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/**
   \ingroup PkgPointSetProcessing3Algorithms
   This is an implementation of the Weighted Locally Optimal Projection (WLOP) simplification algorithm.
   The WLOP simplification algorithm can produce a set of
   denoised, outlier-free and evenly distributed particles over the original
   dense point cloud.
   The core of the algorithm is a Weighted Locally Optimal Projection operator
   with a density uniformization term.
   For more details, please refer to \cgalCite{wlop-2009}.

   A parallel version of WLOP is provided and requires the executable to be
   linked against the <a href="https://www.threadingbuildingblocks.org">Intel TBB library</a>.
   To control the number of threads used, the user may use the tbb::task_scheduler_init class.
   See the <a href="https://www.threadingbuildingblocks.org/documentation">TBB documentation</a>
   for more details.

   \tparam ConcurrencyTag enables sequential versus parallel algorithm. Possible values are `Sequential_tag`,
                          `Parallel_tag`, and `Parallel_if_available_tag`.
   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.
   \tparam OutputIterator Type of the output iterator.
   It must accept objects of type `geom_traits::Point_3`.

   \param points input point range.
   \param output iterator where output points are put.
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadWritePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point set `points`}
       \cgalParamType{a model of `ReadWritePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Vector_3`}
     \cgalParamNEnd

     \cgalParamNBegin{select_percentage}
       \cgalParamDescription{percentage of points to retain}
       \cgalParamType{floating scalar value}
       \cgalParamDefault{`5`}
     \cgalParamNEnd

     \cgalParamNBegin{neighbor_radius}
       \cgalParamDescription{the spherical neighborhood radius}
       \cgalParamType{floating scalar value}
       \cgalParamDefault{8 times the average spacing of the point set}
       \cgalParamExtra{This is a key parameter that needs to be finely tuned.
                       The result will be irregular if too small, but a larger value will impact the runtime.
                       In practice, choosing a radius such that the neighborhood of each sample point
                       includes at least two rings of neighboring sample points gives satisfactory result.}
     \cgalParamNEnd

     \cgalParamNBegin{number_of_iterations}
       \cgalParamDescription{number of iterations to solve the optimsation problem}
       \cgalParamType{unsigned int}
       \cgalParamDefault{`35`}
       \cgalParamExtra{More iterations give a more regular result but increase the runtime}
     \cgalParamNEnd

     \cgalParamNBegin{require_uniform_sampling}
       \cgalParamDescription{If `true`, an optional preprocessing is applied, which will give
                             better results if the distribution of the input points is highly non-uniform.}
       \cgalParamType{Boolean}
       \cgalParamDefault{`35`}
       \cgalParamExtra{More iterations give a more regular result but increase the runtime}
     \cgalParamNEnd

     \cgalParamNBegin{callback}
       \cgalParamDescription{a mechanism to get feedback on the advancement of the algorithm
                             while it's running and to interrupt it if needed}
       \cgalParamType{an instance of `std::function<bool(double)>`.}
       \cgalParamDefault{unused}
       \cgalParamExtra{It is called regularly when the
                       algorithm is running: the current advancement (between 0. and
                       1.) is passed as parameter. If it returns `true`, then the
                       algorithm continues its execution normally; if it returns
                       `false`, the algorithm is stopped, no output points are
                       generated.}
       \cgalParamExtra{The callback will be copied and therefore needs to be lightweight.}
       \cgalParamExtra{When `CGAL::Parallel_tag` is used, the `callback` mechanism is called asynchronously
                       on a separate thread and shouldn't access or modify the variables that are parameters of the algorithm.}
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
wlop_simplify_and_regularize_point_set(
  PointRange& points,
  OutputIterator output,
  const NamedParameters& np
)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // basic geometric types
  typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;

  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));
  double select_percentage = choose_parameter(get_parameter(np, internal_np::select_percentage), 5.);
  double radius = choose_parameter(get_parameter(np, internal_np::neighbor_radius), -1);
  unsigned int iter_number = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 35);
  bool require_uniform_sampling = choose_parameter(get_parameter(np, internal_np::require_uniform_sampling), false);
  const std::function<bool(double)>& callback = choose_parameter(get_parameter(np, internal_np::callback),
                                                                 std::function<bool(double)>());

  typedef typename Kernel::Point_3   Point;
  typedef typename Kernel::FT        FT;

  // types for K nearest neighbors search structure
  typedef simplify_and_regularize_internal::Kd_tree_element<Kernel> Kd_tree_element;
  typedef simplify_and_regularize_internal::Kd_tree_traits<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Kd_Tree;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(points.begin() != points.end());
  CGAL_point_set_processing_precondition(select_percentage >= 0
                                         && select_percentage <= 100);

  // Random shuffle
  CGAL::cpp98::random_shuffle (points.begin(), points.end());

  // Computes original(input) and sample points size
  std::size_t number_of_original = std::distance(points.begin(), points.end());
  std::size_t number_of_sample = (std::size_t)(FT(number_of_original) *
                                 (select_percentage / FT(100.0)));
  std::size_t first_index_to_sample = number_of_original - number_of_sample;

  // The first point iter of original and sample points
  typename PointRange::iterator it;                             // point iterator
  typename PointRange::iterator first_original_iter = points.begin();
  typename PointRange::iterator first_sample_iter = points.begin();
  std::advance(first_sample_iter, first_index_to_sample);

  //Copy sample points
  std::vector<Point> sample_points;
  sample_points.reserve(number_of_sample);
  unsigned int i;

  for(it = first_sample_iter; it != points.end(); ++it)
  {
    sample_points.push_back(get(point_map, *it));
  }

  //compute default neighbor_radius, if no radius in
  if (radius < 0)
  {
    const unsigned int nb_neighbors = 6; // 1 ring
    FT average_spacing = CGAL::compute_average_spacing<ConcurrencyTag>(points, nb_neighbors, np);
    radius = average_spacing * 8.0;

#ifdef CGAL_PSP3_VERBOSE
    std::cout << "The estimated radius size is: " << radius << std::endl;
    std::cout << "Be careful! Using this radius estimation may not be able to have good performance/result for different input" << std::endl;
#endif
  }

  CGAL_point_set_processing_precondition(radius > 0);

  // Initiate a KD-tree search for original points
  std::vector<Kd_tree_element> original_treeElements;
  for (it = first_original_iter, i=0 ; it != points.end() ; ++it, ++i)
    original_treeElements.push_back( Kd_tree_element(get(point_map, *it), i) );
  Kd_Tree original_kd_tree(original_treeElements.begin(),
                           original_treeElements.end());


  std::vector<Point> update_sample_points(number_of_sample);
  typename std::vector<Point>::iterator sample_iter;

  // Compute original density weight for original points if user needed
  std::vector<FT> original_density_weights;

  if (require_uniform_sampling)//default value is false
  {
    //todo: this part could also be parallelized if needed
    for (it = first_original_iter, i = 0; it != points.end() ; ++it, ++i)
    {
      FT density = simplify_and_regularize_internal::
                   compute_density_weight_for_original_point<Kernel, Kd_Tree>
                                         (
                                           get(point_map, *it),
                                           original_kd_tree,
                                           radius);

      original_density_weights.push_back(density);
    }
  }

  for (unsigned int iter_n = 0; iter_n < iter_number; ++iter_n)
  {
    // Initiate a KD-tree search for sample points
    std::vector<Kd_tree_element> sample_treeElements;
    for (i=0 ; i < sample_points.size(); i++)
    {
      Point& p0 = sample_points[i];
      sample_treeElements.push_back(Kd_tree_element(p0,i));
    }
    Kd_Tree sample_kd_tree(sample_treeElements.begin(), sample_treeElements.end());

    // Compute sample density weight for sample points
    std::vector<FT> sample_density_weights;

    for (sample_iter = sample_points.begin();
         sample_iter != sample_points.end(); ++sample_iter)
    {
      FT density = simplify_and_regularize_internal::
                   compute_density_weight_for_sample_point<Kernel, Kd_Tree>
                   (*sample_iter,
                    sample_kd_tree,
                    radius);

      sample_density_weights.push_back(density);
    }

    typedef boost::zip_iterator<boost::tuple<typename std::vector<Point>::iterator,
                                             typename std::vector<Point>::iterator> > Zip_iterator;

    Point_set_processing_3::internal::Callback_wrapper<ConcurrencyTag>
      callback_wrapper (callback, iter_number * number_of_sample, iter_n * number_of_sample);

    CGAL::for_each<ConcurrencyTag>
      (CGAL::make_range (boost::make_zip_iterator (boost::make_tuple (sample_points.begin(), update_sample_points.begin())),
                         boost::make_zip_iterator (boost::make_tuple (sample_points.end(), update_sample_points.end()))),
       [&](const typename Zip_iterator::reference& t)
       {
         if (callback_wrapper.interrupted())
           return false;

         get<1>(t) = simplify_and_regularize_internal::
           compute_update_sample_point<Kernel, Kd_Tree, typename PointRange::iterator>(
             get<0>(t),
             original_kd_tree,
             sample_kd_tree,
             radius,
             original_density_weights,
             sample_density_weights);
         ++ callback_wrapper.advancement();

         return true;
       });

    bool interrupted = callback_wrapper.interrupted();

    // We interrupt by hand as counter only goes halfway and won't terminate by itself
    callback_wrapper.interrupted() = true;
    callback_wrapper.join();

    // If interrupted during this step, nothing is computed, we return NaN
    if (interrupted)
      return output;

    sample_iter = sample_points.begin();
    for (std::size_t i = 0; i < sample_points.size(); ++ i)
      sample_points[i] = update_sample_points[i];
  }

  // final output
  std::copy (sample_points.begin(), sample_points.end(), output);

  return output;
}


/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename ConcurrencyTag,
          typename PointRange,
          typename OutputIterator>
OutputIterator
wlop_simplify_and_regularize_point_set(
  PointRange& points,
  OutputIterator output)       ///< output iterator where output points are put.
{
  return wlop_simplify_and_regularize_point_set<ConcurrencyTag>
    (points, output, CGAL::Point_set_processing_3::parameters::all_default(points));
}
/// \endcond

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_wlop_simplify_and_regularize_point_set_H
