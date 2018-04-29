// Copyright (c) 2013-06  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
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

#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/algorithm.h>
#include <iterator>
#include <set>
#include <algorithm>
#include <cmath>
#include <ctime>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#endif // CGAL_LINKED_WITH_TBB

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
  Kd_tree_element(const Kd_tree_element& other)
    : Base(other), index(other.index)
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

#ifdef CGAL_LINKED_WITH_TBB
/// \cond SKIP_IN_MANUAL
/// This is for parallelization of function: compute_denoise_projection()
template <typename Kernel, typename Tree, typename RandomAccessIterator>
class Sample_point_updater 
{
  typedef typename Kernel::Point_3   Point;
  typedef typename Kernel::FT        FT;

  std::vector<Point> &update_sample_points;
  std::vector<Point> &sample_points;
  const Tree &original_kd_tree;            
  const Tree &sample_kd_tree;              
  const typename Kernel::FT radius;  
  const std::vector<typename Kernel::FT> &original_densities;
  const std::vector<typename Kernel::FT> &sample_densities; 

public:
  Sample_point_updater(
    std::vector<Point> &out,
    std::vector<Point> &in,
    const Tree &_original_kd_tree,            
    const Tree &_sample_kd_tree,              
    const typename Kernel::FT _radius,
    const std::vector<typename Kernel::FT> &_original_densities,
    const std::vector<typename Kernel::FT> &_sample_densities): 
  update_sample_points(out), 
    sample_points(in),
    original_kd_tree(_original_kd_tree),
    sample_kd_tree(_sample_kd_tree),
    radius(_radius),
    original_densities(_original_densities),
    sample_densities(_sample_densities){} 


  void operator() ( const tbb::blocked_range<size_t>& r ) const 
  { 
    for (size_t i = r.begin(); i != r.end(); ++i) 
    {
      update_sample_points[i] = simplify_and_regularize_internal::
        compute_update_sample_point<Kernel, Tree, RandomAccessIterator>(
        sample_points[i], 
        original_kd_tree,
        sample_kd_tree,
        radius, 
        original_densities,
        sample_densities);
    }
  }
};
/// \endcond  
#endif // CGAL_LINKED_WITH_TBB


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/**
   \ingroup PkgPointSetProcessingAlgorithms
   This is an implementation of the Weighted Locally Optimal Projection (WLOP) simplification algorithm.
   The WLOP simplification algorithm can produce a set of 
   denoised, outlier-free and evenly distributed particles over the original 
   dense point cloud. 
   The core of the algorithm is a Weighted Locally Optimal Projection operator
   with a density uniformization term. 
   For more details, please refer to \cgalCite{wlop-2009}.

   A parallel version of WLOP is provided and requires the executable to be 
   linked against the <a href="http://www.threadingbuildingblocks.org">Intel TBB library</a>.
   To control the number of threads used, the user may use the tbb::task_scheduler_init class.
   See the <a href="http://www.threadingbuildingblocks.org/documentation">TBB documentation</a> 
   for more details.

   \tparam ConcurrencyTag enables sequential versus parallel algorithm.
   Possible values are `Sequential_tag`
   and `Parallel_tag`.
   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.
   \tparam OutputIterator Type of the output iterator. 
   It must accept objects of type `geom_traits::Point_3`.

   \param points input point range.
   \param output iterator where output points are put.
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadWritePropertyMap` with value type `geom_traits::Point_3`.
     If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{normal_map} a model of `ReadWritePropertyMap` with value type
     `geom_traits::Vector_3`.\cgalParamEnd
     \cgalParamBegin{select_percentage} percentage of points to retain. The default value is set to 
     5 (\%).\cgalParamEnd
     \cgalParamBegin{neighbor_radius} spherical neighborhood radius. This is a key parameter that needs to be
     finely tuned. The result will be irregular if too small, but a larger value will impact the runtime. In 
     practice, choosing a radius such that the neighborhood of each sample point includes at least two rings 
     of neighboring sample points gives satisfactory result. The default value is set to 8 times the average 
     spacing of the point set.\cgalParamEnd
     \cgalParamBegin{number_of_iterations} number of iterations to solve the optimsation problem. The default
     value is 35. More iterations give a more regular result but increase the runtime.\cgalParamEnd
     \cgalParamBegin{require_uniform_sampling} an optional preprocessing, which will give better result if the
     distribution of the input points is highly non-uniform. The default value is `false`. \cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
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
  using boost::choose_param;
  
  // basic geometric types
  typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;

  PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
  double select_percentage = choose_param(get_param(np, internal_np::select_percentage), 5.);
  double radius = choose_param(get_param(np, internal_np::neighbor_radius), -1);
  unsigned int iter_number = choose_param(get_param(np, internal_np::number_of_iterations), 35);
  bool require_uniform_sampling = choose_param(get_param(np, internal_np::require_uniform_sampling), false);

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

  FT radius2 = radius * radius;
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
    

    typename std::vector<Point>::iterator update_iter = update_sample_points.begin();
#ifndef CGAL_LINKED_WITH_TBB
  CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
			     "Parallel_tag is enabled but TBB is unavailable.");
#else
    //parallel
    if (boost::is_convertible<ConcurrencyTag, Parallel_tag>::value)
    {
      tbb::blocked_range<size_t> block(0, number_of_sample);
      Sample_point_updater<Kernel, Kd_Tree, typename PointRange::iterator> sample_updater(
                           update_sample_points,
                           sample_points,          
                           original_kd_tree,
                           sample_kd_tree,
                           radius2,
                           original_density_weights,
                           sample_density_weights);

       tbb::parallel_for(block, sample_updater);
    }else
#endif
    {
      //sequential
      for (sample_iter = sample_points.begin();
        sample_iter != sample_points.end(); ++sample_iter, ++update_iter)
      {
        *update_iter = simplify_and_regularize_internal::
          compute_update_sample_point<Kernel,
                                      Kd_Tree,
                                      typename PointRange::iterator>
                                      (*sample_iter,
                                       original_kd_tree,
                                       sample_kd_tree,
                                       radius2,
                                       original_density_weights,
                                       sample_density_weights);
      }
    }
    
    sample_iter = sample_points.begin();
    for (update_iter = update_sample_points.begin();
         update_iter != update_sample_points.end();
         ++update_iter, ++sample_iter)
    {
      *sample_iter = *update_iter;
    }
  }

  // final output
  for(sample_iter = sample_points.begin(); 
      sample_iter != sample_points.end(); ++sample_iter)
  {
    *output++ = *sample_iter;
  }

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

#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated API
template <typename ConcurrencyTag,
          typename OutputIterator,
          typename RandomAccessIterator,
          typename PointMap,
          typename Kernel>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::wlop_simplify_and_regularize_point_set(), please update your code")
OutputIterator
wlop_simplify_and_regularize_point_set(
  RandomAccessIterator first,  ///< random-access iterator to the first input point.
  RandomAccessIterator beyond, ///< past-the-end iterator.
  OutputIterator output,       ///< output iterator where output points are put.
  PointMap point_map,        ///< point property map.
  double select_percentage,    ///< percentage of points to retain. 
                               ///< The default value is set to 5 (\%).
  double radius,               ///< spherical neighborhood radius.
                               ///< This is a key parameter that needs to be finely tuned.  
                               ///< The result will be irregular if too small, but a larger
                               ///< value will impact the runtime.
                               ///< In practice, choosing a radius such that the neighborhood of each sample point
                               ///< includes at least two rings of neighboring sample points
                               ///< gives satisfactory result.
                               ///< The default value is set to 8 times the average spacing of the point set.
  unsigned int iter_number,    ///< number of iterations to solve the optimsation problem. The default value is 35.
                               ///< More iterations give a more regular result but increase the runtime.
  bool require_uniform_sampling,///< an optional preprocessing, which will give better result
                               ///< if the distribution of the input points is highly non-uniform. 
                               ///< The default value is `false`. 
  const Kernel&                ///< geometric traits.
)
{
  CGAL::Iterator_range<RandomAccessIterator> points (first, beyond);
  return wlop_simplify_and_regularize_point_set<ConcurrencyTag>
    (points, output,
     CGAL::parameters::point_map (point_map).
     select_percentage (select_percentage).
     neighbor_radius (radius).
     number_of_iterations(iter_number).
     require_uniform_sampling (require_uniform_sampling).
     geom_traits (Kernel()));
}
  

// deprecated API
template <typename ConcurrencyTag,
          typename OutputIterator,     
          typename RandomAccessIterator, 
          typename PointMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::wlop_simplify_and_regularize_point_set(), please update your code")
OutputIterator 
wlop_simplify_and_regularize_point_set(
  RandomAccessIterator  first,  ///< iterator over the first input point
  RandomAccessIterator  beyond, ///< past-the-end iterator
  OutputIterator output,        ///< add back-inserter
  PointMap point_map, ///< property map RandomAccessIterator  -> Point_3
  const double select_percentage,     ///< percentage of points to retain
  double neighbor_radius,       ///< size of neighbors.
  const unsigned int max_iter_number, ///< number of iterations.
  const bool require_uniform_sampling     ///< if needed to compute density 
                                      ///  to generate more rugularized result.                                 
) 
{
  CGAL::Iterator_range<RandomAccessIterator> points (first, beyond);
  return wlop_simplify_and_regularize_point_set<ConcurrencyTag>
    (points, output,
     CGAL::parameters::point_map (point_map).
     select_percentage (select_percentage).
     neighbor_radius (neighbor_radius).
     number_of_iterations(max_iter_number).
     require_uniform_sampling (require_uniform_sampling));
}

// deprecated API
template <typename ConcurrencyTag, 
          typename OutputIterator,     
          typename RandomAccessIterator >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::wlop_simplify_and_regularize_point_set(), please update your code")
OutputIterator
wlop_simplify_and_regularize_point_set(
  RandomAccessIterator  first,  ///< iterator to the first input point.
  RandomAccessIterator  beyond, ///< past-the-end iterator.
  OutputIterator output,        ///< add back-inserter.
  const double select_percentage = 5, ///< percentage of points to retain
  double neighbor_radius = -1,  ///< size of neighbors.
  const unsigned int max_iter_number = 35, ///< number of iterations.
  const bool require_uniform_sampling = false ///< if needed to compute density   
                                           ///to generate a more uniform result. 
)
{
  CGAL::Iterator_range<RandomAccessIterator> points (first, beyond);
  return wlop_simplify_and_regularize_point_set<ConcurrencyTag>
    (points, output,
     CGAL::parameters::select_percentage (select_percentage).
     neighbor_radius (neighbor_radius).
     number_of_iterations(max_iter_number).
     require_uniform_sampling (require_uniform_sampling));
}
#endif // CGAL_NO_DEPRECATED_CODE
/// @endcond

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_wlop_simplify_and_regularize_point_set_H
