// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Pierre Alliez and Laurent Saboret and Marc Pouget and Frederic Cazals

#ifndef CGAL_JET_ESTIMATE_NORMALS_H
#define CGAL_JET_ESTIMATE_NORMALS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/IO/trace.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Point_set_processing_3/internal/neighbor_query.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Memory_sizer.h>
#include <functional>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iterator>
#include <list>

#ifdef CGAL_LINKED_WITH_TBB
#include <CGAL/Point_set_processing_3/internal/Parallel_callback.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>  
#endif // CGAL_LINKED_WITH_TBB

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL
namespace internal {


/// Estimates normal direction using jet fitting
/// on the k nearest neighbors.
///
/// \pre `k >= 2`
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return Computed normal. Orientation is random.
template < typename Kernel,
           typename SvdTraits,
           typename Tree
>
typename Kernel::Vector_3
jet_estimate_normal(const typename Kernel::Point_3& query, ///< point to compute the normal at
                    Tree& tree, ///< KD-tree
                    unsigned int k, ///< number of neighbors
                    typename Kernel::FT neighbor_radius,
                    unsigned int degree_fitting)
{
  // basic geometric types
  typedef typename Kernel::Point_3  Point;

  // types for jet fitting
  typedef Monge_via_jet_fitting< Kernel,
                                 Simple_cartesian<double>,
                                 SvdTraits> Monge_jet_fitting;
  typedef typename Monge_jet_fitting::Monge_form Monge_form;

  std::vector<Point> points;
  CGAL::Point_set_processing_3::internal::neighbor_query
    (query, tree, k, neighbor_radius, points);
  
  // performs jet fitting
  Monge_jet_fitting monge_fit;
  const unsigned int degree_monge = 1; // we seek for normal and not more.
  Monge_form monge_form = monge_fit(points.begin(), points.end(),
                                    degree_fitting, degree_monge);

  // output normal vector (already normalized in monge form)
  return monge_form.normal_direction();
}

#ifdef CGAL_LINKED_WITH_TBB
  template <typename Kernel, typename SvdTraits, typename Tree>
  class Jet_estimate_normals {
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    const Tree& tree;
    const unsigned int k;
    const FT neighbor_radius;
    const unsigned int degree_fitting;
    const std::vector<Point>& input;
    std::vector<Vector>& output;
    cpp11::atomic<std::size_t>& advancement;
    cpp11::atomic<bool>& interrupted;

  public:
    Jet_estimate_normals(Tree& tree, unsigned int k, FT neighbor_radius,
                         std::vector<Point>& points,
                         unsigned int degree_fitting, std::vector<Vector>& output,
                         cpp11::atomic<std::size_t>& advancement,
                         cpp11::atomic<bool>& interrupted)
      : tree(tree), k (k), neighbor_radius (neighbor_radius)
      , degree_fitting (degree_fitting), input (points), output (output)
      , advancement (advancement)
      , interrupted (interrupted)
    { }
    
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for( std::size_t i = r.begin(); i != r.end(); ++i)
      {
        if (interrupted)
          break;
        output[i] = CGAL::internal::jet_estimate_normal<Kernel,SvdTraits>(input[i], tree, k, neighbor_radius, degree_fitting);
        ++ advancement;
      }
    }

  };
#endif // CGAL_LINKED_WITH_TBB


  
} /* namespace internal */
/// \endcond



// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/**
   \ingroup PkgPointSetProcessing3Algorithms
   Estimates normal directions of the range of `points`
   using jet fitting on the nearest neighbors.
   The output normals are randomly oriented.

   \pre `k >= 2`

   \tparam ConcurrencyTag enables sequential versus parallel algorithm. Possible values are `Sequential_tag`,
                          `Parallel_tag`, and `Parallel_if_available_tag`.
   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param k number of neighbors
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`.
     If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{normal_map} a model of `ReadWritePropertyMap` with value type
     `geom_traits::Vector_3`.\cgalParamEnd
     \cgalParamBegin{neighbor_radius} spherical neighborhood radius. If
     provided, the neighborhood of a query point is computed with a fixed spherical
     radius instead of a fixed number of neighbors. In that case, the parameter
     `k` is used as a limit on the number of points returned by each spherical
     query (to avoid overly large number of points in high density areas). If no
     limit is wanted, use `k=0`.\cgalParamEnd
     \cgalParamBegin{degree_fitting} degree of jet fitting.\cgalParamEnd
     \cgalParamBegin{svd_traits} template parameter for the class `Monge_via_jet_fitting`. If
     \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available and `CGAL_EIGEN3_ENABLED` is defined,
     then `CGAL::Eigen_svd` is used.\cgalParamEnd
     \cgalParamBegin{callback} an instance of
      `std::function<bool(double)>`. It is called regularly when the
      algorithm is running: the current advancement (between 0. and
      1.) is passed as parameter. If it returns `true`, then the
      algorithm continues its execution normally; if it returns
      `false`, the algorithm is stopped and the remaining normals are
      left unchanged.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
   \cgalNamedParamsEnd
*/
template <typename ConcurrencyTag,
	  typename PointRange,
          typename NamedParameters
>
void
jet_estimate_normals(
  PointRange& points,
  unsigned int k,
  const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  
  CGAL_TRACE("Calls jet_estimate_normals()\n");

  // basic geometric types
  typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::type NormalMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;
  typedef typename Kernel::FT FT;
  typedef typename GetSvdTraits<NamedParameters>::type SvdTraits;

  CGAL_static_assertion_msg(!(boost::is_same<NormalMap,
                              typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::NoMap>::value),
                            "Error: no normal map");
  CGAL_static_assertion_msg(!(boost::is_same<SvdTraits,
                              typename GetSvdTraits<NamedParameters>::NoTraits>::value),
                            "Error: no SVD traits");

  PointMap point_map = choose_parameter(get_parameter(np, internal_np::point_map), PointMap());
  NormalMap normal_map = choose_parameter(get_parameter(np, internal_np::normal_map), NormalMap());
  unsigned int degree_fitting = choose_parameter(get_parameter(np, internal_np::degree_fitting), 2);
  FT neighbor_radius = choose_parameter(get_parameter(np, internal_np::neighbor_radius), FT(0));
  
  const std::function<bool(double)>& callback = choose_parameter(get_parameter(np, internal_np::callback),
                                                               std::function<bool(double)>());

  typedef typename Kernel::Point_3 Point;

  // Input points types
  typedef typename boost::property_traits<NormalMap>::value_type Vector;

  // types for K nearest neighbors search structure
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(points.begin() != points.end());

  // precondition: at least 2 nearest neighbors
  CGAL_point_set_processing_precondition(k >= 2 || neighbor_radius > FT(0));

  std::size_t memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
  CGAL_TRACE("  Creates KD-tree\n");

  typename PointRange::iterator it;

  // Instanciate a KD-tree search.
  // Note: We have to convert each input iterator to Point_3.
  std::vector<Point> kd_tree_points; 
  for(it = points.begin(); it != points.end(); it++)
    kd_tree_points.push_back(get(point_map, *it));
  Tree tree(kd_tree_points.begin(), kd_tree_points.end());

  memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
  CGAL_TRACE("  Computes normals\n");

  // iterate over input points, compute and output normal
  // vectors (already normalized)
#ifndef CGAL_LINKED_WITH_TBB
  CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
			     "Parallel_tag is enabled but TBB is unavailable.");
#else
   if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
   {
     Point_set_processing_3::internal::Parallel_callback
       parallel_callback (callback, kd_tree_points.size());
     
     std::vector<Vector> normals (kd_tree_points.size (),
                                  CGAL::NULL_VECTOR);
     CGAL::internal::Jet_estimate_normals<Kernel, SvdTraits, Tree>
       f (tree, k, neighbor_radius,
          kd_tree_points, degree_fitting, normals,
          parallel_callback.advancement(),
          parallel_callback.interrupted());
     tbb::parallel_for(tbb::blocked_range<size_t>(0, kd_tree_points.size ()), f);
     std::size_t i = 0;
     for(it = points.begin(); it != points.end(); ++ it, ++ i)
       if (normals[i] != CGAL::NULL_VECTOR)
         put (normal_map, *it, normals[i]);

     parallel_callback.join();
   }
   else
#endif
     {
       std::size_t nb = 0;
       for(it = points.begin(); it != points.end(); it++, ++ nb)
	 {
	   Vector normal = internal::jet_estimate_normal<Kernel,SvdTraits,Tree>(
										get(point_map,*it), 
										tree, k, neighbor_radius, degree_fitting);

	   put(normal_map, *it, normal); // normal_map[it] = normal
           if (callback && !callback ((nb+1) / double(kd_tree_points.size())))
             break;
    	 }
     }


  memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
  CGAL_TRACE("End of jet_estimate_normals()\n");
}


/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename ConcurrencyTag,
	  typename PointRange>
void
jet_estimate_normals(
  PointRange& points,
  unsigned int k) ///< number of neighbors.
{
  jet_estimate_normals<ConcurrencyTag>
    (points, k, CGAL::Point_set_processing_3::parameters::all_default(points));
}
/// \endcond

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_JET_ESTIMATE_NORMALS_H
