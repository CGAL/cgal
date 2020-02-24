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

#ifndef CGAL_BILATERAL_SMOOTH_POINT_SET_H
#define CGAL_BILATERAL_SMOOTH_POINT_SET_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/number_type_config.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Point_set_processing_3/internal/neighbor_query.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/squared_distance_3.h>
#include <functional>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iterator>
#include <set>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <CGAL/Real_timer.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/property_map.h>

#ifdef CGAL_LINKED_WITH_TBB

#include <CGAL/Point_set_processing_3/internal/Parallel_callback.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>  
#include <atomic>
#endif // CGAL_LINKED_WITH_TBB

// Default allocator: use TBB allocators if available
#ifdef CGAL_LINKED_WITH_TBB
# define CGAL_PSP3_DEFAULT_ALLOCATOR tbb::scalable_allocator
#else // CGAL_LINKED_WITH_TBB
# define CGAL_PSP3_DEFAULT_ALLOCATOR std::allocator
#endif // CGAL_LINKED_WITH_TBB


//#define CGAL_PSP3_VERBOSE 

namespace CGAL {

// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL

namespace bilateral_smooth_point_set_internal{

// Item in the Kd-tree: position (Point_3) + index
template <typename Kernel>
class Kd_tree_element : public Point_with_normal_3<Kernel>
{
public:
  unsigned int index;

  // basic geometric types
  typedef typename CGAL::Origin Origin;
  typedef CGAL::Point_with_normal_3<Kernel> Base;

  Kd_tree_element(const Origin& o = ORIGIN, unsigned int id=0)
    : Base(o), index(id)
  {}
  Kd_tree_element(const Base& p, unsigned int id=0)
    : Base(p), index(id)
  {}
  Kd_tree_element(const Kd_tree_element& other)
    : Base(other), index(other.index)
  {}
  
  Kd_tree_element& operator=(const Kd_tree_element&)=default;
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


/// Compute bilateral projection for each point
/// according to their KNN neighborhood points
/// 
/// \pre `k >= 2`, radius > 0 , sharpness_angle > 0 && sharpness_angle < 90
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return 

template <typename Kernel>
CGAL::Point_with_normal_3<Kernel>
compute_denoise_projection(
  const CGAL::Point_with_normal_3<Kernel>& query, ///< 3D point to project
  const std::vector<CGAL::Point_with_normal_3<Kernel>,
  CGAL_PSP3_DEFAULT_ALLOCATOR<CGAL::Point_with_normal_3<Kernel> > >& neighbor_pwns,  //
  typename Kernel::FT radius,                   ///< accept neighborhood radius
  typename Kernel::FT sharpness_angle           ///< control sharpness(0-90)
)
{
  CGAL_point_set_processing_precondition(radius > 0);
  CGAL_point_set_processing_precondition(sharpness_angle > 0
                                         && sharpness_angle < 90);

  // basic geometric types
  typedef typename Kernel::FT FT;
  typedef CGAL::Point_with_normal_3<Kernel> Pwn;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Point_3 Point;

  FT radius2 = radius * radius;

  FT weight = (FT)0.0;
  FT iradius16 = -(FT)4.0/radius2;
  FT project_dist_sum = FT(0.0);
  FT project_weight_sum = FT(0.0);
  Vector normal_sum = CGAL::NULL_VECTOR; 

  FT cos_sigma = cos(sharpness_angle * CGAL_PI / 180.0);
  FT sharpness_bandwidth = std::pow((CGAL::max)(1e-8, 1 - cos_sigma), 2);

  typename std::vector<Pwn,CGAL_PSP3_DEFAULT_ALLOCATOR<Pwn> >::const_iterator 
    pwn_iter = neighbor_pwns.begin();
  for (; pwn_iter != neighbor_pwns.end(); ++pwn_iter)
  {
    const Point& np = pwn_iter->position();
    const Vector& nn = pwn_iter->normal();

    FT dist2 = CGAL::squared_distance(query.position(), np);
    if (dist2 < radius2)
    {
      FT theta = std::exp(dist2 * iradius16);
      FT psi = std::exp(-std::pow(1 - query.normal() * nn, 2)
        / sharpness_bandwidth);

      weight = theta * psi;

      project_dist_sum += ((query.position() - np) * nn) * weight;
      project_weight_sum += weight;
      normal_sum = normal_sum + nn * weight;
    }
  }

  Vector update_normal = normal_sum / project_weight_sum;
  update_normal = update_normal / sqrt(update_normal.squared_length());

  Point update_point = query.position() - update_normal * 
                      (project_dist_sum / project_weight_sum); 

  return Pwn(update_point, update_normal);
}

/// Computes max-spacing of one query point from K nearest neighbors.
///
/// \pre `k >= 2`.
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return max spacing.
template < typename Kernel,
           typename Tree >
typename Kernel::FT
compute_max_spacing(
  const CGAL::Point_with_normal_3<Kernel>& query, ///< 3D point
  Tree& tree,                                     ///< KD-tree
  unsigned int k)                                 ///< number of neighbors
{
  // basic geometric types
  typedef typename Kernel::FT FT;
  typedef CGAL::Point_with_normal_3<Kernel> Pwn;

  // types for K nearest neighbors search
  typedef bilateral_smooth_point_set_internal::Kd_tree_traits<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;

  // performs k + 1 queries (if unique the query point is
  // output first). search may be aborted when k is greater
  // than number of input points
  Neighbor_search search(tree,query,k+1);
  Search_iterator search_iterator = search.begin();
  ++search_iterator;
  FT max_distance = (FT)0.0;
  unsigned int i;
  for(i = 0; i < (k+1) ; ++i)
  {
    if(search_iterator == search.end())
      break; // premature ending

    Pwn pwn = search_iterator->first;
    double dist2 = CGAL::squared_distance(query.position(), pwn.position());
    max_distance = (CGAL::max)(dist2, max_distance);
    ++search_iterator;
  }

  // output max spacing
  return std::sqrt(max_distance);
}

} /* namespace internal */

/// \endcond

#ifdef CGAL_LINKED_WITH_TBB
/// \cond SKIP_IN_MANUAL
/// This is for parallelization of function: bilateral_smooth_point_set()
template <typename Kernel, typename Tree>
class Compute_pwns_neighbors
{
  typedef typename CGAL::Point_with_normal_3<Kernel> Pwn;
  typedef typename std::vector<Pwn,CGAL_PSP3_DEFAULT_ALLOCATOR<Pwn> > Pwns;
  typedef typename std::vector<Pwns,CGAL_PSP3_DEFAULT_ALLOCATOR<Pwns> > 
                                                                Pwns_neighbors;
  typedef typename Kernel::FT FT;

  unsigned int                                              m_k;
  FT                                                        m_neighbor_radius;
  const Tree                                              & m_tree;
  const Pwns                                              & m_pwns;
  Pwns_neighbors                                          & m_pwns_neighbors;
  cpp11::atomic<std::size_t>& advancement;
  cpp11::atomic<bool>& interrupted;

public:
  Compute_pwns_neighbors(unsigned int k, FT neighbor_radius, const Tree &tree,
                         const Pwns &pwns, Pwns_neighbors &neighbors,
                         cpp11::atomic<std::size_t>& advancement,
                         cpp11::atomic<bool>& interrupted)
    : m_k(k), m_neighbor_radius (neighbor_radius), m_tree(tree)
    , m_pwns(pwns), m_pwns_neighbors(neighbors)
    , advancement (advancement), interrupted (interrupted) {} 

  void operator() ( const tbb::blocked_range<size_t>& r ) const 
  {
    for (size_t i = r.begin(); i!=r.end(); i++)
    {
      if (interrupted)
        break;

      CGAL::Point_set_processing_3::internal::neighbor_query
        (m_pwns[i], m_tree, m_k, m_neighbor_radius, m_pwns_neighbors[i]);
      
      ++ advancement;
    }
  }
};
/// \endcond  

/// \cond SKIP_IN_MANUAL
/// This is for parallelization of function: compute_denoise_projection()
template <typename Kernel>
class Pwn_updater 
{
  typedef typename CGAL::Point_with_normal_3<Kernel> Pwn;
  typedef typename std::vector<Pwn,CGAL_PSP3_DEFAULT_ALLOCATOR<Pwn> > Pwns;
  typedef typename Kernel::FT FT;

  FT sharpness_angle;
  FT radius;
  Pwns* pwns;
  Pwns* update_pwns;
  std::vector<Pwns,CGAL_PSP3_DEFAULT_ALLOCATOR<Pwns> >* pwns_neighbors;
  cpp11::atomic<std::size_t>& advancement;
  cpp11::atomic<bool>& interrupted;

public:
  Pwn_updater(FT sharpness, 
    FT r,
    Pwns *in,
    Pwns *out, 
    std::vector<Pwns,CGAL_PSP3_DEFAULT_ALLOCATOR<Pwns> >* neighbors,
    cpp11::atomic<std::size_t>& advancement,
    cpp11::atomic<bool>& interrupted): 
  sharpness_angle(sharpness), 
    radius(r),
    pwns(in),
    update_pwns(out),
    pwns_neighbors(neighbors),
    advancement (advancement),
    interrupted (interrupted) {} 

  void operator() ( const tbb::blocked_range<size_t>& r ) const 
  { 
    for (size_t i = r.begin(); i != r.end(); ++i) 
    {
      if (interrupted)
        break;
      (*update_pwns)[i] = bilateral_smooth_point_set_internal::
        compute_denoise_projection<Kernel>((*pwns)[i], 
        (*pwns_neighbors)[i], 
        radius,
        sharpness_angle);  
      ++ advancement;
    }
  }
};
/// \endcond
#endif // CGAL_LINKED_WITH_TBB


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------


/**
   \ingroup PkgPointSetProcessing3Algorithms
 
   This function smooths an input point set by iteratively projecting each 
   point onto the implicit surface patch fitted over its nearest neighbors.
   Bilateral projection preserves sharp features according to the normal
   (gradient) information. Both point positions and normals will be modified.  
   For more details, please see section 4 in \cgalCite{ear-2013}.  

   A parallel version of this function is provided and requires the executable to be 
   linked against the <a href="https://www.threadingbuildingblocks.org">Intel TBB library</a>.
   To control the number of threads used, the user may use the tbb::task_scheduler_init class.
   See the <a href="https://www.threadingbuildingblocks.org/documentation">TBB documentation</a> 
   for more details.

   \pre Normals must be unit vectors
   \pre k >= 2

   \tparam ConcurrencyTag enables sequential versus parallel algorithm. Possible values are `Sequential_tag`,
                          `Parallel_tag`, and `Parallel_if_available_tag`.
   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param k size of the neighborhood for the implicit surface patch fitting.
   The larger the value is, the smoother the result will be.
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadWritePropertyMap` with value type `geom_traits::Point_3`.
     If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{normal_map} a model of `ReadWritePropertyMap` with value type
     `geom_traits::Vector_3`.\cgalParamEnd
     \cgalParamBegin{neighbor_radius} spherical neighborhood radius. If
     provided, the neighborhood of a query point is computed with a fixed spherical
     radius instead of a fixed number of neighbors. In that case, the parameter
     `k` is used as a limit on the number of points returned by each spherical
     query (to avoid overly large number of points in high density areas). If no
     limit is wanted, use `k=0`.\cgalParamEnd
     \cgalParamBegin{sharpness_angle} controls the sharpness of the result.\cgalParamEnd
     \cgalParamBegin{callback} an instance of
      `std::function<bool(double)>`. It is called regularly when the
      algorithm is running: the current advancement (between 0. and
      1.) is passed as parameter. If it returns `true`, then the
      algorithm continues its execution normally; if it returns
      `false`, the algorithm is stopped, all points are left unchanged
      and the function return `NaN`.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
   \cgalNamedParamsEnd

   \return Average point movement error. It's a convergence criterium for the algorithm.
   This value can help the user to decide how many iterations are
   sufficient.
*/
template <typename ConcurrencyTag,
          typename PointRange,
          typename NamedParameters>
double
bilateral_smooth_point_set(
  PointRange& points,
  unsigned int k,
  const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  
  // basic geometric types
  typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::type NormalMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;

  CGAL_static_assertion_msg(!(boost::is_same<NormalMap,
                              typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::NoMap>::value),
                            "Error: no normal map");
  
  typedef typename CGAL::Point_with_normal_3<Kernel> Pwn;
  typedef typename std::vector<Pwn,CGAL_PSP3_DEFAULT_ALLOCATOR<Pwn> > Pwns;
  typedef typename Kernel::FT FT;
  
  double sharpness_angle = choose_parameter(get_parameter(np, internal_np::sharpness_angle), 30.);
  const std::function<bool(double)>& callback = choose_parameter(get_parameter(np, internal_np::callback),
                                                                 std::function<bool(double)>());
  
  CGAL_point_set_processing_precondition(points.begin() != points.end());
  CGAL_point_set_processing_precondition(k > 1);

  // types for K nearest neighbors search structure
  typedef bilateral_smooth_point_set_internal::
                                       Kd_tree_element<Kernel> Kd_tree_element;
  typedef bilateral_smooth_point_set_internal::Kd_tree_traits<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;

  PointMap point_map = choose_parameter(get_parameter(np, internal_np::point_map), PointMap());
  NormalMap normal_map = choose_parameter(get_parameter(np, internal_np::normal_map), NormalMap());
  FT neighbor_radius = choose_parameter(get_parameter(np, internal_np::neighbor_radius), FT(0));
  
  // copy points and normals
  Pwns pwns;
  for(typename PointRange::iterator it = points.begin(); it != points.end(); ++it)
  {
    typename boost::property_traits<PointMap>::reference p = get(point_map, *it);
    typename boost::property_traits<NormalMap>::reference n = get(normal_map, *it);
    CGAL_point_set_processing_precondition(n.squared_length() > 1e-10);
    
    pwns.push_back(Pwn(p, n));
  }

  std::size_t nb_points = pwns.size();

#ifdef CGAL_PSP3_VERBOSE
   std::cout << "Initialization and compute max spacing: " << std::endl;
#endif
   // initiate a KD-tree search for points
   std::vector<Kd_tree_element,
     CGAL_PSP3_DEFAULT_ALLOCATOR<Kd_tree_element> > treeElements;
   treeElements.reserve(pwns.size());
   typename std::vector<Pwn,CGAL_PSP3_DEFAULT_ALLOCATOR<Pwn> >::iterator 
     pwn_iter = pwns.begin();
   for (unsigned int i = 0; pwn_iter != pwns.end(); ++pwn_iter)
   {
     treeElements.push_back(Kd_tree_element(*pwn_iter, i));
   }
   Tree tree(treeElements.begin(), treeElements.end());
   // Guess spacing
#ifdef CGAL_PSP3_VERBOSE
   CGAL::Real_timer task_timer;
   task_timer.start();
#endif
   FT guess_neighbor_radius = 0.0; 

   for(pwn_iter = pwns.begin(); pwn_iter != pwns.end(); ++pwn_iter)
   {
     FT max_spacing = bilateral_smooth_point_set_internal::
       compute_max_spacing<Kernel,Tree>(*pwn_iter, tree, k);
     guess_neighbor_radius = (CGAL::max)(max_spacing, guess_neighbor_radius); 
   }
   
#ifdef CGAL_PSP3_VERBOSE
   task_timer.stop();
#endif
   guess_neighbor_radius *= 0.95;

#ifdef CGAL_PSP3_VERBOSE
   CGAL::Memory_sizer::size_type memory = CGAL::Memory_sizer().virtual_size();
   std::cout << "done: " << task_timer.time() << " seconds, "
             << (memory>>20) << " Mb allocated" << std::endl;

   std::cout << "Compute all neighbors: " << std::endl;
   task_timer.reset();
   task_timer.start();
#endif
   // compute all neighbors
   std::vector<Pwns,CGAL_PSP3_DEFAULT_ALLOCATOR<Pwns> > pwns_neighbors;
   pwns_neighbors.resize(nb_points);
 
#ifndef CGAL_LINKED_WITH_TBB
  CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
			     "Parallel_tag is enabled but TBB is unavailable.");
#else
   if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
   {
     Point_set_processing_3::internal::Parallel_callback
       parallel_callback (callback, 2 * nb_points);

     Compute_pwns_neighbors<Kernel, Tree> f(k, neighbor_radius, tree, pwns, pwns_neighbors,
                                            parallel_callback.advancement(),
                                            parallel_callback.interrupted());
     tbb::parallel_for(tbb::blocked_range<size_t>(0, nb_points), f);

     bool interrupted = parallel_callback.interrupted();
  
     // We interrupt by hand as counter only goes halfway and won't terminate by itself
     parallel_callback.interrupted() = true;
     parallel_callback.join();       

     // If interrupted during this step, nothing is computed, we return NaN
     if (interrupted)
       return std::numeric_limits<double>::quiet_NaN();
   }
   else
#endif
   {
     typename std::vector<Pwns,CGAL_PSP3_DEFAULT_ALLOCATOR<Pwns> >::iterator 
       pwns_iter = pwns_neighbors.begin();

     std::size_t nb = 0;
     for(pwn_iter = pwns.begin(); pwn_iter != pwns.end(); ++pwn_iter, ++pwns_iter, ++ nb)
     {
       CGAL::Point_set_processing_3::internal::neighbor_query
         (*pwn_iter, tree, k, neighbor_radius, *pwns_iter);

       if (callback && !callback ((nb+1) / double(2. * nb_points)))
         return std::numeric_limits<double>::quiet_NaN();
     }
   }
   
#ifdef CGAL_PSP3_VERBOSE
   task_timer.stop();
   memory = CGAL::Memory_sizer().virtual_size();
   std::cout << "done: " << task_timer.time() << " seconds, "
             << (memory>>20) << " Mb allocated" << std::endl;

   std::cout << "Compute update points and normals: " << std::endl;
   task_timer.reset();
   task_timer.start();
#endif
   // update points and normals
   Pwns update_pwns(nb_points);

#ifdef CGAL_LINKED_WITH_TBB
   if(boost::is_convertible<ConcurrencyTag, CGAL::Parallel_tag>::value)
   {
     Point_set_processing_3::internal::Parallel_callback
       parallel_callback (callback, 2 * nb_points, nb_points);
     
     //tbb::task_scheduler_init init(4);
     tbb::blocked_range<size_t> block(0, nb_points);
     Pwn_updater<Kernel> pwn_updater(sharpness_angle,
                                     guess_neighbor_radius,
                                     &pwns,
                                     &update_pwns,
                                     &pwns_neighbors,
                                     parallel_callback.advancement(),
                                     parallel_callback.interrupted());
     tbb::parallel_for(block, pwn_updater);

     parallel_callback.join();

     // If interrupted during this step, nothing is computed, we return NaN
     if (parallel_callback.interrupted())
       return std::numeric_limits<double>::quiet_NaN();
   }
   else
#endif // CGAL_LINKED_WITH_TBB
   {
     std::size_t nb = nb_points;
     
     typename std::vector<Pwn,CGAL_PSP3_DEFAULT_ALLOCATOR<Pwn> >::iterator 
       update_iter = update_pwns.begin();
     typename std::vector<Pwns,CGAL_PSP3_DEFAULT_ALLOCATOR<Pwns> >::iterator 
       neighbor_iter = pwns_neighbors.begin();
     for(pwn_iter = pwns.begin(); pwn_iter != pwns.end(); 
         ++pwn_iter, ++update_iter, ++neighbor_iter, ++ nb)
     {
       *update_iter = bilateral_smooth_point_set_internal::
         compute_denoise_projection<Kernel>
         (*pwn_iter, 
          *neighbor_iter, 
          guess_neighbor_radius, 
          sharpness_angle);
       if (callback && !callback ((nb+1) / double(2. * nb_points)))
         return std::numeric_limits<double>::quiet_NaN();
     }
   }
#ifdef CGAL_PSP3_VERBOSE
   task_timer.stop(); 
   memory = CGAL::Memory_sizer().virtual_size();
   std::cout << "done: " << task_timer.time() << " seconds, "
             << (memory>>20) << " Mb allocated" << std::endl;
#endif
   // save results
   FT sum_move_error = 0;
   typename PointRange::iterator it = points.begin();
   for(unsigned int i = 0 ; it != points.end(); ++it, ++i)
   {
     typename boost::property_traits<PointMap>::reference p = get(point_map, *it);
     sum_move_error += CGAL::squared_distance(p, update_pwns[i].position());
     put (point_map, *it, update_pwns[i].position());
     put (normal_map, *it, update_pwns[i].normal());
   }
     
   return sum_move_error / nb_points;
}

/// \cond SKIP_IN_MANUAL
// variant with default NP  
template <typename ConcurrencyTag,
          typename PointRange>
double
bilateral_smooth_point_set(
  PointRange& points,
  unsigned int k)           ///< size of the neighborhood for the implicit surface patch fitting.
                            ///< The larger the value is, the smoother the result will be.
{
  return bilateral_smooth_point_set<ConcurrencyTag>
    (points, k, CGAL::Point_set_processing_3::parameters::all_default(points));
}
/// \endcond


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_BILATERAL_SMOOTH_POINT_SET_H
