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
//
// Author(s) : Shihao Wu, Cl¨¦ment Jamin 

#ifndef CGAL_DENOSISE_POINTS_WITH_NORMALS_H
#define CGAL_DENOSISE_POINTS_WITH_NORMALS_H

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>

#include <iterator>
#include <set>
#include <algorithm>
#include <cmath>
#include <ctime>

//#include "tbb/parallel_for.h"
//#include "tbb/blocked_range.h"

// Not sure...
//struct Sequential_tag {};
//struct Parallel_tag : public Sequential_tag {};
//#ifdef CGAL_LINKED_WITH_TBB
//#define CGAL_LINKED_WITH_TBB

/// \cond SKIP_IN_MANUAL

class Timer
{
public:

  void start(const std::string& str)
  {
    std::cout << std::endl;
    starttime = clock();
    mid_start = clock();
   // std::cout << "@@@@@ Time Count Strat For: " << str << std::endl;

    _str = str;
  }

  void end()
  {
    stoptime = clock();
    timeused = stoptime - starttime;
    std::cout << /*endl <<*/ "@@@@ finish	" << _str << "  time used:  " << timeused / double(CLOCKS_PER_SEC) << " seconds." << std::endl;
    std::cout << std::endl;
  }

private:
  int starttime, mid_start, mid_end, stoptime, timeused;
  std::string _str;
};


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace denoise_points_with_normals_internal{

  // Item in the Kd-tree: position (Point_3) + index
  template <typename Kernel>
  class KdTreeElement : public Kernel::Point_3
  {
  public:
    unsigned int index;

    // basic geometric types
    typedef typename CGAL::Origin Origin;
    typedef typename Kernel::Point_3 Point;

    KdTreeElement(const Origin& o = ORIGIN, unsigned int id=0)
      : Point(o), index(id)
    {}
    KdTreeElement(const Point& p, unsigned int id=0)
      : Point(p), index(id)
    {}
    KdTreeElement(const KdTreeElement& other)
      : Point(other), index(other.index)
    {}
  };

  // Helper class for the Kd-tree
  template <typename Kernel>
  class KdTreeGT : public Kernel
  {
  public:
    typedef KdTreeElement<Kernel> Point_3;
  };

  template <typename Kernel>
  class KdTreeTraits : public CGAL::Search_traits_3<KdTreeGT<Kernel> >
  {
  public:
    typedef typename Kernel::Point_3 PointType;
  };

  /// compute rimls projection for each point
  /// according to their KNN neighborhood sample points
  /// 
  /// \pre `k >= 2`, radius > 0
  ///
  /// @tparam Kernel Geometric traits class.
  /// @tparam Tree KD-tree.
  ///
  /// @return 
  template <typename Kernel,
    typename Tree>
    void
    compute_denoise_projection(
    const typename Kernel::Point_3& query, ///< 3D point to project
    const typename Kernel::Vector_3& query_nromal, ///< normal of query point
    Tree& tree, ///< KD-tree
    const unsigned int k, ///< nb neighbors
    const typename Kernel::FT radius, ///< accept neighborhood radius
    const std::vector<typename Kernel::Vector_3>& normal_set, ///< normal set
    typename Kernel::Point_3&  update_point, ///< return point
    typename Kernel::Vector_3& update_normal ///< return normal
    )
  {
    CGAL_point_set_processing_precondition( k > 1);
    CGAL_point_set_processing_precondition(radius > 0);

    // basic geometric types
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::FT FT;

    FT radius2 = radius * radius;

    // types for K nearest neighbors search
    typedef denoise_points_with_normals_internal::KdTreeTraits<Kernel> Tree_traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
    typedef typename Neighbor_search::iterator Search_iterator;


    // Gather set of k neighboring original points.
    std::vector<Point> neighbor_points; 
    std::vector<Vector> neighbor_normals; 

    neighbor_points.reserve(k);
    neighbor_normals.reserve(k);
    Neighbor_search search(tree, query, k+1);
    Search_iterator search_iterator = search.begin();
    ++search_iterator;
    std::vector<FT> dist2_set;
    for(unsigned int i = 0; i < k; i++)
    {
      if(search_iterator == search.end())
        break; // premature ending

      Point& np = search_iterator->first;
      FT dist2 = CGAL::squared_distance(query, np);
      if (dist2 < radius2)
      {
        neighbor_points.push_back(search_iterator->first);
        neighbor_normals.push_back(normal_set[search_iterator->first.index]);
        dist2_set.push_back(dist2);
      }

      ++search_iterator;
    }

    if (neighbor_points.empty())
    {
      update_point = query;
      update_normal = query_nromal;
      return;
    }


    //Compute 
    FT weight = (FT)0.0;
    FT iradius16 = -(FT)4.0/radius2;
    FT project_dist_sum = FT(0.0);
    FT project_weight_sum = FT(0.0);
    Vector normal_sum = CGAL::NULL_VECTOR; 

    FT sigma = 45; // should be a parameter
    //FT sharpness_bandwidth = std::pow(std::max(1e-8,1-cos(sigma/180.0*3.1415926)), 2);
    FT sharpness_bandwidth = std::pow(1-cos(sigma/180.0*3.1415926), 2);
    
    for (unsigned int i = 0; i < neighbor_points.size(); i++)
    {
      Point& np = neighbor_points[i];
      Vector& nn = neighbor_normals[i];

      FT dist2 = dist2_set[i];
      FT theta = std::exp(dist2 * iradius16);
      FT psi = std::exp(-std::pow(1 - query_nromal * nn, 2) / sharpness_bandwidth);
     
      weight = theta * psi;

      project_dist_sum += ((query - np) * nn) * weight;
      project_weight_sum += weight;
      normal_sum = normal_sum + nn;
    }

    // output
    update_normal = normal_sum / project_weight_sum;
    update_normal = update_normal / sqrt(update_normal.squared_length());
    update_point = query - update_normal * (project_dist_sum / project_weight_sum); 
  }


  /// Computes max-spacing of one query point from K nearest neighbors.
  ///
  /// \pre `k >= 2`.
  ///
  /// @tparam Kernel Geometric traits class.
  /// @tparam Tree KD-tree.
  ///
  /// @return average spacing (scalar).
  template < typename Kernel,
    typename Tree >
    typename Kernel::FT
    compute_max_spacing(const typename Kernel::Point_3& query, ///< 3D point whose spacing we want to compute
    Tree& tree,                            ///< KD-tree
    unsigned int k)                        ///< number of neighbors
  {
    // basic geometric types
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_3 Point;

    // types for K nearest neighbors search
    typedef denoise_points_with_normals_internal::KdTreeTraits<Kernel> Tree_traits;
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
    for(i=0;i<(k+1);i++)
    {
      if(search_iterator == search.end())
        break; // premature ending

      Point p = search_iterator->first;
      double dist2 = CGAL::squared_distance(query,p);
      max_distance = dist2 > max_distance ? dist2 : max_distance;// can be simplify, no need to compare..
      ++search_iterator;
    }

    // output average max spacing
    return std::sqrt(max_distance);
  }
}

// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------
namespace CGAL {

//===================================================================================
/// \ingroup PkgPointSetProcessing
/// 
/// A function for bilateral point set denoising (smoothing) with sharp features
/// \pre normals must be unit vectors
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = Point_3<Kernel>.
///        It can be omitted if ForwardIterator value_type is convertible to Point_3<Kernel>.
/// @tparam NormalPMap is a model of `WritablePropertyMap` with a value_type = Vector_3<Kernel>.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from PointPMap value_type.
///
/// @return average point move error.

// This variant requires all parameters.
template <typename ForwardIterator,
  typename PointPMap,
  typename NormalPMap,
  typename Kernel
>
double
denoise_points_with_normals(
ForwardIterator first,  ///< iterator over the first input point.
ForwardIterator beyond, ///< past-the-end iterator over the input points.
PointPMap point_pmap, ///< property map ForwardIterator -> Point_3.
NormalPMap normal_pmap, ///< property map ForwardIterator -> Vector_3.
const unsigned int k, ///< number of neighbors.
const Kernel& /*kernel*/) ///< geometric traits.
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;

  CGAL_point_set_processing_precondition(first != beyond);
  CGAL_point_set_processing_precondition(k > 1);

  // types for K nearest neighbors search structure
  typedef denoise_points_with_normals_internal::KdTreeElement<Kernel> KdTreeElement;
  typedef denoise_points_with_normals_internal::KdTreeTraits<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Search_iterator;

  // copy points and normals
  std::vector<Point> point_set;
  std::vector<Vector> normal_set;
  for(ForwardIterator it = first; it != beyond; ++it)
  {
    point_set.push_back(get(point_pmap, it));
    normal_set.push_back(get(normal_pmap, it));
  }

  // initiate a KD-tree search for points
  unsigned int i;
  std::vector<KdTreeElement> treeElements;
  for (i = 0 ; i < point_set.size(); i++)
  {
    Point& p0 = point_set[i];
    treeElements.push_back(KdTreeElement(p0,i));
  }
  Tree tree(treeElements.begin(), treeElements.end());

  // Guess spacing
  FT guess_neighbor_radius = (FT)(std::numeric_limits<double>::max)(); // Or a better max number: (numeric_limits<double>::max)()?
  for(ForwardIterator it = first; it != beyond ; ++it)
  {
    FT max_spacing = denoise_points_with_normals_internal::compute_max_spacing<Kernel,Tree>(get(point_pmap,it),tree, k);
    guess_neighbor_radius = max_spacing < guess_neighbor_radius ? max_spacing : guess_neighbor_radius;
  }
  guess_neighbor_radius *= 0.95;

  
  std::cout << "Guess Neighborhood Radius:" << guess_neighbor_radius << std::endl;

  // update points and normals
  std::vector<Point> update_point_set(point_set.size());
  std::vector<Vector> update_normal_set(point_set.size());
  for (i = 0 ; i < point_set.size(); i++)
  {
    Point& p0 = point_set[i];
    Vector& n0 = normal_set[i];

    denoise_points_with_normals_internal::compute_denoise_projection<Kernel, Tree>(
      p0, 
      n0,
      tree,
      k,
      guess_neighbor_radius,
      normal_set,
      update_point_set[i],
      update_normal_set[i]);
  }

  for (i = 0 ; i < point_set.size(); i++)
  {
    Point& p0 = point_set[i];
    Vector& n0 = normal_set[i];
    
    p0 = update_point_set[i];
    n0 = update_normal_set[i];
  }

  // save results
  FT sum_move_error = 0;
  ForwardIterator it;
  for(i = 0, it = first; it != beyond; ++it, i++)
  {
    Point& p = get(point_pmap, it);
    Vector& n = get(normal_pmap, it);

    sum_move_error += CGAL::squared_distance(p, point_set[i]);
    p = point_set[i];
    n = normal_set[i];
  }

  return sum_move_error / point_set.size();
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the point property map.
template <typename ForwardIterator,
  typename PointPMap,
  typename NormalPMap
>
double
denoise_points_with_normals(
ForwardIterator first, ///< first input point.
ForwardIterator beyond, ///< past-the-end input point.
PointPMap point_pmap, ///< property map OutputIterator -> Point_3.
NormalPMap normal_pmap,
const unsigned int k ///< number of neighbors.
) ///< property map OutputIterator -> Vector_3.
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return denoise_points_with_normals(
    first, beyond,
    point_pmap,
    normal_pmap,
    k,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Dereference_property_map.
template <typename ForwardIterator,
  typename NormalPMap
>
double
denoise_points_with_normals(
ForwardIterator first, ///< first input point.
ForwardIterator beyond, ///< past-the-end input point.
const unsigned int k, ///< number of neighbors.
NormalPMap normal_pmap) ///< property map OutputIterator -> Vector_3.
{
  return denoise_points_with_normals(
    first, beyond,
    make_dereference_property_map(first),
    normal_pmap, k);
}
      /// @endcond


} //namespace CGAL

#endif // CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_H
