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
// Author(s) : Shihao Wu, Clement Jamin 

#ifndef CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_H
#define CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_H

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
namespace regularize_and_simplify_internal{

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
    typedef regularize_and_simplify_internal::KdTreeTraits<Kernel> Tree_traits;
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


  /// Compute average term for each sample points
  /// According to their KNN neighborhood original points
  /// 
  /// \pre `k >= 2`, radius > 0
  ///
  /// @tparam Kernel Geometric traits class.
  /// @tparam Tree KD-tree.
  ///
  /// @return computed point
  template <typename Kernel,
    typename Tree>
    typename Kernel::Vector_3
    compute_average_term(
    const typename Kernel::Point_3& query, ///< 3D point to project
    Tree& tree, ///< KD-tree
    unsigned int& k, // nb neighbors
    const typename Kernel::FT radius, //accept neighborhood radius
    const std::vector<typename Kernel::FT>& density_weight_set //if user need density
    //const std::vector<unsigned int>& guess_knn_set//guess knn to spp
    )
  {
    CGAL_point_set_processing_precondition( k > 1);
    CGAL_point_set_processing_precondition(radius > 0);
    bool is_density_weight_set_empty = density_weight_set.empty();

    // basic geometric types
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::FT FT;

    FT radius2 = radius * radius;

    // types for K nearest neighbors search
    typedef regularize_and_simplify_internal::KdTreeTraits<Kernel> Tree_traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
    typedef typename Neighbor_search::iterator Search_iterator;


    // Gather set of k neighboring original points.
    std::vector<Point> neighbor_original_points; 
    neighbor_original_points.reserve(k);
    Neighbor_search search(tree, query, k);
    Search_iterator search_iterator = search.begin();
    std::vector<FT> dist2_set;
    std::vector<FT> density_set;
    unsigned int redundant_neighbor_number = 0;
    for(unsigned int i = 0; i < k; i++)
    {
      if(search_iterator == search.end())
        break; // premature ending


      Point& np = search_iterator->first;
      FT dist2 = CGAL::squared_distance(query, np);
      if (dist2 < radius2)
      {
        if (!is_density_weight_set_empty)
        {
          density_set.push_back(density_weight_set[search_iterator->first.index]);
        }
        dist2_set.push_back(dist2);
        neighbor_original_points.push_back(search_iterator->first);
      }
      else
      {
        redundant_neighbor_number++;
      }

      ++search_iterator;
    }
    //std::cout << k << std::endl;
    if (redundant_neighbor_number <= 0)
    {
      k *= 1.1;
    }

    if (neighbor_original_points.empty())
    {
      return query - CGAL::ORIGIN;
    }
    CGAL_point_set_processing_precondition(neighbor_original_points.size() >= 1);

    //Compute average term
    FT weight = (FT)0.0, average_weight_sum = (FT)0.0;
    FT iradius16 = -(FT)4.0/radius2;
    Vector average = CGAL::NULL_VECTOR; 
    for (unsigned int i = 0; i < neighbor_original_points.size(); i++)
    {

      Point& np = neighbor_original_points[i];

      FT dist2 = dist2_set[i];
      weight = exp(dist2 * iradius16);

      if(!is_density_weight_set_empty)
      {
        weight *= density_set[i];
      }

      average_weight_sum += weight;
      average = average + (np - CGAL::ORIGIN) * weight;
    }

    // output
    return average/average_weight_sum;
  }

  /// Compute repulsion term for each sample points
  /// According to their KNN neighborhood sample points
  /// 
  /// \pre `k >= 2`, radius > 0
  ///
  /// @tparam Kernel Geometric traits class.
  /// @tparam Tree KD-tree.
  ///
  /// @return computed point
  template <typename Kernel,
    typename Tree>
    typename Kernel::Vector_3
    compute_repulsion_term(
    const typename Kernel::Point_3& query, ///< 3D point to project
    Tree& tree, ///< KD-tree
    const unsigned int k, // nb neighbors
    const typename Kernel::FT radius, //accept neighborhood radius
    const std::vector<typename Kernel::FT>& density_weight_set //if user need density
    )
  {
    CGAL_point_set_processing_precondition( k > 1);
    CGAL_point_set_processing_precondition(radius > 0);
    bool is_density_weight_set_empty = density_weight_set.empty();

    // basic geometric types
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::FT FT;

    FT radius2 = radius * radius;

    // types for K nearest neighbors search
    typedef regularize_and_simplify_internal::KdTreeTraits<Kernel> Tree_traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
    typedef typename Neighbor_search::iterator Search_iterator;


    // Gather set of k neighboring original points.
    std::vector<Point> neighbor_sample_points; 
    neighbor_sample_points.reserve(k);
    Neighbor_search search(tree, query, k+1);
    Search_iterator search_iterator = search.begin();
    ++search_iterator;
    std::vector<FT> dist2_set;
    std::vector<FT> density_set;
    for(unsigned int i = 0; i < k; i++)
    {
      if(search_iterator == search.end())
        break; // premature ending

      Point& np = search_iterator->first;
      FT dist2 = CGAL::squared_distance(query, np);
      if (dist2 < radius2)
      {
        if (!is_density_weight_set_empty)
        {
          density_set.push_back(density_weight_set[search_iterator->first.index]);
        }
        neighbor_sample_points.push_back(search_iterator->first);
        dist2_set.push_back(dist2);
      }

      ++search_iterator;
    }

    if (neighbor_sample_points.empty())
    {
      return CGAL::NULL_VECTOR; 
    }
    CGAL_point_set_processing_precondition(neighbor_sample_points.size() >= 1);

    //Compute average term
    FT weight = (FT)0.0, repulsion_weight_sum = (FT)0.0;
    FT iradius16 = -(FT)4.0/radius2;

    Vector repulsion = CGAL::NULL_VECTOR; 
    for (unsigned int i = 0; i < neighbor_sample_points.size(); i++)
    {
      Point& np = neighbor_sample_points[i];
      Vector diff = query - np;

      FT dist2 = dist2_set[i];
      FT dist = std::sqrt(dist2);

      weight = std::exp(dist2 * iradius16) * std::pow(FT(1.0)/dist, 2);
      if(!is_density_weight_set_empty)
      {
        weight *= density_set[i];
      }

      repulsion_weight_sum += weight;
      repulsion = repulsion + diff * weight;
    }

    // output
    return repulsion/repulsion_weight_sum;
  }




  /// Compute density weight for each original points,
  /// according to their KNN neighborhood original points
  /// 
  /// \pre `k >= 2`, radius > 0
  ///
  /// @tparam Kernel Geometric traits class.
  /// @tparam Tree KD-tree.
  ///
  /// @return computed point
  template <typename Kernel,
    typename Tree>
    typename Kernel::FT
    compute_density_weight_for_original_point(
    const typename Kernel::Point_3& query, ///< 3D point to project
    Tree& tree, ///< KD-tree
    const unsigned int k, // nb neighbors
    const typename Kernel::FT radius
    )
  {
    CGAL_point_set_processing_precondition( k > 1);
    CGAL_point_set_processing_precondition(radius > 0);

    // basic geometric types
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::FT FT;

    // types for K nearest neighbors search
    typedef regularize_and_simplify_internal::KdTreeTraits<Kernel> Tree_traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
    typedef typename Neighbor_search::iterator Search_iterator;


    //Compute density weight
    FT radius2 = radius * radius;
    FT density_weight = (FT)1.0, weight;
    FT iradius16 = -(FT)4.0/radius2;

    Neighbor_search search(tree, query, k);
    Search_iterator search_iterator = search.begin();
    for(unsigned int i = 0; i < k; i++)
    {
      if(search_iterator == search.end())
        break; // premature ending

      Point& np = search_iterator->first;
      FT dist2 = CGAL::squared_distance(query, np);
      if (dist2 < radius2)
      {
        weight = std::exp(dist2 * iradius16);
        density_weight += weight;
      }

      ++search_iterator;
    }

    // output
    return FT(1.0) / density_weight;
  }


  /// Compute density weight for sample point,
  /// according to their KNN neighborhood sample points
  /// 
  /// \pre `k >= 2`, radius > 0
  ///
  /// @tparam Kernel Geometric traits class.
  /// @tparam Tree KD-tree.
  ///
  /// @return computed point
  template <typename Kernel,
    typename Tree>
    typename Kernel::FT
    compute_density_weight_for_sample_point(
    const typename Kernel::Point_3& query, ///< 3D point to project
    Tree& tree, ///< KD-tree
    const unsigned int k, // nb neighbors
    const typename Kernel::FT radius
    )
  {
    // basic geometric types
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::FT FT;

    // types for K nearest neighbors search
    typedef regularize_and_simplify_internal::KdTreeTraits<Kernel> Tree_traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
    typedef typename Neighbor_search::iterator Search_iterator;

    //Compute density weight
    FT radius2 = radius * radius;
    FT density_weight = (FT)1.0;
    FT iradius16 = -(FT)4.0/radius2;

    Neighbor_search search(tree, query, k+1);
    Search_iterator search_iterator = search.begin();
    ++search_iterator;
    for(unsigned int i = 0; i < k; i++)
    {
      if(search_iterator == search.end())
        break; // premature ending

      Point& np = search_iterator->first;
      FT dist2 = CGAL::squared_distance(query, np);
      if (dist2 < radius2)
      {
        density_weight += std::exp(dist2 * iradius16);
      }

      ++search_iterator;
    }

    // output
    //return std::sqrt(density_weight);
    return density_weight;

  }

  /// Compute guess KNN number to speed up
  /// 
  /// \pre `k >= 2`, radius > 0
  ///
  /// @tparam Kernel Geometric traits class.
  /// @tparam Tree KD-tree.
  ///
  /// @return computed point
  template <typename Kernel,
    typename Tree>
    unsigned int
    guess_KNN_number_for_original_set(
    const typename Kernel::Point_3& query, ///< 3D point to project
    Tree& tree, ///< KD-tree
    const unsigned int k, // nb neighbors
    const typename Kernel::FT radius
    )
  {
    // basic geometric types
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::FT FT;


    // types for K nearest neighbors search
    typedef regularize_and_simplify_internal::KdTreeTraits<Kernel> Tree_traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
    typedef typename Neighbor_search::iterator Search_iterator;

    //Compute guess knn
    FT radius2 = radius * radius;
    Neighbor_search search(tree, query, k);
    Search_iterator search_iterator = search.begin();
    unsigned int guess_knn = 0;

    for(unsigned int i = 0; i < k; i++)
    {
      if(search_iterator == search.end())
        break; // premature ending

      Point& np = search_iterator->first;
      FT dist2 = CGAL::squared_distance(query, np);
      if (dist2 < radius2)
      {
        guess_knn++;
      }

      ++search_iterator;
    }

    // output
    return unsigned int(guess_knn * 1.2);

  }
}

// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------
namespace CGAL {

  /// \ingroup PkgPointSetProcessing
  /// WLOP Algorithm: The simplification algorithm can produces a set of 
  ///	denoised, outlier-free and evenly distributed particles over the original dense point cloud,
  ///	so as to improve the reliability of current normal orientation algorithm. 
  ///	The core of the algorithm is a Weighted Locally Optimal projection operator with a density uniformization term. 
  /// More deatail please see: http://web.siat.ac.cn/~huihuang/WLOP/WLOP_page.html
  ///
  /// @tparam ForwardIterator iterator over input points.
  /// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = Point_3<Kernel>.
  ///        It can be omitted if ForwardIterator value_type is convertible to Point_3<Kernel>.
  /// @tparam Kernel Geometric traits class.
  ///        It can be omitted and deduced automatically from PointPMap value_type.
  ///
  /// @return iterator of the first point to downsampled points.

  // This variant requires all parameters.
  template <typename ForwardIterator,
    typename PointPMap,
    typename Kernel
  >
  ForwardIterator
  regularize_and_simplify_point_set(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  double retain_percentage, ///< percentage of points to retain.
  unsigned int k, ///< number of neighbors.
  const unsigned int iter_number,///< number of iterations.
  const bool need_compute_density, ///< if needed to compute density to generate more rugularized result, 
  ///  especially when the density of input is uneven.
  const Kernel& /*kernel*/) ///< geometric traits.
  {
    CGAL_point_set_processing_precondition(k > 1);
    Timer time;

    // basic geometric types
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::FT FT;

    // types for K nearest neighbors search structure
    typedef regularize_and_simplify_internal::KdTreeElement<Kernel> KdTreeElement;
    typedef regularize_and_simplify_internal::KdTreeTraits<Kernel> Tree_traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
    typedef typename Neighbor_search::Tree Tree;
    typedef typename Neighbor_search::iterator Search_iterator;

    // precondition: at least one element in the container.
    // to fix: should have at least three distinct points
    // but this is costly to check
    CGAL_point_set_processing_precondition(first != beyond);
    CGAL_point_set_processing_precondition(retain_percentage >= 0 && retain_percentage <= 100);

    // Random shuffle
    std::random_shuffle (first, beyond);

    // Computes original(input) and sample points size 
    std::size_t nb_points_original = std::distance(first, beyond);
    std::size_t nb_points_sample = (std::size_t)(FT(nb_points_original) * (retain_percentage/100.0));
    std::size_t first_index_to_sample = nb_points_original - nb_points_sample;

    // The first point iter of original and sample points
    ForwardIterator it;// point iterator
    ForwardIterator first_original_point = first;
    ForwardIterator first_sample_point = first;
    std::advance(first_sample_point, first_index_to_sample);

    //Copy sample points
    std::vector<Point> sample_points(nb_points_sample);
    unsigned int i; // sample point index
    for(it = first_sample_point, i = 0; it != beyond; ++it, i++)
      sample_points[i] = get(point_pmap, it);

    // Initiate a KD-tree search for original points
    time.start("Build Original Neighbor Tree");
    std::vector<KdTreeElement> original_treeElements;
    for (it = first_original_point, i=0 ; it != beyond ; ++it, ++i)
    {
      Point& p0 = get(point_pmap,it);
      original_treeElements.push_back(KdTreeElement(p0,i));
    }
    Tree original_tree(original_treeElements.begin(), original_treeElements.end());
    time.end();

    // Guess spacing
    time.start("Guess Neighborhood Radiuse");
    FT guess_neighbor_radius = (FT)(std::numeric_limits<double>::max)(); // Or a better max number: (numeric_limits<double>::max)()?
    for(it = first_original_point; it != beyond ; ++it)
    {
      FT max_spacing = regularize_and_simplify_internal::compute_max_spacing<Kernel,Tree>(get(point_pmap,it),original_tree,k);
      guess_neighbor_radius = max_spacing < guess_neighbor_radius ? max_spacing : guess_neighbor_radius;
    }
    guess_neighbor_radius *= 0.95;
    time.end();
    std::cout << "Guess Neighborhood Radius:" << guess_neighbor_radius << std::endl;

    // Compute original density weight for original points if user needed
    time.start("Compute Density For Original");
    std::vector<FT> original_density_weight_set;
    if (need_compute_density)
    {
      for (it = first_original_point; it != beyond ; ++it)
      {
        FT density = regularize_and_simplify_internal::compute_density_weight_for_original_point<Kernel, Tree>(get(point_pmap,it), original_tree, k, guess_neighbor_radius * 0.3);
        original_density_weight_set.push_back(density);
      }
    }
    time.end();

    // Compute guess KNN set
    time.start("Compute guess KNN set");
    std::vector<unsigned int> guess_KNN_set;
    for (i=0 ; i < sample_points.size(); i++)
    {
      Point& p0 = sample_points[i];
      unsigned int guess_knn = regularize_and_simplify_internal::guess_KNN_number_for_original_set<Kernel, Tree>(p0, original_tree, k, guess_neighbor_radius);
      guess_KNN_set.push_back(guess_knn);
    }
    time.end();


    for (unsigned int iter_n = 0; iter_n < iter_number; iter_n++)
    {
      // Initiate a KD-tree search for sample points
      time.start("Build Sample Neighbor Tree");
      std::vector<KdTreeElement> sample_treeElements;
      unsigned int k_for_sample = 30; // Or it can be conducted by the "guess_neighbor_radius"

      for (i=0 ; i < sample_points.size(); i++)
      {
        Point& p0 = sample_points[i];
        sample_treeElements.push_back(KdTreeElement(p0,i));
      }
      Tree sample_tree(sample_treeElements.begin(), sample_treeElements.end());
      time.end();

      // Compute sample density weight for sample points if user needed
      std::vector<FT> sample_density_weight_set;
      time.start("Compute Density For Sample");
      if (need_compute_density)
      {
        for (i=0 ; i < sample_points.size(); i++)
        {
          FT density = regularize_and_simplify_internal::compute_density_weight_for_sample_point<Kernel, Tree>(sample_points[i], sample_tree, k_for_sample, guess_neighbor_radius);
          sample_density_weight_set.push_back(density);
        }
      }
      time.end();

      // Compute average term and repulsion term for each sample points separately,
      // then update each sample points
      std::vector<Vector> average_set(nb_points_sample);
      std::vector<Vector> repulsion_set(nb_points_sample);
      time.start("Compute Average Term");
      for (i = 0; i < sample_points.size(); i++)
      {
        Point& p = sample_points[i];
        average_set[i] = regularize_and_simplify_internal::compute_average_term<Kernel>(p, original_tree, k, guess_neighbor_radius, original_density_weight_set); // Before speed up
        //average_set[i] = regularize_and_simplify_internal::compute_average_term<Kernel>(p, original_tree, guess_KNN_set[i], guess_neighbor_radius, original_density_weight_set);

      }
      time.end();

      time.start("Compute Repulsion Term");
      for (i = 0; i < sample_points.size(); i++)
      {
        Point& p = sample_points[i];
        repulsion_set[i] = regularize_and_simplify_internal::compute_repulsion_term<Kernel>(p, sample_tree, k_for_sample, guess_neighbor_radius, sample_density_weight_set);
      }
      time.end();

      for (i = 0; i < sample_points.size(); i++)
      {
        Point& p = sample_points[i];
        p = CGAL::ORIGIN + average_set[i] + (FT)0.5 * repulsion_set[i];
      }

      std::cout << "iterate:	" << iter_n + 1 <<  "	"<< std::endl;
    }

    //Copy back modified sample points to original points for output
    for(it = first_sample_point, i = 0; it != beyond; ++it, i++)
    {
      Point& original_p = get(point_pmap, it);
      const Point& sample_p = sample_points[i];
      original_p = sample_p;
    }

    return first_sample_point;
  }

  /// @cond SKIP_IN_MANUAL
  // This variant deduces the kernel from the iterator type.
  template <typename ForwardIterator,
    typename PointPMap
  >
  ForwardIterator
  regularize_and_simplify_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  double retain_percentage, ///< percentage of points to retain.
  unsigned int k, ///< number of neighbors.
  const unsigned int iter_number, ///< number of iterations.
  const bool need_compute_density  ///< if needed to compute density to generate more rugularized result, 
  ///  especially when the density of input is uneven.
  ) 
  {
    typedef typename boost::property_traits<PointPMap>::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;
    return regularize_and_simplify_point_set(
      first,beyond,
      point_pmap,
      retain_percentage,
      k,
      iter_number,
      need_compute_density,
      Kernel());
  }
  /// @endcond

  /// @cond SKIP_IN_MANUAL
  // This variant creates a default point property map = Dereference_property_map.
  template <typename ForwardIterator 
  >
  ForwardIterator
  regularize_and_simplify_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  double retain_percentage, ///< percentage of points to retain.
  unsigned int k, ///< number of neighbors.
  const unsigned int iter_number, ///< number of iterations.
  const bool need_compute_density ///< if needed to compute density to generate more rugularized result, 
  ///  especially when the density of input is uneven.
  ) 
  {
    return regularize_and_simplify_point_set(
      first,beyond,
      make_dereference_property_map(first),
      retain_percentage, k, iter_number, need_compute_density);
  }
  /// @endcond


} //namespace CGAL

#endif // CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_H
