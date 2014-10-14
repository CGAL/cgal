// Copyright (c) 2014  INRIA Sophia-Antipolis (France)
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Clement Jamin

#ifndef POINT_CLOUD_H
#define POINT_CLOUD_H

#include <CGAL/basic.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Dimension.h>

#ifdef CGAL_TC_USE_NANOFLANN

#include "nanoflann.hpp"

#include <array>
#include <utility>
#include <limits>

namespace CGAL {
namespace Tangential_complex_ {

// "dataset to kd-tree" adaptor class
template <typename K, typename Point_container_>
class Point_cloud_adaptator
{
public:
  typedef typename Point_container_::value_type     Point;
  typedef typename CGAL::Kernel_traits<Point>::type Kernel;
  typedef typename Kernel::FT                       FT;

  /// The constructor that sets the data set source
  Point_cloud_adaptator(Point_container_ &points, Kernel const& k) 
    : m_points(points), m_k(k)
  {}

  /// CRTP helper method
  inline Point_container_ const& points() const 
  { 
    return m_points;
  }  
  inline Point_container_& points() 
  {
    return m_points;
  }

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const 
  {
    return m_points.size();
  }

  // Returns the distance between the vector "p1[0:size-1]" 
  // and the data point with index "idx_p2" stored in the class:
  inline FT kdtree_distance(
    const FT *p1, const size_t idx_p2, size_t size) const
  {
    Point sp(p1, p1 + size);
    return m_k.squared_distance_d_object()(sp, points()[idx_p2]);
  }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an 
  // immediate value, the "if/else's" are actually solved at compile time.
  inline FT kdtree_get_pt(const size_t idx, int dim) const
  {
    return m_k.compute_coordinate_d_object()(points()[idx], dim);
  }

  // Optional bounding-box computation: return false to default to a standard 
  // bbox computation loop.
  // Return true if the BBOX was already computed by the class and returned 
  // in "bb" so it can be avoided to redo it again.
  // Look at bb.size() to find out the expected dimensionality 
  // (e.g. 2 or 3 for point clouds)
  template <class Bbox>
  bool kdtree_get_bbox(Bbox &bb) const
  {
    return false;
  }

  Kernel const& kernel() const
  {
    return m_k;
  }

protected:  
  Point_container_& m_points; //!< A ref to the data set origin
  Kernel const& m_k;      //!< A const ref to the kernel

};

template <typename K, typename Point_container_>
class Point_cloud_data_structure
{
public:
  typedef typename Point_container_::value_type     Point;
  typedef typename K                                Kernel;
  typedef typename Kernel::FT                       FT;

  static const int AMB_DIM = Ambient_dimension<Point>::value; // CJTODO: use Point_dimension_d or similar

  /// Constructor
  Point_cloud_data_structure(Point_container_ &points, Kernel const& k)
  : m_adaptor(points, k),
    m_kd_tree(AMB_DIM, 
              m_adaptor, 
              nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) )
  {
    m_kd_tree.buildIndex();
  }
  
  /*Point_container_ &points()
  {
    return m_adaptor.points();
  }

  const Point_container_ &points() const
  {
    return m_adaptor.points();
  }*/

  void query_ANN(const Point &sp,
    std::size_t k,
    size_t *neighbor_indices,
    FT *squared_distance) const
  {
    /*std::vector<FT> sp_vec(
      m_adaptor.kernel().construct_cartesian_const_iterator_d_object()(sp),
      m_adaptor.kernel().construct_cartesian_const_iterator_d_object()(sp, 0));*/ // CJTODO remettre
    std::vector<FT> sp_vec;
    for (int i = 0 ; i < 4 ; ++i)
      sp_vec.push_back(sp[i]);
    nanoflann::KNNResultSet<FT> result_set(k);
    result_set.init(neighbor_indices, squared_distance);
    m_kd_tree.findNeighbors(result_set,
                            &sp_vec[0], 
                            nanoflann::SearchParams());

    /*std::cout << "knnSearch(nn="<< num_results <<"): \n";
    for (int i = 0 ; i < num_results ; ++i)
    {
      std::cout << "  * neighbor_indices = " << neighbor_indices [i]
                << " (out_dist_sqr = " << squared_distance[i] << ")" 
                << std::endl;
    }*/
  }
  
  void query_ball(const Point &sp, 
                  const FT radius,
                  std::vector<std::pair<std::size_t, FT> > &neighbors,
                  bool sort_output = true)
  {
    /*std::vector<FT> sp_vec(
      m_adaptor.kernel().construct_cartesian_const_iterator_d_object()(sp),
      m_adaptor.kernel().construct_cartesian_const_iterator_d_object()(sp, 0));*/ // CJTODO remettre
    std::vector<FT> sp_vec;
    for (int i = 0 ; i < 4 ; ++i)
      sp_vec.push_back(sp[i]);
    m_kd_tree.radiusSearch(&sp_vec[0],
                           radius,
                           neighbors,
                           nanoflann::SearchParams(32, 0.f, sort_output));

    /*std::cout << "radiusSearch(num="<< neighbors.size() <<"): \n";
    for (const auto idx_and_dist : neighbors)
    {
      std::cout << "  * neighbor_indices = " << idx_and_dist.first
                << " (out_dist_sqr = " << idx_and_dist.second << ")" 
                << std::endl;
    }*/
  }

protected:
  typedef Point_cloud_adaptator<Kernel, Point_container_> Adaptor;
  typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<FT, Adaptor> ,
    Adaptor,
    AMB_DIM // dim
    > Kd_tree;

  Adaptor m_adaptor;
  Kd_tree m_kd_tree;
};

} // namespace Tangential_complex_
} //namespace CGAL

#else // !CGAL_TC_USE_NANOFLANN => use CGAL Spatial searching

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include <utility>
#include <limits>

namespace CGAL {
namespace Tangential_complex_ {

template <typename K, typename Point_container_>
class Point_cloud_data_structure
{
public:
  typedef typename Point_container_::value_type         Point;
  typedef typename K                                    Kernel;
  typedef typename Kernel::FT                           FT;

  typedef CGAL::Search_traits<
    FT, Point,
    typename Kernel::Cartesian_const_iterator_d, 
    typename Kernel::Construct_cartesian_const_iterator_d>  Traits_base;
  // using a pointer as a special property map type
  typedef CGAL::Search_traits_adapter<
    std::ptrdiff_t, Point*, Traits_base>                    STraits;
  
  typedef CGAL::Orthogonal_k_neighbor_search<STraits>       K_neighbor_search;
  typedef typename K_neighbor_search::Tree                  Tree;
  typedef typename K_neighbor_search::Distance              Distance;
  typedef typename K_neighbor_search::iterator              KNS_iterator;
  typedef K_neighbor_search                                 KNS_range;

  typedef CGAL::Orthogonal_incremental_neighbor_search<
    STraits, Distance, CGAL::Sliding_midpoint<STraits>, Tree>      
                                                   Incremental_neighbor_search;
  typedef typename Incremental_neighbor_search::iterator    INS_iterator;
  typedef Incremental_neighbor_search                       INS_range;

  static const int AMB_DIM = Ambient_dimension<Point>::value; // CJTODO: use Point_dimension_d or similar

  /// Constructor
  Point_cloud_data_structure(Point_container_ const& points)
  : m_points(points),
    m_tree(
      boost::counting_iterator<std::ptrdiff_t>(0),
      boost::counting_iterator<std::ptrdiff_t>(points.size()),
      Tree::Splitter(),
      STraits((Point*)&(points[0])) )
  {
  }
  
  /*Point_container_ &points()
  {
    return m_points;
  }

  const Point_container_ &points() const
  {
    return m_points;
  }*/

  KNS_range query_ANN(const 
    Point &sp,
    unsigned int k,
    bool sorted = true) const
  {
    // Initialize the search structure, and search all N points
    // Note that we need to pass the Distance explicitly since it needs to
    // know the property map
    K_neighbor_search search(
      m_tree, 
      sp, 
      k, 
      FT(0), 
      true,
      Distance_adapter<std::ptrdiff_t,Point*,Euclidean_distance<Traits_base> >(
        (Point*)&(m_points[0])),
      sorted);
    
    return search;
  }
  
  INS_range query_incremental_ANN(const Point &sp) const
  {
    // Initialize the search structure, and search all N points
    // Note that we need to pass the Distance explicitly since it needs to
    // know the property map
    Incremental_neighbor_search search(
      m_tree,
      sp,
      FT(0), 
      true,
      Distance_adapter<std::ptrdiff_t,Point*,Euclidean_distance<Traits_base> >(
        (Point*)&(m_points[0])) );
    
    return search;
  }

protected:
  Point_container_ const& m_points;
  Tree m_tree;
};

} // namespace Tangential_complex_
} //namespace CGAL

#endif // CGAL_TC_USE_NANOFLANN

#endif // POINT_CLOUD_H