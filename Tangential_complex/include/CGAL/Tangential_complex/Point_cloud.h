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
#include <CGAL/Kernel_traits.h>
#include <CGAL/Dimension.h>

#include "nanoflann.hpp"

#include <array>
#include <utility>
#include <limits>

namespace CGAL {
namespace Tangential_complex_ {

template<typename Point_>
class Point_cloud
: public std::vector<Point_>
{
public:
  typedef std::vector<Point_>                       Base;
  typedef Base                                      Raw_container;
  typedef Point_                                    Point;
  typedef typename CGAL::Kernel_traits<Point>::type Kernel;
  typedef typename Kernel::FT                       FT;
  
  //typedef typename Base::iterator                   iterator;
  //typedef typename Base::const_iterator             const_iterator;
  
  static const int AMB_DIM = Ambient_dimension<Point>::value;

  Point_cloud(Kernel const& k)
  : m_k(k)
  {
    m_mins.fill(std::numeric_limits<FT>::max());
    m_maxs.fill(std::numeric_limits<FT>::min());
  }

  template <class InputIterator>
  Point_cloud(InputIterator first, InputIterator last, Kernel const& k)
  : Base(first, last), m_k(k)
  {
    m_mins.fill(std::numeric_limits<FT>::max());
    m_maxs.fill(std::numeric_limits<FT>::min());
  }

  void push_back(const Point &point, bool update_bbox = false)
  {
    Base::push_back(point);

    // Adjust bbox?
    if (update_bbox)
    {
      for (int i = 0 ; i < AMB_DIM ; ++i)
      {
        if (point.get_param(i) < m_mins[i])
          m_mins[i] = point.get_param(i);
      
        if (point.get_param(i) > m_maxs[i])
          m_maxs[i] = point.get_param(i);
      }
    }
  }

  void compute_bbox(bool only_if_not_already_done = false)
  {
    if (only_if_not_already_done && m_mins[0] != std::numeric_limits<FT>::max())
      return;

    // Reset
    m_mins.fill(std::numeric_limits<FT>::max());
    m_maxs.fill(std::numeric_limits<FT>::min());

    // Adjust bbox
    for (const auto &point : *this)
    {
      typedef typename Kernel::Compute_coordinate_d Ccd;
      const Ccd ccd = m_k.compute_coordinate_d_object();
      for (int i = 0 ; i < AMB_DIM ; ++i)
      {
        if (ccd(point, i) < m_mins[i])
          m_mins[i] =ccd(point, i);
      
        if (ccd(point, i) > m_maxs[i])
          m_maxs[i] = ccd(point, i);
      }
    }
  }

  FT get_min(int dim) const
  {
    return m_mins[dim];
  }
  
  FT get_max(int dim) const
  {
    return m_maxs[dim];
  }

  FT bbox_diagonal() const
  {
    FT sqdiag = 0;
    for (std::size_t i = 0 ; i < AMB_DIM ; ++i)
    {
      FT d = m_maxs[i] - m_mins[i];
      sqdiag += d*d;
    }
    return CGAL::sqrt(sqdiag);
  }

  void recenter_points_around_origin(
    bool compute_bbox_if_not_already_done = true)
  {
    // If the bounding box has not been computed already
    if (m_mins[0] == std::numeric_limits<FT>::max())
      if (compute_bbox_if_not_already_done)
        compute_bbox();
      else
        return;

    // Compute centre of bbox
    std::array<FT, AMB_DIM> transl_array;
    for (std::size_t i = 0 ; i < AMB_DIM ; ++i)
      transl_array[i] = -0.5*(m_maxs[i] + m_mins[i]);

    Point transl(transl_array);

#ifdef CGAL_LINKED_WITH_TBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, size()),
      [&]( const tbb::blocked_range<size_t>& r )
    {
      for (auto i = r.begin(); i != r.end(); ++i)
#else
    for (auto i = 0; i != size(); ++i)
#endif
    {
        (*this)[i] += transl;
    }
#ifdef CGAL_LINKED_WITH_TBB
  });
#endif

  }

protected:
  Kernel const& m_k;      //!< A const ref to the kernel
  // Bounding box
  std::array<FT, AMB_DIM> m_mins;
  std::array<FT, AMB_DIM> m_maxs;
};


// And this is the "dataset to kd-tree" adaptor class:
template <typename Point_cloud_>
class Point_cloud_adaptator
{
public:
  typedef typename Point_cloud_::Kernel Kernel;
  typedef typename Point_cloud_::Point  Point;
  typedef typename Point_cloud_::FT     FT;

  /// The constructor that sets the data set source
  Point_cloud_adaptator(Point_cloud_ &point_cloud, Kernel const& k) 
    : m_points(point_cloud), m_k(k)
  {}

  /// CRTP helper method
  inline Point_cloud_ const& point_cloud() const 
  { 
    return m_points;
  }  
  inline Point_cloud_& point_cloud() 
  {
    return m_points;
  }

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const 
  {
    return point_cloud().size();
  }

  // Returns the distance between the vector "p1[0:size-1]" 
  // and the data point with index "idx_p2" stored in the class:
  inline FT kdtree_distance(
    const FT *p1, const size_t idx_p2, size_t size) const
  {
    Point sp(p1, p1 + size);
    return m_k.squared_distance_d_object()(sp, point_cloud()[idx_p2]);
  }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an 
  // immediate value, the "if/else's" are actually solved at compile time.
  inline FT kdtree_get_pt(const size_t idx, int dim) const
  {
    return m_k.compute_coordinate_d_object()(point_cloud()[idx], dim);
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
    for (int i = 0 ; i < bb.size() ; ++i)
    {
      bb[i].low = m_points.get_min(i);
      bb[i].high = m_points.get_max(i);
    }

    return true;
  }

  Kernel const& kernel() const
  {
    return m_k;
  }

protected:  
  Point_cloud_& m_points; //!< A ref to the data set origin
  Kernel const& m_k;      //!< A const ref to the kernel

}; // end of PointCloudAdaptor

template <typename Point_cloud_>
class Point_cloud_data_structure
{
public:
  typedef typename Point_cloud_::Kernel Kernel;
  typedef typename Point_cloud_::Point  Point;
  typedef typename Point_cloud_::FT     FT;

  static const int AMB_DIM = Ambient_dimension<Point>::value;

  /// Constructor
  Point_cloud_data_structure(Point_cloud_ &cloud, Kernel const& k)
  : m_adaptor(cloud, k),
    m_kd_tree(AMB_DIM, 
              m_adaptor, 
              nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) )
  {
    cloud.compute_bbox(true);
    //cloud.recenter_points_around_origin();
    m_kd_tree.buildIndex();
  }
  
  Point_cloud_ &point_cloud()
  {
    return m_adaptor.point_cloud();
  }

  const Point_cloud_ &point_cloud() const
  {
    return m_adaptor.point_cloud();
  }

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
  typedef Point_cloud_adaptator<Point_cloud_> Adaptor;
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

#endif // POINT_CLOUD_H