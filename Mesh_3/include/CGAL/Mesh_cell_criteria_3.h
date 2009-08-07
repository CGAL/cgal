// Copyright (c) 2004-2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
//
// Author(s)     : Laurent RINEAU, Stephane Tayeb


#ifndef CGAL_MESH_CELL_CRITERIA_3_H
#define CGAL_MESH_CELL_CRITERIA_3_H

#include <iostream>
#include <CGAL/Mesh_3/mesh_standard_cell_criteria.h>

namespace CGAL {
  
template <class Tr>
class Mesh_default_cells_criteria_3
{
public:
  typedef typename Tr::Geom_traits::FT FT;
  
  /**
   * @struct Quality
   * @brief Represents the quality of a cell
   */
  struct Cell_quality
  {
    Cell_quality(const FT& aspect, const FT& sq_size)
    : aspect_(aspect)
    , sq_size_(sq_size) {};
    
    // q1<q2 means q1 is prioritized over q2
    // ( q1 = *this, q2 = q )
    bool operator<(const Cell_quality& q) const
    {
      if (sq_size_ > 1)
        if (q.sq_size_ > 1)
          return (sq_size_ > q.sq_size_);
        else
          return true; // *this is big but not q
        else if (q.sq_size_ > 1)
          return false; // q is big but not *this
      
      return (aspect_ > q.aspect_);
    }
    
  private:
    FT aspect_;
    FT sq_size_;
  };
  
  typedef typename Tr::Cell_handle Cell_handle;
  typedef boost::optional<Cell_quality> Cell_badness;
  
  /**
   * @brief Constructor
   * @param radius_edge_bound the radius-edge bound
   * @param radius_bound the radius bound (tet sizing)
   */
  Mesh_default_cells_criteria_3(const FT radius_edge_bound,
                                const FT radius_bound)
  : radius_edge_bound_(radius_edge_bound)
  , squared_radius_bound_(radius_bound*radius_bound) { };
  
  /// Destructor
  ~Mesh_default_cells_criteria_3() {};
  
  /**
   * @brief returns the badness of cell \c cell
   * @param cell the cell
   * @return the badness of \c cell
   */
  Cell_badness operator()(const Cell_handle& cell) const;
  
private:
  FT radius_edge_bound_;
  FT squared_radius_bound_;
  
};  // end class Mesh_default_cells_criteria


template <class Tr>
typename Mesh_default_cells_criteria_3<Tr>::Cell_badness
Mesh_default_cells_criteria_3<Tr>::operator()(const Cell_handle& cell) const
{
  typedef typename Tr::Point Point_3;
  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Geom_traits::Compute_squared_radius_3 Radius;
  typedef typename Geom_traits::Compute_squared_distance_3 Distance;
  
  const Point_3& p = cell->vertex(0)->point();
  const Point_3& q = cell->vertex(1)->point();
  const Point_3& r = cell->vertex(2)->point();
  const Point_3& s = cell->vertex(3)->point();
  
  Radius radius = Geom_traits().compute_squared_radius_3_object();
  Distance distance = Geom_traits().compute_squared_distance_3_object();
  
  double size = to_double(radius(p, q, r, s));
  
  FT qual_sq_size = 0;
  if( 0 != squared_radius_bound_ )
  {
    qual_sq_size = size / squared_radius_bound_;
    // normalized by size bound to deal with size field
    if( qual_sq_size > 1 )
    {
      // (do not compute aspect)
#ifdef CGAL_MESH_3_DEBUG_CELL_CRITERIA
      std::cerr << "bad cell (radius bound): size[" << size
      << "] bound[" << squared_radius_bound_ << "]\n" ;
#endif
      return Cell_badness(Cell_quality(1,qual_sq_size));
    }
  }
  
  if( 0 == radius_edge_bound_ )
  {
    return Cell_badness();
  }
  
  double min_sq_length = to_double(distance(p, q));
  min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(p, r)));
  min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(p, s)));
  min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(q, r)));
  min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(q, s)));
  min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(r, s)));
  
  const FT qual_aspect = size / min_sq_length;
  
  if ( qual_aspect > radius_edge_bound_ )
  {
#ifdef CGAL_MESH_3_DEBUG_CELL_CRITERIA
    std::cerr << "bad cell (radius-edge bound): radius-edge["
    << qual_aspect << "] bound[" << radius_edge_bound_
    << "]\n" ;
#endif
    return Cell_badness(Cell_quality(qual_aspect,qual_sq_size));
  }
  else
  {
    return Cell_badness();
  }
}
  
  
  
template <typename Tr, typename Visitor_ = Mesh_3::Cell_criterion_visitor<Tr> >
class Mesh_cell_criteria_3
{
  typedef Visitor_ Visitor;
  typedef Mesh_3::Criteria<Tr,Visitor> Criteria;
  
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Geom_traits::FT FT;
  
  typedef Mesh_cell_criteria_3<Tr> Self;
  
public:
  typedef typename Visitor::Cell_quality Cell_quality;
  typedef typename Visitor::Cell_badness Cell_badness;
  
  /**
   * @brief Constructor
   * @param radius_edge_bound the radius-edge bound
   * @param radius_bound the radius bound (tet sizing)
   */
  Mesh_cell_criteria_3(const FT& radius_edge_bound,
                       const FT& radius_bound)
  {
    typedef Mesh_3::Cell_radius_criterion<Tr,Visitor> Radius_criterion;
    typedef Mesh_3::Cell_radius_edge_criterion<Tr,Visitor> Radius_edge_criterion;
    
    if ( 0 != radius_bound )
      criteria_.add(new Radius_criterion(radius_bound));
    
    if ( 0 != radius_edge_bound )
      criteria_.add(new Radius_edge_criterion(radius_edge_bound));
  }
  
  /// Destructor
  ~Mesh_cell_criteria_3() { }
  
  /**
   * @brief returns the badness of cell \c cell
   * @param cell the cell
   * @return the badness of \c cell
   */
  Cell_badness operator()(const Cell_handle& cell) const
  {
    return criteria_(cell);
  }
  
private:
  Criteria criteria_;
  
};  // end class Mesh_cell_criteria_3

}  // end namespace CGAL


#endif // CGAL_MESH_CELL_CRITERIA_3_H

