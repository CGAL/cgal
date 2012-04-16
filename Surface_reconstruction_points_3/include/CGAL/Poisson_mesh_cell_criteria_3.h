// Copyright (c) 2004-2005  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_POISSON_MESH_CRITERIA_3_H
#define CGAL_POISSON_MESH_CRITERIA_3_H

#include <iostream>

namespace CGAL {

template <class Tr>
class Poisson_mesh_cell_criteria_3 
{
  double squared_radius_bound_;
  double radius_edge_bound_;
public:
  struct Cell_quality : public std::pair<double, double>
  {
    typedef std::pair<double, double> Base;

    Cell_quality() : Base() {}
    Cell_quality(double _aspect, double _sq_size) : Base(_aspect, _sq_size) {};

    double sq_size() const { return second; }
    double aspect() const { return first; }

    // q1<q2 means q1 is prioritised over q2
    // ( q1 == *this, q2 == q )
    bool operator<(const Cell_quality& q) const
    {
      if( sq_size() > 1 )
	if( q.sq_size() > 1 )
	  return ( sq_size() > q.sq_size() );
	else
	  return true; // *this is big but not q
      else
	if( q.sq_size() >  1 )
	  return false; // q is big but not *this
      return( aspect() > q.aspect() );
    }
  };

  inline
  double squared_radius_bound() const 
  {
    return squared_radius_bound_; 
  }

  typedef typename Tr::Cell_handle Cell_handle;

  Poisson_mesh_cell_criteria_3(const double radius_edge_bound = 2, //< radius edge ratio bound (ignored if zero)
		  const double radius_bound = 0) //< cell radius bound (ignored if zero)
    : squared_radius_bound_(radius_bound*radius_bound),
      radius_edge_bound_(radius_edge_bound)
  {}

  inline 
  void set_squared_radius_bound(const double squared_radius_bound) 
  { 
    squared_radius_bound_ = squared_radius_bound;
  }

  inline
  double radius_edge_bound() const 
  {
    return radius_edge_bound_; 
  }

  inline 
  void set_radius_edge_bound(const double radius_edge_bound) 
  { 
    radius_edge_bound_ = radius_edge_bound;
  }

  class Is_bad
  {
  protected:
    const double radius_edge_bound_;
    const double squared_radius_bound_;
  public:
    typedef typename Tr::Point Point_3;
      
    Is_bad(const double radius_edge_bound, 
	   const double squared_radius_bound)
      : radius_edge_bound_(radius_edge_bound),
	squared_radius_bound_(squared_radius_bound) {}
      
    bool operator()(const Cell_handle& c,
                    Cell_quality& qual) const
    {
      const Point_3& p = c->vertex(0)->point();
      const Point_3& q = c->vertex(1)->point();
      const Point_3& r = c->vertex(2)->point();
      const Point_3& s = c->vertex(3)->point();

      typedef typename Tr::Geom_traits Geom_traits;
      typedef typename Geom_traits::Compute_squared_radius_3 Radius;
      typedef typename Geom_traits::Compute_squared_distance_3 Distance;
      typedef typename Geom_traits::FT FT;

      Radius radius = Geom_traits().compute_squared_radius_3_object();
      Distance distance = Geom_traits().compute_squared_distance_3_object();

      double size = to_double(radius(p, q, r, s));

      if( squared_radius_bound_ != 0 )
        {
          qual.second = size / squared_radius_bound_;
            // normalized by size bound to deal
            // with size field
          if( qual.sq_size() > 1 )
            {
              qual.first = 1; // (do not compute aspect)
              return true;
            }
        }
      if( radius_edge_bound_ == 0 )
	{
	  qual = Cell_quality(0,1);
	  return false;
	}

      double min_sq_length = CGAL::to_double(distance(p, q));
      min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(p, r)));
      min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(p, s)));
      min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(q, r)));
      min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(q, s)));
      min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(r, s)));

      qual.first = size / min_sq_length;

      return (qual.first > radius_edge_bound_);
    }
    
  }; // end Is_bad
    

  Is_bad is_bad_object() const
  { return Is_bad(radius_edge_bound_, squared_radius_bound_); }

}; // end Poisson_mesh_cell_criteria_3

  template <typename Tr>
  std::ostream& operator<<(std::ostream& os,
                           const typename Poisson_mesh_cell_criteria_3<Tr>::Cell_quality& q)
  {
    return os << q.sq_size() << ", " << q.aspect();
  }

} // end namespace CGAL

#endif // CGAL_POISSON_MESH_CRITERIA_3_H
