// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_CRITERIA_3_H
#define CGAL_MESH_CRITERIA_3_H

namespace CGAL {

template <typename Tr>
class Mesh_criteria_3 
{
  double size_bound;
  double shape_bound;
public:
  typedef typename Tr::Cell_handle Cell_handle;

  Mesh_criteria_3(const double radius_edge_bound = 2, 
		  const double squared_radius_bound = 0)
    : size_bound(squared_radius_bound),
      shape_bound(radius_edge_bound)
  {
  };

  typedef double Quality;

  inline
  double squared_radius_bound() const 
  {
    return size_bound; 
  };

  inline 
  void set_squared_radius_bound(const double squared_radius_bound) 
  { 
    size_bound = squared_radius_bound;
  };

  inline
  double radius_edge_bound() const 
  {
    return shape_bound; 
  };

  inline 
  void set_radius_edge_bound(const double radius_edge_bound) 
  { 
    shape_bound = radius_edge_bound;
  };

  class Is_bad
  {
  protected:
    const double shape_bound;
    const double size_bound;
  public:
    typedef typename Tr::Point Point_3;
      
    Is_bad(const double radius_edge_bound, 
	   const double squared_radius_bound)
      : shape_bound(radius_edge_bound),
	size_bound(squared_radius_bound) {};
      
    bool operator()(const Cell_handle& c,
                    Quality& qual) const
    {
      typedef typename Tr::Geom_traits Geom_traits;
      typedef typename Geom_traits::Compute_squared_radius_3 Radius;

      if( size_bound==0 )
        {
          qual = 1;
          return false;
        }

      const Point_3& p = c->vertex(0)->point();
      const Point_3& q = c->vertex(1)->point();
      const Point_3& r = c->vertex(2)->point();
      const Point_3& s = c->vertex(3)->point();

      Radius radius = Geom_traits().compute_squared_radius_3_object();
      qual = size_bound / to_double(radius(p, q, r, s));
      return qual < 1.;
    };
  };

  Is_bad is_bad_object() const
    { return Is_bad(shape_bound, size_bound); }

}; // end Mesh_criteria_3
  

} // end namespace CGAL

#endif // CGAL_MESH_CRITERIA_3_H
