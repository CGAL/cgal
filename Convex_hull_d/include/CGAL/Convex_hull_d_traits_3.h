// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_CONVEX_HULL_D_TRAITS_3_H
#define CGAL_CONVEX_HULL_D_TRAITS_3_H

#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/predicates_on_points_3.h>
#include <CGAL/intersection_3.h>
#include <vector>

namespace CGAL {

template <class R_> struct Convex_hull_d_traits_3
{
    typedef R_                    R;
    typedef typename R_::RT       RT;
    typedef typename R_::FT       FT;
    typedef typename R_::Point_3  Point_d;
    typedef typename R_::Plane_3  Hyperplane_d;
    typedef typename R_::Vector_3 Vector_d;
    typedef typename R_::Ray_3    Ray_d;

  typedef typename R::Oriented_side_3 Oriented_side_d;
  Oriented_side_d oriented_side_d_object() const
  { return Oriented_side_d(); }

  typedef typename R::Has_on_positive_side_3 Has_on_positive_side_d;
  Has_on_positive_side_d has_on_positive_side_d_object() const
  { return Has_on_positive_side_d(); }

  struct Orientation_d {
    template <class I>
    Orientation operator()(I s, I e) {
      Point_d A[4];
      CGAL_assertion(s != e);
      A[0] = *s;
      ++s;
      CGAL_assertion(s != e);
      A[1] = *s;
      ++s;
      CGAL_assertion(s != e);
      A[2] = *s;
      ++s;
      CGAL_assertion(s != e);
      A[3] = *s;
      ++s;
      CGAL_assertion(s == e);
	
      return orientation(A[0],A[1],A[2],A[3]);
    }
  };
  Orientation_d orientation_d_object() const
  { return Orientation_d(); }

  struct Affinely_independent_d {
    template <class I>
    bool operator()(I s, I e)
    { 
      Point_d A[4];
      CGAL_assertion(s != e);
      A[0] = *s;
      ++s;
      if(s == e){
	return true;
      }
      A[1] = *s;
      ++s;
      if (s == e){
	return A[0] != A[1];
      }
      A[2] = *s;
      ++s;
      if (s == e){
	return ! collinear( A[0], A[1], A[2] );
      }
      A[3] = *s;
      ++s;
      if (s == e){
	return !coplanar( A[0], A[1], A[2], A[3] );
      } 
      return false;
    }
  };


  Affinely_independent_d 
  affinely_independent_d_object() const
  { 
    return Affinely_independent_d(); 
  }


  struct Contained_in_simplex_d {

    template <class I>
    bool 
    operator()(I s, I e, const Point_d& p)
    { 
      Point_d A[4];
      CGAL_assertion(s != e);
      A[0] = *s;
      ++s;
      if(s == e){
	return A[0] == p;
      }
      A[1] = *s;
      ++s;
      if (s == e){
        typename R_::Segment_3 s( A[0], A[1] );
        return s.has_on(p);
      }
      A[2] = *s;
      ++s;
      if (s == e){
        typename R_::Triangle_3 t( A[0], A[1], A[2] );
        return t.has_on(p);
      }
      A[3] = *s;
      ++s;
      if (s == e){
        typename R_::Tetrahedron_3 t( A[0], A[1], A[2], A[3] );
        return !t.has_on_unbounded_side(p);
      }
      return false;  // should be unreachable !
    }
  };

  Contained_in_simplex_d 
  contained_in_simplex_d_object() const
  { 
    return Contained_in_simplex_d(); 
  }
 
  struct Contained_in_affine_hull_d {
    template <class I>
    bool operator()(I s, I e, const Point_d& p)
    { 
      CGAL_assertion_code( Affinely_independent_d affinely_independent; )
      CGAL_assertion( affinely_independent(s,e) );
      Point_d A[3];
      A[0] = *s;
      ++s;
      if(s == e){
	return  p == A[0];
      }
      A[1] = *s;
      ++s;
      if (s == e){
	return collinear( p, A[0], A[1] );
      }
      A[2] = *s;
      ++s;
      if (s == e){
        return coplanar( p, A[0], A[1], A[2] );
      }
      return false;
    }
  };


  Contained_in_affine_hull_d 
  contained_in_affine_hull_d_object() const
  { 
    return Contained_in_affine_hull_d(); 
  }


  struct Component_accessor_d {
    template <typename C>
    int dimension(const C& c) const { return c.dimension(); }
    template <typename C>
    RT homogeneous(const C& c, int i) { return c.homogeneous(i); }
    template <typename C>
    FT cartesian(const C& c, int i) { return c.cartesian(i); }
  };
  Component_accessor_d component_accessor_d_object() const
  { return Component_accessor_d(); }

  struct Orthogonal_vector_d {
    Vector_d operator()(const Hyperplane_d& h) const
    { return h.orthogonal_vector(); }
  };
  Orthogonal_vector_d orthogonal_vector_d_object() const
  { return Orthogonal_vector_d(); }

  struct Point_to_vector_d {
    Vector_d operator()(const Point_d& p) const
    { return p-CGAL::ORIGIN; }
  };
  Point_to_vector_d point_to_vector_d_object() const
  { return Point_to_vector_d(); }

  struct Vector_to_point_d {
    Point_d operator()(const Vector_d& v) const
    { return CGAL::ORIGIN+v; }
  };
  Vector_to_point_d vector_to_point_d_object() const
  { return Vector_to_point_d(); }

  struct Construct_vector_d {
    Vector_d operator()(int, Null_vector) const
    { return Vector_d(NULL_VECTOR); }
  };
  Construct_vector_d construct_vector_d_object() const
  { return Construct_vector_d(); }

  struct Construct_hyperplane_d {
    template <class I>
    Hyperplane_d operator()(I s, I e, const Point_d& p, 
                            Oriented_side side)
    { 
      Hyperplane_d pl;
      Point_d A[3];
      A[0] = *s;
      ++s;
      if(s == e){
	pl = Hyperplane_d( A[0], A[0] - p);
      } else { 
	A[1] = *s;
	++s;
	if(s == e){
	  typename R_::Point_3 hp =
	    A[0] + cross_product( p - A[0], A[1] - A[0] );
	  pl = Hyperplane_d( A[0], A[1], hp );
	} else {
	  A[2] = *s;
	  ++s;
	  if(s == e){
	    pl = Hyperplane_d( A[0], A[1], A[2] );
	  } else {
	    CGAL_error();
	  }
	}
      }
      if (side != 0) {
        if ( pl.oriented_side(p) !=  side ) { pl = pl.opposite(); }
      }
      return pl;
    }
  };
  Construct_hyperplane_d construct_hyperplane_d_object() const
  { return Construct_hyperplane_d(); }

  typedef typename R::Intersect_3 Intersect_d;
  Intersect_d intersect_d_object() const 
  { return Intersect_d(); }

};

} //namespace CGAL

#endif // CGAL_CONVEX_HULL_D_TRAITS_3_H
