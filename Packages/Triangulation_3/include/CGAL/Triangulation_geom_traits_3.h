// ============================================================================
//
// Copyright (c) 1998-1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Triangulation_geom_traits_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (Mariette Yvinec)
//
// ============================================================================
//
// geometric traits for a <=3 D triangulation
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_GEOM_TRAITS_3_H
#define CGAL_TRIANGULATION_GEOM_TRAITS_3_H

#include <CGAL/Point_3.h>
#include <CGAL/Point_2.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Tetrahedron_3.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_3.h>

template < class R >
class CGAL_Triangulation_geom_traits_3 
{
public:

  typedef CGAL_Point_3<R>  Point;
  typedef CGAL_Point_2< CGAL_Cartesian< CGAL_Quotient<typename R::RT> > >  Point2;
  typedef CGAL_Segment_3<R> Segment;
  typedef CGAL_Triangle_3<R> Triangle;
  typedef CGAL_Tetrahedron_3<R> Tetrahedron;

  CGAL_Triangulation_geom_traits_3()
    {}

  CGAL_Triangulation_geom_traits_3(const CGAL_Triangulation_geom_traits_3 & gt)
     {}

  CGAL_Triangulation_geom_traits_3 & 
  operator=(const CGAL_Triangulation_geom_traits_3 & gt)
    {return *this;}

  // PREDICATES ON POINTS

  bool equal(const Point & p, const Point & q) const
  {
    return (p == q);
  }
  
  CGAL_Comparison_result compare_x(const Point & p, const Point & q) const
    {
      return CGAL_compare_x(p, q);
    }


  CGAL_Comparison_result compare_y(const Point & p, const Point & q) const
    {
      return CGAL_compare_y(p, q);
    }

  CGAL_Comparison_result compare_z(const Point & p, const Point & q) const
    {
      return CGAL_compare_z(p, q);
    }

  CGAL_Orientation orientation(const Point & p,
			       const Point & q,
			       const Point & r,
			       const Point & s) const
  {
    return CGAL_orientation(p, q, r, s);
  }

  CGAL_Orientation orientation_in_plane(const Point & q,
					const Point & r,
					const Point & s,
					const Point & p) const
    // p,q,r,s supposed to be coplanar
    // q,r,s supposed to be non collinear
    // tests whether p is on the same side of q,r as s
    // returns :
    // CGAL_COLLINEAR if pqr collinear
    // CGAL_POSITIVE if qrp and qrs have the same orientation
    // CGAL_NEGATIVE if qrp and qrs have opposite orientations
  {
    // should be used only when p,q,r,s are coplanar
    CGAL_triangulation_precondition( ! collinear(q,r,s) );
    CGAL_triangulation_precondition( orientation(p,q,r,s) == CGAL_COPLANAR );
    // projection on the x,y-plane
    Point2 P(p.hx(), p.hy(), p.hw());
    Point2 Q(q.hx(), q.hy(), q.hw());
    Point2 R(r.hx(), r.hy(), r.hw());
    Point2 S(s.hx(), s.hy(), s.hw());
    CGAL_Orientation oxy_qrs = CGAL_orientation(Q,R,S);

    if ( oxy_qrs != CGAL_COLLINEAR )
      // the projection on x,y is OK
      return CGAL_Orientation( oxy_qrs * CGAL_orientation(Q,R,P) );

    // else : must project on another plane
    // tests on which plane :

    if ( ( Q.x() != R.x() ) || 
	 ( Q.x() != S.x() ) ) {
      // projection on x,z-plane is ok
      P = Point2(p.hx(), p.hz(), p.hw());
      Q = Point2(q.hx(), q.hz(), q.hw());
      R = Point2(r.hx(), r.hz(), r.hw());
      S = Point2(s.hx(), s.hz(), s.hw());
    }
    else
    { // projection on y,z-plane
      P = Point2(p.hy(), p.hz(), p.hw());
      Q = Point2(q.hy(), q.hz(), q.hw());
      R = Point2(r.hy(), r.hz(), r.hw());
      S = Point2(s.hy(), s.hz(), s.hw());
    }

    return CGAL_Orientation ( CGAL_orientation(Q,R,S)
                            * CGAL_orientation(Q,R,P) );
  }

  bool collinear(const Point & p,
		 const Point & q,
		 const Point & r) const
    {
      return CGAL_collinear(p,q,r);
    }

  // DELAUNAY

  CGAL_Oriented_side 
  side_of_oriented_sphere(const Point & p,
			  const Point & q,
			  const Point & r,
			  const Point & s,
			  const Point & test) const
    {
      return CGAL_side_of_oriented_sphere(p, q, r, s, test);
    }

  CGAL_Oriented_side 
  side_of_oriented_circle(const Point & p,
			  const Point & q,
			  const Point & r,
			  const Point & test) const
    {
      CGAL_triangulation_precondition( orientation(p,q,r,test) ==
				       CGAL_COPLANAR );

      // test belongs to the circle if and only if it belongs to a
      // sphere passing through pqr
      CGAL_Orientation or;
      Point O(0,0,0), A(1,0,0), B(0,1,0), C(0,0,1);

      Point P = ((or = orientation(p,q,r,O)) != CGAL_ZERO) ? O:
                ((or = orientation(p,q,r,A)) != CGAL_ZERO) ? A:
                ((or = orientation(p,q,r,B)) != CGAL_ZERO) ? B:
                ((or = orientation(p,q,r,C)) != CGAL_ZERO) ? C: C;

      return CGAL_Oriented_side( or *
	      CGAL_side_of_oriented_sphere(p, q, r, P, test));
    }

};


#endif CGAL_TRIANGULATION_GEOM_TRAITS_3_H
