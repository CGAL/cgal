// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// 
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================
//
// geometric traits for a <=3 D triangulation
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_GEOM_TRAITS_3_H
#define CGAL_TRIANGULATION_GEOM_TRAITS_3_H

#include <CGAL/basic.h>

#include <CGAL/Point_3.h>
#include <CGAL/Point_2.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Tetrahedron_3.h>

#include <CGAL/triangulation_assertions.h>

#include <CGAL/Triangulation_short_names_3.h>

CGAL_BEGIN_NAMESPACE

template < class Repres >
class Triangulation_geom_traits_3 
{
public:

  typedef Point_3<Repres>  Point;
  typedef Point_2< Cartesian< typename Repres::FT> >  Point2;
  typedef Segment_3<Repres> Segment;
  typedef Triangle_3<Repres> Triangle;
  typedef Tetrahedron_3<Repres> Tetrahedron;

  Triangulation_geom_traits_3()
    {}

  Triangulation_geom_traits_3(const Triangulation_geom_traits_3 & )
     {}

  Triangulation_geom_traits_3 & 
  operator=(const Triangulation_geom_traits_3 & )
    {return *this;}

  // PREDICATES ON POINTS

  bool equal(const Point & p, const Point & q) const
  {
    return ( CGAL::compare_x(p, q)== EQUAL &&  
	     CGAL::compare_y(p, q)== EQUAL &&
	     CGAL::compare_z(p, q)== EQUAL );
  }
  
  Comparison_result compare_x(const Point & p, const Point & q) const
    {
      return CGAL::compare_x(p, q);
    }


  Comparison_result compare_y(const Point & p, const Point & q) const
    {
      return CGAL::compare_y(p, q);
    }

  Comparison_result compare_z(const Point & p, const Point & q) const
    {
      return CGAL::compare_z(p, q);
    }

  Orientation orientation(const Point & p,
			  const Point & q,
			  const Point & r,
			  const Point & s) const
  {
    return CGAL::orientation(p, q, r, s);
  }

  Orientation orientation_in_plane(const Point & q,
					const Point & r,
					const Point & s,
					const Point & p) const
    // p,q,r,s supposed to be coplanar
    // q,r,s supposed to be non collinear
    // tests whether p is on the same side of q,r as s
    // returns :
    // COLLINEAR if pqr collinear
    // POSITIVE if qrp and qrs have the same orientation
    // NEGATIVE if qrp and qrs have opposite orientations
  {
    // should be used only when p,q,r,s are coplanar
    CGAL_triangulation_precondition( ! CGAL::collinear(q,r,s) );
    CGAL_triangulation_precondition( CGAL::orientation(p,q,r,s) == COPLANAR );
    // projection on the x,y-plane
    Point2 P(p.x(), p.y());
    Point2 Q(q.x(), q.y());
    Point2 R(r.x(), r.y());
    Point2 S(s.x(), s.y());
    Orientation oxy_qrs = CGAL::orientation(Q,R,S);

    if ( oxy_qrs != COLLINEAR )
      // the projection on x,y is OK
      return Orientation( oxy_qrs * CGAL::orientation(Q,R,P) );

    // else : must project on another plane
    // tests on which plane :

    if ( ( Q.x() != R.x() ) || 
	 ( Q.x() != S.x() ) ) {
      // projection on x,z-plane is ok
      P = Point2(p.x(), p.z());
      Q = Point2(q.x(), q.z());
      R = Point2(r.x(), r.z());
      S = Point2(s.x(), s.z());
    }
    else
    { // projection on y,z-plane
      P = Point2(p.y(), p.z());
      Q = Point2(q.y(), q.z());
      R = Point2(r.y(), r.z());
      S = Point2(s.y(), s.z());
    }

    return Orientation ( CGAL::orientation(Q,R,S) * CGAL::orientation(Q,R,P) );
  }

  bool collinear(const Point & p,
		 const Point & q,
		 const Point & r) const
    {
      return CGAL::collinear(p,q,r);
    }

  // DELAUNAY

  Oriented_side 
  side_of_oriented_sphere(const Point & p,
			  const Point & q,
			  const Point & r,
			  const Point & s,
			  const Point & test) const
    {
      return CGAL::side_of_oriented_sphere(p, q, r, s, test);
    }

  Oriented_side 
  side_of_oriented_circle(const Point & p,
			  const Point & q,
			  const Point & r,
			  const Point & test) const
    {
      CGAL_triangulation_precondition( CGAL::orientation(p,q,r,test) ==
				       COPLANAR );
      CGAL_triangulation_precondition( ! CGAL::collinear(p,q,r) );

      // test belongs to the circle if and only if it belongs to a
      // sphere passing through pqr
      Orientation or;
      Point O(0,0,0), A(1,0,0), B(0,1,0), C(0,0,1);

      Point P = ((or = CGAL::orientation(p,q,r,O)) != ZERO) ? O:
                ((or = CGAL::orientation(p,q,r,A)) != ZERO) ? A:
                ((or = CGAL::orientation(p,q,r,B)) != ZERO) ? B:
                ((or = CGAL::orientation(p,q,r,C)) != ZERO) ? C: C;

      return Oriented_side( or *
	      CGAL::side_of_oriented_sphere(p, q, r, P, test));
    }

};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_GEOM_TRAITS_3_H
