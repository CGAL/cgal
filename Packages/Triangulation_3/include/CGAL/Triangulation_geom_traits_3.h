// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
//   CGAL_Triangulation_geom_traits_3(const CGAL_Triangulation_geom_traits_3 & gt)
//     {}

//   CGAL_Triangulation_geom_traits_3 & 
//   operator=(const CGAL_Triangulation_geom_traits_3 & gt)
//     {return *this;}

  // PREDICATES ON POINTS

  bool equal(const Point & p,
	     const Point & q) const
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
    CGAL_triangulation_precondition( orientation(p,q,r,s) ==
				     CGAL_COPLANAR );
    // projection on the x,y-plane
    Point2 pxy(p.hx(), p.hy(), p.hw());
    Point2 qxy(q.hx(), q.hy(), q.hw());
    Point2 rxy(r.hx(), r.hy(), r.hw());
    Point2 sxy(s.hx(), s.hy(), s.hw());
    CGAL_Orientation oxy_qrs = CGAL_orientation(qxy,rxy,sxy);

    if ( oxy_qrs !=  CGAL_COLLINEAR ) {
      // the projection on x,y is OK
      // tests whether pxy is on the same side of qxy, rxy as sxy
      CGAL_Orientation oxy_qrp = CGAL_orientation(qxy,rxy,pxy);
      if ( oxy_qrp == oxy_qrs) {
	return CGAL_POSITIVE;
      }
      else {
	if ( oxy_qrp == CGAL_COLLINEAR ) { return CGAL_COLLINEAR; }
	else { return CGAL_NEGATIVE; }
      }
    }

    // else : must project on another plane
    // tests on which plane :

    if ( ( qxy.x() != rxy.x() ) || 
	 ( qxy.x() != sxy.x() ) ) {
      // projection on x,z-plane is ok
      Point2 pxz(p.hx(), p.hz(), p.hw());
      Point2 qxz(q.hx(), q.hz(), q.hw());
      Point2 rxz(r.hx(), r.hz(), r.hw());
      Point2 sxz(s.hx(), s.hz(), s.hw());
      // tests whether pxz is on the same side of qxz, rxz as sxz
      CGAL_Orientation oxz_qrs = CGAL_orientation(qxz,rxz,sxz);
      CGAL_Orientation oxz_qrp = CGAL_orientation(qxz,rxz,pxz);
      if ( oxz_qrp == oxz_qrs) {
	return CGAL_POSITIVE;
      }
      else {
	if ( oxz_qrp == CGAL_COLLINEAR ) { return CGAL_COLLINEAR;	}
	else { return CGAL_NEGATIVE; }
      }
    }
   
    // else : projection on y,z-plane
    Point2 pyz(p.hy(), p.hz(), p.hw());
    Point2 qyz(q.hy(), q.hz(), q.hw());
    Point2 ryz(r.hy(), r.hz(), r.hw());
    Point2 syz(s.hy(), s.hz(), s.hw());
    // tests whether pyz is on the same side of qyz, ryz as syz
    CGAL_Orientation oyz_qrs = CGAL_orientation(qyz,ryz,syz);
    CGAL_Orientation oyz_qrp = CGAL_orientation(qyz,ryz,pyz);
    if ( oyz_qrp == oyz_qrs) {
      return CGAL_POSITIVE;
    }
    else {
      if ( oyz_qrp == CGAL_COLLINEAR ) { return CGAL_COLLINEAR;	}
      else { return CGAL_NEGATIVE; }
    }
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
      Point O(0,0,0);
      switch ( orientation( p,q,r,O ) ) {
      case CGAL_POSITIVE:
	return CGAL_side_of_oriented_sphere(p, q, r, O, test); 
      case CGAL_NEGATIVE:
	return CGAL_side_of_oriented_sphere(O, p, q, r, test);
      default: 
	break;
      }
      // if O coplanar, use A
      Point A(1,0,0);
      switch ( orientation( p,q,r,A ) ) {
      case CGAL_POSITIVE:
	return CGAL_side_of_oriented_sphere(p, q, r, A, test); 
      case CGAL_NEGATIVE:
	return CGAL_side_of_oriented_sphere(A, p, q, r, test);
      default: 
	break;
      }
      // if A is coplanar, use B
      Point B(0,1,0);
      switch ( orientation( p,q,r,B ) ) {
      case CGAL_POSITIVE:
	return CGAL_side_of_oriented_sphere(p, q, r, B, test); 
      case CGAL_NEGATIVE:
	return CGAL_side_of_oriented_sphere(B, p, q, r, test);
      default: 
	break;
     }
      // if B also coplanar, use C
      Point C(0,0,1);
      switch ( orientation( p,q,r,C ) ) {
      case CGAL_POSITIVE:
	return CGAL_side_of_oriented_sphere(p, q, r, C, test); 
      case CGAL_NEGATIVE:
	return CGAL_side_of_oriented_sphere(C, p, q, r, test);
      default: 
	break;
      }
      // impossible, only to avoid compilation warnings :
      CGAL_triangulation_assertion(false);
      return CGAL_ON_POSITIVE_SIDE;
    }

};


#endif CGAL_TRIANGULATION_GEOM_TRAITS_3_H
