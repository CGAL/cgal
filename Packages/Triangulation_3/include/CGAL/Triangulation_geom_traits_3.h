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
// coordinator   : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
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
  typedef CGAL_Point_2<R>  Point2;
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

  CGAL_Orientation orientation_in_plane(const Point & p,
					const Point & q,
					const Point & r,
					const Point & s) const
    // p,q,r,s supposed to be coplanar
    // q,r,s supposed to be a positively oriented triangle
    // tests whether p is on the same side of q,r as s
    // returns :
    // CGAL_COLLINEAR if pqr collinear
    // CGAL_POSITIVE if qrp and qrs have the same orientation
    // CGAL_NEGATIVE if qrp and qrs have opposite orientations
  {
    // no precondition but should be used only when p,q,r,s are coplanar

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
	if ( oxy_qrp == CGAL_COLLINEAR ) { return CGAL_COLLINEAR;	}
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
};


#endif CGAL_TRIANGULATION_GEOM_TRAITS_3_H
