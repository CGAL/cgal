// ============================================================================
//
// Copyright (c) 1997  The CGAL Consortium
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
// file          : include/CGAL/Triangulation_euclidean_traits_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
#define CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_2_H

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/basic_constructions_2.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Distance_2.h>

CGAL_BEGIN_NAMESPACE 

template < class R >
class Triangulation_euclidean_traits_2 {
public:
  typedef R Rep;
  typedef Point_2<R>  Point;
  typedef Segment_2<R> Segment;
  typedef Triangle_2<R> Triangle;
  typedef Line_2<R> Line;
  typedef Direction_2<R> Direction;
  typedef Ray_2<R> Ray;

  typedef Distance_2<Triangulation_euclidean_traits_2<R> > Distance;

  Triangulation_euclidean_traits_2() {}
  Triangulation_euclidean_traits_2(const Triangulation_euclidean_traits_2 &) {}
  Triangulation_euclidean_traits_2 &operator=
      (const Triangulation_euclidean_traits_2 &)
  {return *this;}
 
    Comparison_result compare_x(const Point &p, const Point &q) const
    {
        return CGAL::compare_x(p, q);
    }


    Comparison_result compare_y(const Point &p, const Point &q) const
    {
        return CGAL::compare_y(p, q);
    }

  bool compare(const Point &p, const Point &q) const
  {
    return (compare_x(p, q)== EQUAL &&  
	    compare_y(p, q)== EQUAL);
  }

  Orientation orientation(const Point &p,
			  const Point &q,
			  const Point &r) const
    {
        return CGAL::orientation(p, q, r);
    }


  Oriented_side side_of_oriented_circle(const Point &p,
					const Point &q,
					const Point &r,
					const Point &s) const
    {
      return CGAL::side_of_oriented_circle(p, q, r, s);
    }

  Point circumcenter(const Point &p, const Point &q, const Point &r) const
    {
        return CGAL::circumcenter(p, q, r);
    }

     // Cette fonction devrait être accessible à tout CGAL.
  Line bisector(const Segment &s) const
  {
    Point p = midpoint (s.source(), s.target());
    Line l(s.source(), s.target());
    return l.perpendicular(p);
  }
};

CGAL_END_NAMESPACE 

#endif // CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
