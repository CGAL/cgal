// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.0-I-12 $
// release_date  : $CGAL_Date: 1999/04/28 $
//
// file          : include/CGAL/Triangulation_euclidean_traits_xz_3.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================


#ifndef CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_XZ_3_H
#define CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_XZ_3_H

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/triangulation_assertions.h>

#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/predicates/kernel_ftC2.h>

CGAL_BEGIN_NAMESPACE

template <class R>
class Orientation_xz_3 
{
public:
  typedef Point_3<R>     Point; 
  typename R::FT x(const Point &p) const { return p.x(); }
  typename R::FT y(const Point &p) const { return p.z(); }

  CGAL::Orientation operator()(const Point& p,
			       const Point& q,
			       const Point& r)
    {
      return orientationC2(x(p), y(p), x(q), y(q), x(r), y(r));
    }
};

template <class R>
class Side_of_oriented_circle_xz_3 
{
public:
  typedef Point_3<R>     Point; 
  typename R::FT x(const Point &p) const { return p.x(); }
  typename R::FT y(const Point &p) const { return p.z(); }

  CGAL::Oriented_side operator() (const Point &p, 
				  const Point &q,
				  const Point &r, 
				  const Point &s) 
    {
      return side_of_oriented_circleC2(x(p), y(p),
				       x(q), y(q),
				       x(r), y(r),
				       x(s), y(s));
    }
};

template < class R >
class Triangulation_euclidean_traits_xz_3 {
public:
  typedef Triangulation_euclidean_traits_xz_3<R> Traits;
  typedef R Rp;
  typedef Point_3<R>    Point_2;
  typedef Segment_3<R>  Segment_2;
  typedef Triangle_3<R> Triangle_2;

  typedef typename Rp::Compare_x_3          Compare_x_2;
  typedef typename Rp::Compare_z_3          Compare_y_2;
  typedef Orientation_xz_3<Rp>              Orientation_2;
  typedef Side_of_oriented_circle_xz_3<Rp>  Side_of_oriented_circle_2;
  typedef typename Rp::Construct_segment_3   Construct_segment_2;
  typedef typename Rp::Construct_triangle_3  Construct_triangle_2;

  // for compatibility with previous versions
  typedef Point_2      Point;
  typedef Segment_2    Segment;
  typedef Triangle_2   Triangle;
  
  Triangulation_euclidean_traits_xz_3(){}
  Triangulation_euclidean_traits_xz_3(
	      const Triangulation_euclidean_traits_xz_3& ){}
  Triangulation_euclidean_traits_xz_3 &operator=(
	 const Triangulation_euclidean_traits_xz_3& ){return *this;}

  typename Rp::FT x(const Point_2 &p) const { return p.x(); }
  typename Rp::FT y(const Point_2 &p) const { return p.z(); }

  Compare_x_2
  compare_x_2_object() const
    { return Compare_x_2();}

  Compare_y_2
  compare_y_2_object() const
    { return Compare_y_2();}

  Orientation_2
  orientation_2_object() const
    { return Orientation_2();}

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
    {return Side_of_oriented_circle_2();}

  Construct_segment_2  construct_segment_2_object() const
    {return Construct_segment_2();}

  Construct_triangle_2  construct_triangle_2_object() const
    {return Construct_triangle_2();}
};

CGAL_END_NAMESPACE


#endif // CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_XZ_3_H
