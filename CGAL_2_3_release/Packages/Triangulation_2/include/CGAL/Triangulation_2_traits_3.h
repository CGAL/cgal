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
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Triangulation_2_traits_3.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>
//
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_2_TRAITS_3_H
#define CGAL_TRIANGULATION_2_TRAITS_3_H


#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_3.h>


#include <CGAL/triangulation_assertions.h>

CGAL_BEGIN_NAMESPACE 

template<class R>
class Compare_yz_3
{
public:
  typedef Point_3<R>     Point;

  Comparison_result operator() (Point p, Point q){
    Comparison_result r;
    r = CGAL::compare_y(p,q);
    if (r == EQUAL) r = CGAL::compare_z(p,q);
    return r;
   }
};
    

template <class R>
class Orientation_2_3 
{
public:
  typedef Point_3<R>     Point; 
  Orientation operator()(const Point& p,
			 const Point& q,
			 const Point& r)
    {
//       Orientation ori;
//       Point O(0.1111,0.1111,0.1111); 
//       Point A(1.1111,0,0);
//       Point B(0,1.1111,0);
//       Point C(0,0,1.1111);

//       Point P = ((ori = CGAL::orientation(p,q,r,O)) != ZERO) ? O:
//                 ((ori = CGAL::orientation(p,q,r,A)) != ZERO) ? A:
//                 ((ori = CGAL::orientation(p,q,r,B)) != ZERO) ? B:
//                 ((ori = CGAL::orientation(p,q,r,C)) != ZERO) ? C: C;
//       return CGAL::orientation(p,q,r,P);
    }
  
};

template <class R>
class Side_of_oriented_circle_2_3 
{
public:
  typedef Point_3<R>     Point; 
  CGAL::Oriented_side operator() (const Point &p, 
				  const Point &q,
				  const Point &r, 
				  const Point &s)
    {
//       //CGAL_triangulation_precondition( 
//       //              CGAL::orientation(p,q,r,s) == COPLANAR );
//       CGAL_triangulation_precondition( !CGAL::collinear(p,q,r) );

//       // test belongs to the circle if and only if it belongs to a
//       // sphere passing through pqr
//       Orientation ori;
//       Point O(0.1111,0.1111,0.1111); 
//       Point A(1.1111,0,0);
//       Point B(0,1.1111,0);
//       Point C(0,0,1.1111);

//       Point P = ((ori = CGAL::orientation(p,q,r,O)) != ZERO) ? O:
//                 ((ori = CGAL::orientation(p,q,r,A)) != ZERO) ? A:
//                 ((ori = CGAL::orientation(p,q,r,B)) != ZERO) ? B:
//                 ((ori = CGAL::orientation(p,q,r,C)) != ZERO) ? C: C;

//       return Oriented_side( ori *
// 	      CGAL::side_of_oriented_sphere(p, q, r, P, s));
      return CGAL::coplanar_side_of_bounded_circle(p,q,r,s);
    }
};



template < class R >
class Triangulation_2_traits_3 {
public:
  typedef R Rep;
  typedef Point_3<R>  Point_2;
  typedef Segment_3<R> Segment_2;
  typedef Triangle_3<R> Triangle_2;
 
  typedef typename R::Compare_x_3         Compare_x_2;
  typedef Compare_yz_3<R>                 Compare_y_2;
  typedef Orientation_2_3<R>              Orientation_2;
  typedef R::Coplanar_side_of_bounded_circle_3 
                                          Side_of_oriented_circle_2;
  
  typedef typename R::Construct_segment_3        Construct_segment_2;
  typedef typename R::Construct_triangle_3       Construct_triangle_2;

  // for compatibility with previous versions
  typedef Point_2      Point;
  typedef Segment_2    Segment;
  typedef Triangle_2   Triangle;

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
#endif // CGAL_TRIANGULATION_2_TRAITS_3_H
