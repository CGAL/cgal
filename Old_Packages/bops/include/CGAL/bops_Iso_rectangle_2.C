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
// release       : 
// release_date  : 1999, October 01
//
// file          : include/CGAL/bops_Iso_rectangle_2.C
// package       : bops (2.2)
// source        : include/CGAL/bops_Iso_rectangle_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL_BOPS_ISO_RECTANGLE_2_C
#define CGAL_BOPS_ISO_RECTANGLE_2_C

//#define CGAL__BOPS_USE_STANDARD_CODE
   // if CGAL__BOPS_USE_STANDARD_CODE is defined, then the
   // code for simple polygons will be used for computing 
   // BOPS on iso-rectangles and triangles.

#ifdef CGAL__BOPS_USE_STANDARD_CODE
#  include <CGAL/bops_simple_polygons_2.h>
#  include <CGAL/bops_Container_Polygon_2.h>
#endif

#include <CGAL/bops_Iso_rectangle_2.h>
#include <CGAL/bops_assertions.h>

CGAL_BEGIN_NAMESPACE

template < class R, class OutputIterator >
OutputIterator Union( const Iso_rectangle_2<R>& A,
                           const Iso_rectangle_2<R>& B,
		           OutputIterator result)
{
  CGAL_bops_precondition_msg( !A.is_degenerate(),
                              "Iso_rectangle_2<R> A is degenerated");
  CGAL_bops_precondition_msg( !B.is_degenerate(),
                              "Iso_rectangle_2<R> B is degenerated");

#ifdef CGAL__BOPS_USE_STANDARD_CODE
  Bops_default_I<R> default_traits;
  Bops_default_I<R>::Input_polygon TA(A), TB(B);
  return Union(
         TA.vertices_begin(), TA.vertices_end(),
         TB.vertices_begin(), TB.vertices_end(),
	 default_traits, result );
#else
  typedef Bops_default_I<R> I;
  I traits;
  typedef typename I::Point  Point;
  std::list<Point> point_list;
  if( do_intersect(A, B) ) {
      const Iso_rectangle_2<R>* r[2]= {&A,&B};
      short a, b;

      /*
                                    o---o 
         0          3               |b  | 
          o--------o             o--+---+--o
          |        |             |a |   |  |
          |        |             |  |   |  |
          o--------o             o--+---+--o
         1          2               |   | 
                                    o---o 

      */

 
      if( A.min().x() < B.min().x() ) { a= 0; b= 1; } else { a= 1; b= 0; }
      /* 0 */
      if( r[b]->max().y() > r[a]->max().y() ) {
        point_list.push_back( Point(r[b]->min().x(), r[b]->max().y()) );
        point_list.push_back( Point(r[b]->min().x(), r[a]->max().y()) );
      }
      point_list.push_back( Point(r[a]->min().x(), r[a]->max().y()) );
      /* 1 */
      point_list.push_back( Point(r[a]->min().x(), r[a]->min().y()) );
      if( r[b]->min().y() < r[a]->min().y() ) {
        point_list.push_back( Point(r[b]->min().x(), r[a]->min().y()) );
        point_list.push_back( Point(r[b]->min().x(), r[b]->min().y()) );
      }
 
      if( A.max().x() > B.max().x() ) { a= 0; b= 1; } else { a= 1; b= 0; }
      /* 2 */
      if( r[b]->min().y() < r[a]->min().y() ) {
        point_list.push_back( Point(r[b]->max().x(), r[b]->min().y()) );
        point_list.push_back( Point(r[b]->max().x(), r[a]->min().y()) );
        point_list.push_back( Point(r[a]->max().x(), r[a]->min().y()) );
      }
      point_list.push_back( Point(r[a]->max().x(), r[a]->min().y()) );
      /* 3 */
      point_list.push_back( Point(r[a]->max().x(), r[a]->max().y()) );
      if( r[b]->max().y() > r[a]->max().y() ) {
        point_list.push_back( Point(r[b]->max().x(), r[a]->max().y()) );
        point_list.push_back( Point(r[b]->max().x(), r[b]->max().y()) );
      }
 
      if( point_list.size() == 4 ) {
        std::list<Point>::const_iterator it= point_list.begin();
        *result++= traits.Make_object( I::Iso_rectangle(*++it++, *++it) );
      }
      else {
        *result++= traits.Make_object(
                     I::Output_polygon(point_list.begin(), point_list.end())
                   );
      }
  }
  else {
      *result++= traits.Make_object(A);
      *result++= traits.Make_object(B);
  }
  
  return result;
#endif // CGAL__BOPS_USE_STANDARD_CODE
}
 
 
template < class R, class OutputIterator >
OutputIterator difference( const Iso_rectangle_2<R>& A,
			        const Iso_rectangle_2<R>& B,
			        OutputIterator result)
{
  CGAL_bops_precondition_msg( !A.is_degenerate(),
                              "Iso_rectangle_2<R> A is degenerated");
  CGAL_bops_precondition_msg( !B.is_degenerate(),
                              "Iso_rectangle_2<R> B is degenerated");

#ifdef CGAL__BOPS_USE_STANDARD_CODE
  Bops_default_I<R> default_traits;
  Bops_default_I<R>::Input_polygon TA(A), TB(B);
  return difference(
         TA.vertices_begin(), TA.vertices_end(),
         TB.vertices_begin(), TB.vertices_end(),
	 default_traits, result );
#else // CGAL__BOPS_USE_STANDARD_CODE
  typedef Bops_default_I<R> I;
  typedef typename I::Point  Point;
  I traits;
  std::list<Point> point_list;

  if( do_intersect(A, B) ) {
   /* only for testing
    Point P[8];
    P[0]= Point(B.min().x(),A.max().y());
    P[1]= Point(B.max().x(),A.max().y());
    P[2]= Point(A.min().x(),B.max().y());
    P[3]= Point(A.min().x(),B.min().y());
    P[4]= Point(A.max().x(),B.max().y());
    P[5]= Point(A.max().x(),B.min().y());
    P[6]= Point(B.min().x(),A.min().y());
    P[7]= Point(B.max().x(),A.min().y());
   */
    
    bool A_min_x_in_B= A.min().x() > B.min().x() && A.min().x() < B.max().x();
    bool A_max_x_in_B= A.max().x() > B.min().x() && A.max().x() < B.max().x();
    bool A_min_y_in_B= A.min().y() > B.min().y() && A.min().y() < B.max().y();
    bool A_max_y_in_B= A.max().y() > B.min().y() && A.max().y() < B.max().y();
    bool B_min_x_in_A= B.min().x() > A.min().x() && B.min().x() < A.max().x();
    bool B_max_x_in_A= B.max().x() > A.min().x() && B.max().x() < A.max().x();
    bool B_min_y_in_A= B.min().y() > A.min().y() && B.min().y() < A.max().y();
    bool B_max_y_in_A= B.max().y() > A.min().y() && B.max().y() < A.max().y();

    bool P_defined[8];
    P_defined[0]= A_max_y_in_B && B_min_x_in_A; // P(d,a)
    P_defined[1]= A_max_y_in_B && B_max_x_in_A; // P(d,c)
    P_defined[2]= A_min_x_in_B && B_max_y_in_A; // P(a,d)
    P_defined[3]= A_min_x_in_B && B_min_y_in_A; // P(a,b)
    P_defined[4]= A_max_x_in_B && B_max_y_in_A; // P(c,d)
    P_defined[5]= A_max_x_in_B && B_min_y_in_A; // P(c,b)
    P_defined[6]= A_min_y_in_B && B_min_x_in_A; // P(b,a)
    P_defined[7]= A_min_y_in_B && B_max_x_in_A; // P(b,c)

    if( P_defined[0] ) {
      if( P_defined[6] ) {
        *result++= traits.Make_object(
                       I::Iso_rectangle(
                          A.min(),
                          Point(B.min().x(),A.max().y())
                       )
                   );
        if( P_defined[1] )  {                 /* case (a) + (c'') */
          //result= [p0, p1, P(b,a), P(d,a)] + [P(d,c),P(b,c),p2,p3];
          *result++= traits.Make_object(
                        I::Iso_rectangle(
                           I::Point(B.max().x(),A.min().y()),
                           A.max()
                        )
                     );
        }
        else {}                               /* case (a) */
          //result= [p0, p1, P(b,a), P(d,a)];
      }
      else {
        point_list.push_back(Point(A.min().x(), A.max().y())); // p0
        point_list.push_back(A.min());                         // p1
        point_list.push_back(Point(A.max().x(), A.min().y())); // p2
        if ( P_defined[5] ) {              /* case (b) */
          //result= [p0, p1, p2, P(c,b), r2(p1), P(d,a)]
          point_list.push_back(Point(A.max().x(), B.min().y())); // P(c,b)=P(5)
        }
        else if ( P_defined[1] ) {              /* case (c') + (c'') */
          //result= [p0, p1, p2, p3, P(d,c), r2(p2), r2(p1), P(d,a)];
          point_list.push_back(A.max());                         // p3
          point_list.push_back(Point(B.max().x(), A.max().y())); // P(d,c)=P(1)
          point_list.push_back(Point(B.max().x(), B.min().y())); // r2(p2)
        }
        point_list.push_back(B.min());                         // r2(p1)
        point_list.push_back(Point(B.min().x(), A.max().y())); // P(d,a)=P(0)
        *result++= traits.Make_object(
                   I::Output_polygon(point_list.begin(), point_list.end()));
     }
    }
    else if( P_defined[1] ) {
      if( P_defined[7] ) {                     /* case (a) + (c'') */
        //result= [ P(d,c), P(b,c), p2, p3];
        *result++= traits.Make_object(
                   I::Iso_rectangle(I::Point(B.max().x(), A.min().y()), A.max())
                   );
      }
      else if ( P_defined[3] ) {               /* case (b) */
        //result= [ P(d,c), r2(p2), P(a,b), p1, p2, p3];
        point_list.push_back(Point(B.max().x(), A.max().y())); // P(d,c)=P(1)
        point_list.push_back(Point(B.max().x(), B.min().y())); // r2(p2)
        point_list.push_back(Point(A.min().x(), B.min().y())); // P(a,b)=P(3)
        point_list.push_back(A.min());                         // p1
        point_list.push_back(Point(A.max().x(), A.min().y())); // p2
        point_list.push_back(A.max());                         // p3
        *result++= traits.Make_object(
                   I::Output_polygon(point_list.begin(), point_list.end()));
      }
    }
    else if( P_defined[2] ) {
      if( P_defined[4] ) {
        *result++= traits.Make_object(
                   I::Iso_rectangle(
                      I::Point(A.min().x(), B.max().y()),
                      A.max()) 
                   );
        if( P_defined[3] )           /* case (a) + (c'') */
          //result= [p0,P(a,d),P(c,d),p3] + [P(a,b),p1,p2,P(c,b)];
          *result++= traits.Make_object(
                     I::Iso_rectangle(A.min(), Point(A.max().x(), B.min().y()))
                     );
        //else                       /* case (a) */
           //result= [p0,P(a,d),P(c,d),p3];
      }
      else if ( P_defined[7] ) {     /* case (b) */
        //result= [p0,P(a,d),r2(p3),P(b,c),p2,p3];
        point_list.push_back(Point(A.min().x(), A.max().y())); // p0
        point_list.push_back(Point(A.min().x(), B.max().y())); // P(a,d)=P(2)
        point_list.push_back(B.max());                         // r2(p3)
        point_list.push_back(Point(B.max().x(), A.min().y())); // P(b,c)=P(7)
        point_list.push_back(Point(A.max().x(), A.min().y())); // p2
        point_list.push_back(A.max());                         // p3
        *result++= traits.Make_object(
                   I::Output_polygon(point_list.begin(), point_list.end()));
      }
      else if( P_defined[3] ) {   /* case (c') + (c'') */
        //result= [p0,P(a,d),r2(p3),r2(p2),P(a,b),p1,p2,p3];
        point_list.push_back(Point(A.min().x(), A.max().y())); // p0
        point_list.push_back(Point(A.min().x(), B.max().y())); // P(a,d)=P(2)
        point_list.push_back(B.max());                         // r2(p3)
        point_list.push_back(Point(B.max().x(), B.min().y())); // r2(p2)
        point_list.push_back(Point(A.min().x(), B.min().y())); // P(a,b)=P(3)
        point_list.push_back(A.min());                         // p1
        point_list.push_back(Point(A.max().x(), A.min().y())); // p2
        point_list.push_back(A.max());                         // p3
        *result++= traits.Make_object(
                   I::Output_polygon(point_list.begin(), point_list.end()));
      }

    }
    else if( P_defined[3] ) {
      if( P_defined[5] )           /* case (a) */
        //result= [P(a,b),p1,p2,P(c,b)];
        *result++= traits.Make_object(
                   I::Iso_rectangle(A.min(), Point(A.max().x(), B.min().y())) );

    }
    else if( P_defined[4] ) {
      point_list.push_back(Point(A.min().x(), A.max().y())); // p0
      point_list.push_back(A.min());                         // p1
      if( P_defined[5] )  {
        //result= [p0,p1,p2,P(c,b),r2(p1),r2(p0),P(c,d),p3);
        point_list.push_back(Point(A.max().x(), A.min().y())); // p2
        point_list.push_back(Point(A.max().x(), B.min().y())); // P(c,b)=P(5)
        point_list.push_back(B.min());                         // r2(p1)
      }
      else if( P_defined[6] ) {
        //result= [p0,p1, P(b,a),r2(p0),P(c,d),p3];
        point_list.push_back(Point(B.min().x(), A.min().y())); // P(b,a)=P(6)
      }
      point_list.push_back(Point(B.min().x(), B.max().y())); // r2(p0)
      point_list.push_back(Point(A.max().x(), B.max().y())); // P(c,d)=P(4)
      point_list.push_back(A.max());                         // p3
      *result++= traits.Make_object(
                 I::Output_polygon(point_list.begin(), point_list.end()));
    }
    else if( P_defined[6] ) {
      if( P_defined[7] ) {
        //result= [p0,p1,P(b,a),r2(p0),r2(p3),P(b,c),p2,p3];
        point_list.push_back(Point(A.min().x(), A.max().y())); // p0
        point_list.push_back(A.min());                         // p1
        point_list.push_back(Point(B.min().x(), A.min().y())); // P(b,a)=P(6)
        point_list.push_back(Point(B.min().x(), B.max().y())); // r2(p0)
        point_list.push_back(B.max());                         // r2(p3)
        point_list.push_back(Point(B.max().x(), A.min().y())); // P(b,c)=P(7)
        point_list.push_back(Point(A.max().x(), A.min().y())); // p2
        point_list.push_back(A.max());                         // p3
        *result++= traits.Make_object(
                   I::Output_polygon(point_list.begin(), point_list.end()));
      }
    }
    
  }
  else
    *result++= traits.Make_object(A);
  
  return result;
#endif // CGAL__BOPS_USE_STANDARD_CODE
} 
 
CGAL_END_NAMESPACE

#endif // CGAL_BOPS_ISO_RECTANGLE_2_C
