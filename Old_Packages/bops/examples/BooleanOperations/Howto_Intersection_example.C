//  -*- Mode: c++ -*-
// ============================================================================
// 
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-1.0 $
// release_date  : $CGAL_Date: 1998/09/12 $
//
// file          : examples/BooleanOperations/Howto_Intersection_example.C
// source        : examples/BooleanOperations/Howto_Intersection_example.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :                  Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ============================================================================

/*
  This is an example program for the usage of the boolean operation
  Intersection with cartesian coordiantes.

  The number type it uses is a non-exact one:    Quotient<int>
                               ---------
*/

#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
#include <CGAL/basic.h>

#include <list>
#include <CGAL/boolean_operations_2.h>

using std::list;

typedef CGAL::Quotient<long int>          TestNum;
typedef CGAL::Homogeneous<TestNum>        R_type;

typedef CGAL::Point_2<R_type>             Point;
typedef CGAL::Segment_2<R_type>           Segment;

typedef list<Point>                  Container;
typedef CGAL::Polygon_traits_2<R_type>    Polygon_traits;
typedef CGAL::Polygon_2< Polygon_traits, Container >  Polygon;

typedef std::vector<Point>                  Input_container;
using CGAL::Object;
using CGAL::assign;

int example_intersection(
  const Input_container& container_A,
  const Input_container& container_B
) {
  /* possible results */
  Point   point;
  Segment segment;
  Polygon polygon;

  /* instantiate Polygon A and B with containers */
  Polygon A(container_A.begin(), container_A.end());
  Polygon B(container_B.begin(), container_B.end());

  /* declaration of the result container */
  list<Object> result;

  /* performing intersection of A and B */
  CGAL::intersection(A, B, std::back_inserter(result));
  
  /* do something with the result */

  /* declaration of an iterator on the result container */
  list<Object>::const_iterator it;
  for( it= result.begin(); it != result.end(); it++) {
    if( assign( polygon, *it) ) {
        /* polygon detected */
    }
    else if( assign( segment, *it) ) {
        /* segment detected */
    }
    else if( assign( point, *it) )  {  
        /* point detected */
    }
    else {
        /* nothing detected */
    }
  }
  
  /* return size of result */
  return result.size();
}



int main(void)
{
  Input_container container_A(6), container_B(4);

  container_A[0]= Point(2,4);
  container_A[1]= Point(0,3);
  container_A[2]= Point(1,1);
  container_A[3]= Point(2,3);
  container_A[4]= Point(3,1);
  container_A[5]= Point(4,3);
  container_B[0]= Point(0,2);
  container_B[1]= Point(0,0);
  container_B[2]= Point(5,0);
  container_B[3]= Point(5,2);

  example_intersection( container_A, container_B);
  return 0;
}
