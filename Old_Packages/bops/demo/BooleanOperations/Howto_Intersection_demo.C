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
// file          : demo/BooleanOperations/Howto_Intersection_demo.C
// source        : demo/BooleanOperations/Howto_Intersection_demo.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :                        Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
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

//#include <iostream.h>

#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
#include <CGAL/basic.h>

#include <list>
#include <CGAL/boolean_operations_2.h>

using std::list;
using std::cout;
using std::endl;

typedef CGAL::Quotient<long int> TestNum;
typedef CGAL::Cartesian<TestNum>              R_type;
//typedef CGAL::Homogeneous<TestNum>            R_type;
typedef CGAL::Point_2<R_type>                 Point_2;
typedef CGAL::Segment_2<R_type>               Segment_2;

typedef list< Point_2 >                      Container;
typedef CGAL::Polygon_traits_2<R_type>        Polygon_traits_2;
typedef CGAL::Polygon_2< Polygon_traits_2, Container >  Polygon_2;


typedef std::vector<Point_2>                      Input_container;

int example_intersection(
  const Input_container& container_A,
  const Input_container& container_B
) {
  /* possible results */
  Point_2   point;
  Segment_2 segment;
  Polygon_2 polygon;

  /* instantiate Polygon A and B with containers */
  Polygon_2 A(container_A.begin(), container_A.end());
  Polygon_2 B(container_B.begin(), container_B.end());

  /* declaration of the result container */
  list<Object> result;

  /* performing intersection of A and B */
  CGAL::intersection(A, B, back_inserter(result));
  
  /* print out the result */
  cout << "size=" << result.size() << endl;

  /* declaration of an iterator on the result container */
  list<Object>::const_iterator it;
  for( it= result.begin(); it != result.end(); it++) {
    if( assign( polygon, *it) ) {
      cout << polygon << endl;               /* polygon detected */
    }
    else if( assign( segment, *it) ) {
      cout << segment << endl;               /* segment detected */
    }
    else if( assign( point, *it) )  {  
      cout << point << endl;                 /* point detected */
    }
    else {
      cout << "undefined object " << endl;   /* nothing detected */
    }
  }
  
  /* return size of result */
  return result.size();
}



int main(void)
{
  Input_container container_A(6), container_B(4);

  container_A[0]= Point_2(2,4);
  container_A[1]= Point_2(0,3);
  container_A[2]= Point_2(1,1);
  container_A[3]= Point_2(2,3);
  container_A[4]= Point_2(3,1);
  container_A[5]= Point_2(4,3);
  container_B[0]= Point_2(0,2);
  container_B[1]= Point_2(0,0);
  container_B[2]= Point_2(5,0);
  container_B[3]= Point_2(5,2);

  example_intersection( container_A, container_B);
  return 0;
}
