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
// file          : demo/BooleanOperations/b-ops-2D-example-1.C
// source        : demo/BooleanOperations/b-ops-2D-example-1.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :                        Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ============================================================================

#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/bops_Polygon_2.h>
#include <vector>
#include <list>

typedef float TestNum;

#ifdef USE_CARTESIAN_COORDINATES
      typedef CGAL::Cartesian<TestNum>        R_type;
#else
      typedef CGAL::Homogeneous<TestNum>      R_type;
#endif

using std::list;
using std::cout;
using std::endl;

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
  /* instantiate Polygon A and B with containers */
  Polygon_2 A(container_A.begin(), container_A.end());
  Polygon_2 B(container_B.begin(), container_B.end());

  /* declaration of the result container */
  list<Object> result;

  /* performing intersection of A and B */
  CGAL::intersection(A, B, back_inserter(result));
  
  cout << "result size=" << result.size() << endl;

  /* possible results */
  Point_2   point;
  Segment_2 segment;
  Polygon_2 polygon;

  list<Object>::const_iterator it;
  for( it= result.begin(); it != result.end(); it++) {
    if( assign( polygon, *it) ) {
      cout << "PGN: " << polygon << endl;    /* polygon detected */
    }
    else if( assign( segment, *it) ) {
      cout << "SEG: " << segment << endl;    /* segment detected */
    }
    else if( assign( point, *it) )  {  
      cout << "PNT:" << point << endl;       /* point detected */
    }
    else {
      cout << "undefined object " << endl;   /* nothing detected */
    }
  }
  
  return result.size();
}



int main(void)
{
  Input_container container_A(6), container_B(4);

  container_A[0]= Point_2(2,4); /* description of polygon A */
  container_A[1]= Point_2(0,3);
  container_A[2]= Point_2(1,1);
  container_A[3]= Point_2(2,3);
  container_A[4]= Point_2(3,1);
  container_A[5]= Point_2(4,3);

  container_B[0]= Point_2(0,2); /* description of polygon B */
  container_B[1]= Point_2(0,0);
  container_B[2]= Point_2(5,0);
  container_B[3]= Point_2(5,2);

  example_intersection( container_A, container_B);
  return 0;
}
