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
// release       : $CGAL_Revision: CGAL-wip $
// release_date  : $CGAL_Date$
//
// file          : examples/Optimisation/example_Min_circle_2.C
// revision      : $Revision$
// revision_date : $Date$
// package       : $CGAL_Package: Optimisation WIP $
// author(s)     : Sven Schönherr <sven@inf.fu-berlin.de>
//                 Bernd Gärtner
//
// coordinator   : ETH Zurich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// example progr.: 2D Smallest Enclosing Circle
// ============================================================================

// includes
#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Point_2.h>
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>
#include <iostream.h>

// typedefs
typedef  CGAL_Gmpz                      NT;
typedef  CGAL_Homogeneous<NT>           R;
typedef  CGAL_Point_2<R>                Point;
typedef  CGAL_Min_circle_2_traits_2<R>  Traits;
typedef  CGAL_Min_circle_2<Traits>      Min_circle;

// main
int
main( int, char**)
{
    int     n = 100;
    Point*  P = new Point[ n];

    for ( int i = 0; i < n; ++i)
	P[ i] = Point( (i%2 == 0 ? i : -i), 0);
    // (0,0), (-1,0), (2,0), (-3,0), ...

    Min_circle  mc1( P, P+n);           // very slow
    Min_circle  mc2( P, P+n, true);     // fast

    CGAL_set_pretty_mode( cout);
    cout << mc2;

    delete[] P;

    return( 0);
}

// ===== EOF ==================================================================
