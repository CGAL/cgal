// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : examples/Optimisation/example_Min_ellipse_2.C
// package       : $CGAL_Package: Min_ellipse_2 $
// chapter       : Geometric Optimisation
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>, Bernd Gärtner
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// example progr.: 2D Smallest Enclosing Ellipse
// ============================================================================

// includes
#include <CGAL/Homogeneous.h>
#include <CGAL/Point_2.h>
#include <CGAL/Min_ellipse_2.h>
#include <CGAL/Min_ellipse_2_traits_2.h>
#include <CGAL/Gmpz.h>

// typedefs
typedef  CGAL::Gmpz                       NT;
typedef  CGAL::Homogeneous<NT>            K;
typedef  CGAL::Point_2<K>                 Point;
typedef  CGAL::Min_ellipse_2_traits_2<K>  Traits;
typedef  CGAL::Min_ellipse_2<Traits>      Min_ellipse;

// main
int
main( int, char**)
{
    int     n = 100;
    Point*  P = new Point[ n];

    for ( int i = 0; i < n; ++i)
	P[ i] = Point( (i%2 == 0 ? i : -i), 0);
    // (0,0), (-1,0), (2,0), (-3,0), ...

    Min_ellipse  me1( P, P+n, false);    // very slow
    Min_ellipse  me2( P, P+n, true);     // fast

    CGAL::set_pretty_mode( std::cout);
    std::cout << me2;

    delete[] P;

    return( 0);
}

// ===== EOF ==================================================================
