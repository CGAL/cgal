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
// file          : demo/Min_ellipse_2/demo_Min_ellipse_2.C
// package       : $CGAL_Package: Min_ellipse_2 $
// chapter       : Geometric Optimisation
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// demo program  : 2D Smallest Enclosing Ellipse
// ============================================================================

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// NOTE: In this release correct results are only guaranteed if exact
// arithmetic is used, so this demo (using inexact floating-point
// arithmetic) is only intended to illustrate the techniques. However,
// the program will terminate, but the computed ellipse may neither
// contain all points nor be the smallest one.
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// includes
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Min_ellipse_2.h>
#include <CGAL/Min_ellipse_2_traits_2.h>
#include <CGAL/IO/Window_stream.h>

// typedefs
typedef  CGAL::Cartesian< double >          K;
typedef  CGAL::Point_2< K >                 Point;
typedef  CGAL::Min_ellipse_2_traits_2< K >  Traits;
typedef  CGAL::Min_ellipse_2< Traits >      Min_ellipse;

// main
int
main( int, char**)
{
    using namespace std;

    cerr << "  left button: insert point" << endl;
    cerr << "middle button: clear points" << endl;
    cerr << " right button: exit"         << endl;

    // create empty min_ellipse;
    Min_ellipse  me;

    // open window
    CGAL::Window_stream ws( "CGAL Demo: Smallest Enclosing Ellipse in 2D");
    ws.set_icon_label("CGAL");
    ws.set_node_width( 5);
    ws.init( -100.0, 100.0, -100.0);
    ws.display();

    // main loop
    double  x, y;
    int     button, i;
    do {
	// get mouse click
	button = ws.read_mouse( x, y);

	switch ( button) {

	  case MOUSE_BUTTON( 1):                        // left button
	    ws << CGAL::WHITE << me.ellipse();
	    me.insert( Point( x, y));
	    ws << CGAL::BLACK << me;
	    ws << CGAL::BLUE  << me.ellipse();
	    ws << CGAL::RED;
	    for ( i = 0; i < me.number_of_support_points(); ++i)
		ws << me.support_point( i);
	    break;

	  case MOUSE_BUTTON( 2):                        // middle button
	    ws.clear();
	    me.clear();
	    break; } }

    while ( button != MOUSE_BUTTON( 3));                // right button

    return( 0);
}

// ===== EOF ==================================================================
