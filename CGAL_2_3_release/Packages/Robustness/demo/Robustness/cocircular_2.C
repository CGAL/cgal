// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : demo/Robustness/cocircular_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 
#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>

#include <cassert>
#include <vector>
#include <algorithm>

#include <CGAL/point_generators_2.h>
#include <CGAL/copy_n.h>
#include <CGAL/IO/Ostream_iterator.h>
#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_real.h>
typedef leda_real exact_NT;
#else
#  include <CGAL/MP_Float.h>
#  include <CGAL/Quotient.h>
typedef CGAL::Quotient<CGAL::MP_Float> exact_NT;
#endif

#include <CGAL/Interval_arithmetic.h>

#include <CGAL/Min_circle_2_traits_2.h>
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_ellipse_2_traits_2.h>
#include <CGAL/Min_ellipse_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Window_stream.h>

#if defined(CGAL_USE_CGAL_WINDOW)
#define leda_window  CGAL::window
#define leda_string  std::string
#define leda_blue    CGAL::blue
#define leda_pink    CGAL::pink
#define leda_grey1   CGAL::grey1
#endif

#include <CGAL/orientation_test_statistics.h>


typedef CGAL::Cartesian<double>                          CartesianDouble;
typedef CartesianDouble::Point_2                         Point;
typedef CGAL::Creator_uniform_2<double,Point>            Pt_creator;
typedef std::vector<Point>                               Vector;
typedef CGAL::Cartesian<exact_NT>                        CartesianLedaReal;
typedef CGAL::Min_circle_2_traits_2<CartesianDouble>     Traits;

typedef CGAL::Min_circle_2_traits_2<CartesianDouble>     Traits;
typedef CGAL::Min_circle_2<Traits>                       Min_circle;
typedef CGAL::Min_ellipse_2_traits_2<CartesianDouble>    ETraits;
typedef CGAL::Min_ellipse_2<ETraits>                     Min_ellipse;
typedef CGAL::Delaunay_triangulation_2<CartesianDouble>  DT_double;

int
main(int argc, char** argv)
{
    int N = (argc > 1) ? CGAL_CLIB_STD::atoi(argv[1]) : 30;
    
    CGAL::Window_stream W( 500, 550);
    CGAL::Window_stream W1( 400, 400);
    CGAL::cgalize(W);
    CGAL::cgalize(W1);
    
    W.init( 0, 500, 0);
    W1.init(-0.12, 0.12, -0.12);
    W.display();
    W1.display(W,50,50);
    
    Vector points1;
    typedef CGAL::Random_points_on_circle_2<Point,Pt_creator>   P_on_circle;
    P_on_circle  pc1( 0.1);
    CGAL::copy_n( pc1, N, std::back_inserter(points1));
    
    W1.set_fg_color(leda_pink);
    std::copy( points1.begin(), points1.end(),
               CGAL::Ostream_iterator< Point, CGAL::Window_stream>( W1));

    Min_circle  mc2( points1.begin(), points1.end(), true);
    W1 << mc2.circle();
    W.draw_ctext(250, 80, "min circle computed");

    Min_ellipse  me2( points1.begin(), points1.end(), true);
    W1.set_fg_color(leda_blue);
    W1 << me2.ellipse();
    W.draw_ctext(250, 60, "min ellipse computed");

    DT_double DTD;
    std::copy( points1.begin(), points1.end(), std::back_inserter( DTD ));
    W1.set_fg_color(leda_grey1);
    W1 << DTD;
    W1.set_fg_color(leda_pink);
    W1 << mc2.circle();
    W1.set_fg_color(leda_blue);
    W1 << me2.ellipse();
    W1.set_fg_color(leda_pink);
    std::copy( points1.begin(), points1.end(),
               CGAL::Ostream_iterator< Point, CGAL::Window_stream>( W1));
    W.draw_ctext(250, 40, "Delaunay triangulation computed");

    W.read_mouse();

    return 0;
}
