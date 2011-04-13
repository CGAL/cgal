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
// file          : demo/Robustness/orientation_2.C
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
#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Ostream_iterator.h>
#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_real.h>
typedef leda_real exact_NT;
#else
#  include <CGAL/Quotient.h>
#  include <CGAL/MP_Float.h>
typedef CGAL::Quotient<CGAL::MP_Float>   exact_NT;
#endif
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Cartesian_converter.h>

#if defined(CGAL_USE_CGAL_WINDOW)
#define leda_window  CGAL::window
#define leda_string  std::string
#endif

#include <CGAL/orientation_test_statistics.h>

typedef CGAL::Cartesian<double>                CartesianDouble;
typedef CartesianDouble::Point_2               Point;
typedef CGAL::Creator_uniform_2<double,Point>  Pt_creator;
typedef std::vector<Point>                     Vector;
typedef CGAL::Cartesian<exact_NT>              CartesianLedaReal;

int
main(int argc, char** argv)
{
    int N = (argc > 1) ? CGAL_CLIB_STD::atoi(argv[1]) : 30;
    
    typedef leda_window  CGAL_Stream;
    CGAL_Stream W( 950, 520);
    CGAL_Stream W1( 400, 400);
    CGAL_Stream W2( 400, 400);
    CGAL::cgalize(W);
    CGAL::cgalize(W1);
    CGAL::cgalize(W2);
    
    W.init( 0, 950, 0);
    W1.init(-1.1, 1.1, -1.1);
    W2.init(-1.1, 1.1, -1.1);
    W.display();
    W1.display(W,50,50);
    W2.display(W,500,50);
    
    W.draw_ctext(475,495,"Orientation computation");
    W1 << CGAL::RED;
    W2 << CGAL::RED;
    
    Vector points1;
    Vector points2;
    typedef CGAL::Random_points_in_square_2<Point,Pt_creator>   P_in_square;
    typedef CGAL::Random_points_on_segment_2<Point,Pt_creator>  P_on_seg;
    P_in_square  pc1;
    Point p1 = Point(-0.8, 0.3);
    Point p2 = Point( 0.8,-0.4);
    P_on_seg     pc2(p1, p2);
    CGAL::copy_n( pc1, N, std::back_inserter(points1));
    CGAL::copy_n( pc2, N, std::back_inserter(points2));
    
    std::copy( points1.begin(), points1.end(),
               CGAL::Ostream_iterator< Point, CGAL_Stream>( W1));
    std::copy( points2.begin(), points2.end(),
               CGAL::Ostream_iterator< Point, CGAL_Stream>( W2));

    CGAL::Orientation* C = new CGAL::Orientation[N*N*N];
    std::vector< CartesianLedaReal::Point_2>  S;
    CGAL::Cartesian_converter<CartesianDouble, CartesianLedaReal> converter;
    leda_string s1;
    leda_string s2;


    std::transform( points1.begin(), points1.end(),
                    std::back_inserter( S), converter);
    fill_control_field( S.begin(), S.end(), C, CartesianLedaReal());
    orientation_statistics(points1.begin(), points1.end(),
                           C, s1, s2, CartesianDouble());
    W.draw_ctext( 250, 50, s1);
    W.draw_ctext( 250, 30, s2);


    std::transform( points2.begin(), points2.end(),
                    S.begin(), converter);
    fill_control_field( S.begin(), S.end(), C, CartesianLedaReal());
    orientation_statistics(points2.begin(), points2.end(),
                           C, s1, s2, CartesianDouble());
    W.draw_ctext( 700, 50, s1);
    W.draw_ctext( 700, 30, s2);

    W.read_mouse();
    delete[] C;

    return 0;
}
