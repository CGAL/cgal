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
// file          : test/WindowStream/test_window_stream_xy_3.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/basic.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/Cartesian.h>
#include <cstddef>

#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/window_stream_xy_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Direction_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/IO/window_stream_xy_3.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Bbox_3.h>
// Check that multiple inclusion works.
#include <CGAL/IO/window_stream_xy_3.h>

typedef CGAL::Cartesian<double>                    REP;
typedef CGAL::Point_3<REP>                         Point3;
typedef CGAL::Vector_3<REP>                        Vector3;
typedef CGAL::Direction_3<REP>                     Direction3;
typedef CGAL::Line_3<REP>                          Line3;
typedef CGAL::Ray_3<REP>                           Ray3;
typedef CGAL::Segment_3<REP>                       Segment3;
typedef CGAL::Triangle_3<REP>                      Triangle3;
typedef CGAL::Tetrahedron_3<REP>                   Tetrahedron3;
typedef CGAL::Bbox_3                               Bbox3;

void test_window_stream_xy_3() {
    CGAL::Window_stream W(512, 512);
    CGAL::cgalize(W);
    W.init(-256.0, 255.0, -256.0);
    W.display(0,0);

    Point3         p1(10,10,10),
                   p2(10,-10,20),
                   p3(-10,10,30),
                   p4(-10,-10,23890);
    W << p1 << p2 << p3 << p4;

    // No output available for 2d vectors at the moment.
    // Vector_3       v(50,10,-300);
    // W << CGAL::GREEN << v;

    // No output available for 2d directions at the moment.
    // Direction_3    d(40,-20,5);
    // W << CGAL::GREEN << d;

    p1 = Point3(-200,100,32);
    p2 = Point3( 200,120,32);
    Line3         l( p1, p2);
    W << CGAL::BLUE << l;
    p1 = Point3(-200,140,32);
    p2 = Point3( 200,160,32);
    Ray3          ray( p1, p2);
    W << ray;
    p1 = Point3(-200,180,32);
    p2 = Point3( 200,200,32);
    Segment3      segment( p1, p2);
    W << segment;
    p1 = Point3(-125, -50,32);
    p2 = Point3( -50,-200,32);
    p3 = Point3(-200,-200,32);
    Triangle3     triangle( p1, p2, p3);
    W << CGAL::RED << triangle;
    p1 = Point3( 125,- 50,32);
    p2 = Point3(  50,-200,32);
    p3 = Point3( 200,-200,32);
    p4 = Point3( 125,-150,32);
    Tetrahedron3  tetra( p1, p2, p3, p4);
    W << tetra;
    Bbox3         bbox(-240,-240,-240, 240, 240, 240);
    W << CGAL::BLACK << bbox;

    W >> p1;
    W.clear();

    // Test stream input operators.
    W << CGAL::BLACK;
    W >> p1;
    // No input available for 2d vectors and directions at the moment.
    // W >> v >> d;
    W << CGAL::BLUE;
    W >> l >> ray >> segment;
    W << CGAL::RED;
    W >> triangle >> tetra;
    W << CGAL::BLACK;
    W >> bbox;

    W >> p1;
}

int
main( int argc, char *argv[])
{
  if (argc >= 2) { return 0; }
  test_window_stream_xy_3();

  return 0;
}
#else
int main() { return 0; }
#endif // CGAL_USE_LEDA
