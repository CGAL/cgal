// Copyright (c) 2001
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL__TEST_IO_H
#define CGAL__TEST_IO_H

#include <sstream>
#include <cassert>

template <class T>
void
_test_io_for(const T& t)
{
    std::stringstream ss;
    ss << t << std::endl;

    T u = t;
    ss >> u;
    assert(! ss.fail() );
    assert(u == t);
}

template <class R>
bool
_test_io(const R&)
{
    std::cout << "Testing IO for :" << std::endl;

    // 2D

    std::cout << "Point_2" << std::endl;
    typename R::Point_2 p2(2, 6, 2);
    typename R::Point_2 p22(2, -6, 2);
    typename R::Point_2 p222(-2, 6, 2);
    _test_io_for(p2);
    _test_io_for(p22);
    _test_io_for(p222);

    std::cout << "Weighted_point_2" << std::endl;
    typename R::Weighted_point_2 wp2(p2);
    typename R::Weighted_point_2 wp22(p22, 2);
    typename R::Weighted_point_2 wp222(p222, -2);
    _test_io_for(p2);
    _test_io_for(p22);
    _test_io_for(p222);

    std::cout << "Vector_2" << std::endl;
    typename R::Vector_2 v2(2, 6, 2);
    typename R::Vector_2 v22(2, -6, 2);
    typename R::Vector_2 v222(-2, 6, 2);
    _test_io_for(v2);
    _test_io_for(v22);
    _test_io_for(v222);

    std::cout << "Direction_2" << std::endl;
    typename R::Direction_2 d2(2, 6);
    typename R::Direction_2 d22(2, -6);
    typename R::Direction_2 d222(-2, 6);
    _test_io_for(d2);
    _test_io_for(d22);
    _test_io_for(d222);

    std::cout << "Segment_2" << std::endl;
    typename R::Segment_2 s2(p2, p22);
    typename R::Segment_2 s22(p2, p222);
    typename R::Segment_2 s222(p22, p22);
    _test_io_for(s2);
    _test_io_for(s22);
    _test_io_for(s222);

    std::cout << "Line_2" << std::endl;
    typename R::Line_2 l2(p2, p22);
    typename R::Line_2 l22(p2, p222);
    typename R::Line_2 l222(p22, p222);
    _test_io_for(l2);
    _test_io_for(l22);
    _test_io_for(l222);

    std::cout << "Ray_2" << std::endl;
    typename R::Ray_2 r2(p2, p22);
    typename R::Ray_2 r22(p2, p222);
    typename R::Ray_2 r222(p22, p222);
    _test_io_for(r2);
    _test_io_for(r22);
    _test_io_for(r222);

    std::cout << "Triangle_2" << std::endl;
    typename R::Triangle_2 t2(p2, p22, p222);
    _test_io_for(t2);

    std::cout << "Circle_2" << std::endl;
    typename R::Circle_2 c2(p2, p22, p222);
    _test_io_for(c2);

    std::cout << "Iso_rectangle_2" << std::endl;
    typename R::Iso_rectangle_2 i2(p2, p22);
    _test_io_for(i2);

    // 3D

    std::cout << "Point_3" << std::endl;
    typename R::Point_3 p3(2, 6, 2, 2);
    typename R::Point_3 p33(2, -6, 2, 2);
    typename R::Point_3 p333(-2, 6, 2, 2);
    typename R::Point_3 p3333(-2, 6, -2, 2);
    _test_io_for(p3);
    _test_io_for(p33);
    _test_io_for(p333);
    _test_io_for(p3333);

    std::cout << "Weighted_point_3" << std::endl;
    typename R::Weighted_point_3 wp3(p3);
    typename R::Weighted_point_3 wp33(p33, 3);
    typename R::Weighted_point_3 wp333(p333, -3);
    _test_io_for(p3);
    _test_io_for(p33);
    _test_io_for(p333);

    std::cout << "Vector_3" << std::endl;
    typename R::Vector_3 v3(2, 6, 2, 2);
    typename R::Vector_3 v33(2, -6, 2, 2);
    typename R::Vector_3 v333(-2, 6, 2, 2);
    _test_io_for(v3);
    _test_io_for(v33);
    _test_io_for(v333);

    std::cout << "Direction_3" << std::endl;
    typename R::Direction_3 d3(2, 6, 2);
    typename R::Direction_3 d33(2, -6, 2);
    typename R::Direction_3 d333(-2, 6, 2);
    _test_io_for(d3);
    _test_io_for(d33);
    _test_io_for(d333);

    std::cout << "Segment_3" << std::endl;
    typename R::Segment_3 s3(p3, p33);
    typename R::Segment_3 s33(p3, p333);
    typename R::Segment_3 s333(p33, p33);
    _test_io_for(s3);
    _test_io_for(s33);
    _test_io_for(s333);

    std::cout << "Line_3" << std::endl;
    typename R::Line_3 l3(p3, p33);
    typename R::Line_3 l33(p3, p333);
    typename R::Line_3 l333(p33, p333);
    _test_io_for(l3);
    _test_io_for(l33);
    _test_io_for(l333);

    std::cout << "Ray_3" << std::endl;
    typename R::Ray_3 r3(p3, p33);
    typename R::Ray_3 r33(p3, p333);
    typename R::Ray_3 r333(p33, p333);
    _test_io_for(r3);
    _test_io_for(r33);
    _test_io_for(r333);

    std::cout << "Triangle_3" << std::endl;
    typename R::Triangle_3 t3(p3, p33, p333);
    _test_io_for(t3);

    std::cout << "Sphere_3" << std::endl;
    typename R::FT squ_rad = 2;
    typename R::Sphere_3 S3(p3, squ_rad);
    _test_io_for(S3);

    std::cout << "Iso_cuboid_3" << std::endl;
    typename R::Iso_cuboid_3 i3(p3, p33);
    _test_io_for(i3);

    std::cout << "Tetrahedron_3" << std::endl;
    typename R::Tetrahedron_3 T3(p3, p33, p333, p3333);
    _test_io_for(T3);

    std::cout << "Plane_3" << std::endl;
    typename R::Plane_3 P3(p3, p33, p333);
    _test_io_for(P3);

    return true;
}

#endif // CGAL__TEST_IO_H
