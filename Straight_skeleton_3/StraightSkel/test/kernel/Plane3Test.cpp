#include <boost/test/unit_test.hpp>

#include "kernel/Plane3.h"

using kernel::Plane3;
using kernel::Point3;

BOOST_AUTO_TEST_SUITE(Plane3Test)

BOOST_AUTO_TEST_CASE(testConstructor) {
    Point3 p1(1.0, 1.0, 1.0);
    Point3 p2(2.0, 1.0, 1.0);
    Point3 p3(1.0, 2.0, 1.0);
    Plane3 h(p1, p2, p3);

    // 0 = a*x + b*y + c*z + d
    // z = - a/c * x - b/c * y - d/c
    // all points on same height:
    // z = -d/c

    const double e = 0.001;
    BOOST_CHECK_CLOSE(h.getA(), 0.0, e);
    BOOST_CHECK_CLOSE(h.getB(), 0.0, e);
    BOOST_CHECK_CLOSE(h.getC(), 1.0, e);
    BOOST_CHECK_CLOSE(h.getD(), -1.0, e);

    p2 = Point3(2.0, 1.0, 2.0);
    p3 = Point3(1.0, 2.0, 1.0);
    h = Plane3(p1, p2, p3);

    BOOST_CHECK_CLOSE(h.getA(), -1.0, e);
    BOOST_CHECK_CLOSE(h.getB(), 0.0, e);
    BOOST_CHECK_CLOSE(h.getC(), 1.0, e);
    BOOST_CHECK_CLOSE(h.getD(), 0.0, e);

    // 3 collinear points
    p2 = Point3(1.0, 1.0, 2.0);
    p3 = Point3(1.0, 1.0, 3.0);

    h = Plane3(p1, p2, p3);
    BOOST_CHECK_CLOSE(h.getA(), 0.0, e);
    BOOST_CHECK_CLOSE(h.getB(), 0.0, e);
    BOOST_CHECK_CLOSE(h.getC(), 0.0, e);
    BOOST_CHECK_CLOSE(h.getD(), 0.0, e);
}

BOOST_AUTO_TEST_CASE(testOpposite) {
    Point3 p1(1.0, 1.0, 1.0);
    Point3 p2(2.0, 2.0, 2.0);
    Point3 p3(1.0, 1.0, 3.0);
    Plane3 plane(p1, p2, p3);
    Plane3 expected(p1, p3, p2);
    Plane3 result = plane.opposite();
    BOOST_CHECK(expected == result);
}

BOOST_AUTO_TEST_CASE(testSide) {
    Point3 p1(1.0, 1.0, 1.0);
    Point3 p2(2.0, 1.0, 1.0);
    Point3 p3(1.0, 2.0, 1.0);
    Plane3 plane(p1, p2, p3);
    Point3 point(2.0, 2.0, 2.0);
    int result = plane.side(point);
    BOOST_CHECK_EQUAL(1, result);
}

BOOST_AUTO_TEST_SUITE_END()

