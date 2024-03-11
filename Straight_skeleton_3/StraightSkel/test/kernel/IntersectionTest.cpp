#include <boost/test/unit_test.hpp>

#include "kernel/intersection.h"
#include "kernel/Point2.h"
#include "kernel/Line2.h"
#include "kernel/Point3.h"
#include "kernel/Plane3.h"
#include "kernel/Line3.h"
#include "kernel/Vector3.h"

using kernel::intersection;
using kernel::Point2;
using kernel::Line2;
using kernel::Point3;
using kernel::Plane3;
using kernel::Line3;
using kernel::Vector3;

BOOST_AUTO_TEST_SUITE(IntersectionTest)

BOOST_AUTO_TEST_CASE(testIntersectionLines) {
    Point2 p1(1.0, 1.0);
    Point2 q1(3.0, 3.0);
    Line2 l1(p1, q1);

    Point2 p2(1.0, 3.0);
    Point2 q2(3.0, 1.0);
    Line2 l2(p2, q2);

    Point2 expected(2.0, 2.0);
    Point2* result = intersection(&l1, &l2);
    BOOST_CHECK(result);
    if (result) {
        BOOST_CHECK(expected == *result);
        delete result;
    }
    result = intersection(&l2, &l1);
    BOOST_CHECK(result);
    if (result) {
        BOOST_CHECK(expected == *result);
        delete result;
    }
}

BOOST_AUTO_TEST_CASE(testIntersectionLinesInf) {
    Point2 p1(2.0, 2.0);
    Point2 q1(2.0, 4.0);
    Line2 l1(p1, q1);

    Point2 p2(1.0, 3.0);
    Point2 q2(3.0, 3.0);
    Line2 l2(p2, q2);

    Point2 expected(2.0, 3.0);
    Point2* result = intersection(&l1, &l2);
    BOOST_CHECK(result);
    if (result) {
        BOOST_CHECK(expected == *result);
        delete result;
    }
    result = intersection(&l2, &l1);
    BOOST_CHECK(result);
    if (result) {
        BOOST_CHECK(expected == *result);
        delete result;
    }
}

BOOST_AUTO_TEST_CASE(testIntersectionLinesNone) {
    Point2 p1(1.0, 1.0);
    Point2 q1(2.0, 2.0);
    Line2 l1(p1, q1);

    Point2 p2(1.0, 2.0);
    Point2 q2(2.0, 3.0);
    Line2 l2(p2, q2);

    Point2* result = intersection(&l1, &l2);
    BOOST_CHECK(0 == result);
    result = intersection(&l2, &l1);
    BOOST_CHECK(0 == result);
}


BOOST_AUTO_TEST_CASE(testIntersection3Planes) {
    Point3 p11(1.0, 1.0, 3.0);
    Point3 p12(1.0, 2.0, 3.0);
    Point3 p13(2.0, 2.0, 3.0);
    Plane3 plane1(p11, p12, p13);

    Point3 p21(1.0, 2.0, 1.0);
    Point3 p22(1.0, 2.0, 2.0);
    Point3 p23(2.0, 2.0, 2.0);
    Plane3 plane2(p21, p22, p23);

    Point3 p31(1.0, 1.0, 1.0);
    Point3 p32(1.0, 1.0, 2.0);
    Point3 p33(1.0, 2.0, 2.0);
    Plane3 plane3(p31, p32, p33);

    Point3 expected(1.0, 2.0, 3.0);
    Point3* result = intersection(&plane1, &plane2, &plane3);
    BOOST_CHECK(result);
    if (result) {
        BOOST_CHECK(expected == *result);
        delete result;
    }
}

BOOST_AUTO_TEST_CASE(testIntersectionPlanes) {
    Point3 p11(2.0, 1.0, 1.0);
    Point3 p12(5.0, 1.0, 3.0);
    Point3 p13(2.0, 2.0, 1.0);
    Plane3 plane1(p11, p12, p13);

    Point3 p21(2.0, 1.0, 3.0);
    Point3 p22(5.0, 1.0, 1.0);
    Point3 p23(2.0, 2.0, 3.0);
    Plane3 plane2(p21, p22, p23);

    Point3 p(3.5, 0.0, 2.0);
    Vector3 dir(0.0, 1.0, 0.0);
    Line3 expected(p, dir);
    Line3* result = intersection(&plane1, &plane2);
    BOOST_CHECK(result);
    if (result) {
        BOOST_CHECK(expected == *result);
        delete result;
    }
}

BOOST_AUTO_TEST_CASE(testIntersectionPlanesNone) {
    Point3 p11(1.0, 1.0, 1.0);
    Point3 p12(2.0, 1.0, 2.0);
    Point3 p13(1.0, 2.0, 1.0);
    Plane3 plane1(p11, p12, p13);

    Point3 p21(1.0, 1.0, 2.0);
    Point3 p22(2.0, 1.0, 3.0);
    Point3 p23(1.0, 2.0, 2.0);
    Plane3 plane2(p21, p22, p23);

    Line3* result = intersection(&plane1, &plane2);
    BOOST_CHECK(0 == result);
    result = intersection(&plane2, &plane1);
    BOOST_CHECK(0 == result);
}

BOOST_AUTO_TEST_SUITE_END()
