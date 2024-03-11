#include <boost/test/unit_test.hpp>

#include "kernel/bisector.h"
#include "kernel/Point2.h"
#include "kernel/Line2.h"
#include "kernel/Point3.h"
#include "kernel/Plane3.h"

using kernel::bisector;
using kernel::Point2;
using kernel::Line2;
using kernel::Point3;
using kernel::Plane3;

BOOST_AUTO_TEST_SUITE(BisectorTest)

BOOST_AUTO_TEST_CASE(testBisectorLines) {
    Point2 p1(1.0, 1.0);
    Point2 q1(3.0, 1.0);
    Line2 l1(p1, q1);

    Point2 p2(1.0, 1.0);
    Point2 q2(3.0, 3.0);
    Line2 l2(p2, q2);

    Line2* result = bisector(&l1, &l2);
    BOOST_CHECK(result);
    if (result) {
        const double e = 0.001;
        BOOST_CHECK_CLOSE(result->getA(), -0.7071, e);
        BOOST_CHECK_CLOSE(result->getB(), 1.7071, e);
        BOOST_CHECK_CLOSE(result->getC(), -1.0, e);
        delete result;
    }
}

BOOST_AUTO_TEST_CASE(testBisectorLinesParallel) {
    Point2 p1(0.0, 1.0);
    Point2 q1(2.0, 2.0);
    Line2 l1(p1, q1);

    Point2 p2(0.0, 2.0);
    Point2 q2(4.0, 4.0);
    Line2 l2(p2, q2);

    Point2 p(0.0, 1.5);
    Point2 q(2.0, 2.5);
    Line2 expected(p, q);

    Line2* result = bisector(&l1, &l2);
    BOOST_CHECK(result);
    if (result) {
        const double e = 0.001;
        BOOST_CHECK_CLOSE(expected.getA(), result->getA(), e);
        BOOST_CHECK_CLOSE(expected.getB(), result->getB(), e);
        BOOST_CHECK_CLOSE(expected.getC(), result->getC(), e);
        delete result;
    }
}

BOOST_AUTO_TEST_CASE(testBisectorPlanes) {
    Point3 p11(0.0, 0.0, 1.0);
    Point3 p12(2.0, 0.0, 1.0);
    Point3 p13(0.0, 2.0, 1.0);
    Plane3 plane1(p11, p12, p13);

    Point3 p21(0.0, 0.0, 0.0);
    Point3 p22(2.0, 0.0, 2.0);
    Point3 p23(0.0, 2.0, 0.0);
    Plane3 plane2(p21, p22, p23);

    Plane3* result = bisector(&plane1, &plane2);
    BOOST_CHECK(result);
    if (result) {
        const double e = 0.001;
        BOOST_CHECK_CLOSE(result->getA(), -0.7071, e);
        BOOST_CHECK_CLOSE(result->getB(), 0.0, e);
        BOOST_CHECK_CLOSE(result->getC(), 1.7071, e);
        BOOST_CHECK_CLOSE(result->getD(), -1.0, e);
        delete result;
    }
}

BOOST_AUTO_TEST_SUITE_END()
