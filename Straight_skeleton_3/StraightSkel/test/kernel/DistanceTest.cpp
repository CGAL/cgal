#include <boost/test/unit_test.hpp>

#include "kernel/distance.h"
#include "kernel/Point2.h"
#include "kernel/Line2.h"
#include "kernel/Point3.h"
#include "kernel/Plane3.h"

using kernel::distance;
using kernel::Point2;
using kernel::Line2;
using kernel::Point3;
using kernel::Plane3;

BOOST_AUTO_TEST_SUITE(DistanceTest)

BOOST_AUTO_TEST_CASE(testDistanceLinePoint) {
    Point2 p(1.0, 1.0);
    Point2 q(3.0, 1.0);
    Line2 line(p, q);

    Point2 point(2.0, 2.0);

    double result = distance(&line, &point);
    BOOST_CHECK_EQUAL(1.0, result);
}

BOOST_AUTO_TEST_CASE(testDistancePlanePoint) {
    Point3 p1(1.0, 1.0, 1.0);
    Point3 p2(3.0, 1.0, 1.0);
    Point3 p3(1.0, 3.0, 1.0);
    Plane3 plane(p1, p2, p3);

    Point3 point(2.0, 2.0, 2.0);

    double result = distance(&plane, &point);
    BOOST_CHECK_EQUAL(1.0, result);
}

BOOST_AUTO_TEST_SUITE_END()
