#include <boost/test/unit_test.hpp>

#include "kernel/Line3.h"
#include "kernel/Point3.h"
#include "kernel/Vector3.h"

using kernel::Line3;
using kernel::Point3;
using kernel::Vector3;

BOOST_AUTO_TEST_SUITE(Line3Test)

BOOST_AUTO_TEST_CASE(testConstructor) {
    Point3 p(1.0, 1.0, 1.0);
    Point3 q(1.0, 2.0, 3.0);

    Line3 l(p, q);

    Vector3 expected(0.0, 1.0, 2.0);
    BOOST_CHECK(p == l.point());
    BOOST_CHECK(expected == l.direction());
}

BOOST_AUTO_TEST_CASE(testHasOn) {
    Point3 p(1.0, 1.0, 1.0);
    Point3 q(1.0, 2.0, 3.0);

    Line3 l(p, q);

    Point3 ison(1.0, 3.0, 5.0);
    Point3 isnot(4.0, 2.0, 1.0);
    BOOST_CHECK(l.hasOn(p));
    BOOST_CHECK(l.hasOn(q));
    BOOST_CHECK(l.hasOn(ison));
    BOOST_CHECK(!l.hasOn(isnot));
}

BOOST_AUTO_TEST_SUITE_END()

