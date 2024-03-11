#include <boost/test/unit_test.hpp>

#include "kernel/Point2.h"

using kernel::Point2;
using kernel::Vector2;

BOOST_AUTO_TEST_SUITE(Point2Test)

BOOST_AUTO_TEST_CASE(testConstructor) {
    Point2 p(1.0, 2.0);

    const double e = 0.001;
    BOOST_CHECK_CLOSE(p.getX(), 1.0, e);
    BOOST_CHECK_CLOSE(p.getY(), 2.0, e);
}

BOOST_AUTO_TEST_CASE(testOperator) {
    Point2 p(1.0, 2.0);
    Point2 q(2.0, 4.0);

    Vector2 v = (q-p);

    const double e = 0.001;
    BOOST_CHECK_CLOSE(v[0], 1.0, e);
    BOOST_CHECK_CLOSE(v[1], 2.0, e);
}

BOOST_AUTO_TEST_SUITE_END()

