#include <boost/test/unit_test.hpp>

#include "kernel/Line2.h"

using kernel::Line2;
using kernel::Point2;
using kernel::Vector2;

BOOST_AUTO_TEST_SUITE(Line2Test)

BOOST_AUTO_TEST_CASE(testConstructor) {
    Point2 p(1.0, 1.0);
    Point2 q(2.0, 1.0);

    Line2 l(p, q);

    const double e = 0.001;
    BOOST_CHECK_CLOSE(l.getA(), 0.0, e);
    BOOST_CHECK_CLOSE(l.getB(), 1.0, e);
    BOOST_CHECK_CLOSE(l.getC(), -1.0, e);
}

BOOST_AUTO_TEST_CASE(testDirection) {
    Point2 p(1.0, 1.0);
    Point2 q(2.0, 3.0);

    Line2 l(p, q);
    Vector2 direction = l.direction();

    const double e = 0.001;
    BOOST_CHECK_CLOSE(direction[0], 1.0, e);
    BOOST_CHECK_CLOSE(direction[1], 2.0, e);
}

BOOST_AUTO_TEST_CASE(testOpposite) {
    Point2 p(1.0, 1.0);
    Point2 q(2.0, 2.0);
    Line2 l(p, q);
    Line2 expected(q, p);
    Line2 result = l.opposite();
    BOOST_CHECK(expected == result);
}

BOOST_AUTO_TEST_SUITE_END()
