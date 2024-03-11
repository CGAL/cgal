#include <boost/test/unit_test.hpp>

#include "kernel/Segment2.h"

using kernel::Segment2;
using kernel::Point2;
using kernel::Line2;

BOOST_AUTO_TEST_SUITE(Segment2Test)

BOOST_AUTO_TEST_CASE(testConstructor) {
    Point2 p(1.0, 2.0);
    Point2 q(2.0, 2.0);

    Segment2 s(p,q);

    BOOST_CHECK_EQUAL(&p, &s.getP());
    BOOST_CHECK_EQUAL(&q, &s.getQ());
}

BOOST_AUTO_TEST_CASE(testLine) {
    Point2 p(1.0, 1.0);
    Point2 q(2.0, 1.0);

    Segment2 s(p,q);
    Line2 l = s.line();

    const double e = 0.001;
    BOOST_CHECK_CLOSE(l.getA(), 0.0, e);
    BOOST_CHECK_CLOSE(l.getB(), 1.0, e);
    BOOST_CHECK_CLOSE(l.getC(), -1.0, e);
}

BOOST_AUTO_TEST_SUITE_END()

