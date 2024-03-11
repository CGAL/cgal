#include <boost/test/unit_test.hpp>

#include "kernel/Point3.h"

using kernel::Point3;
using kernel::Vector3;

BOOST_AUTO_TEST_SUITE(Point3Test)

BOOST_AUTO_TEST_CASE(testConstructor) {
    Point3 p(1.0, 2.0, 3.0);

    const double e = 0.001;
    BOOST_CHECK_CLOSE(p.getX(), 1.0, e);
    BOOST_CHECK_CLOSE(p.getY(), 2.0, e);
    BOOST_CHECK_CLOSE(p.getZ(), 3.0, e);
}

BOOST_AUTO_TEST_CASE(testOperator) {
    Point3 p(1.0, 2.0, 3.0);
    Point3 q(2.0, 4.0, 6.0);

    Vector3 v = (q-p);

    const double e = 0.001;
    BOOST_CHECK_CLOSE(v[0], 1.0, e);
    BOOST_CHECK_CLOSE(v[1], 2.0, e);
    BOOST_CHECK_CLOSE(v[2], 3.0, e);
}

BOOST_AUTO_TEST_SUITE_END()

