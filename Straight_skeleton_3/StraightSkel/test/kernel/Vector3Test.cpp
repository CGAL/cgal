#include <boost/test/unit_test.hpp>

#include "kernel/Vector3.h"

using kernel::Vector3;

BOOST_AUTO_TEST_SUITE(Vector3Test)

BOOST_AUTO_TEST_CASE(testConstructor) {
    Vector3 v(1.0, 2.0, 3.0);

    const double e = 0.001;
    BOOST_CHECK_CLOSE(v[0], 1.0, e);
    BOOST_CHECK_CLOSE(v[1], 2.0, e);
    BOOST_CHECK_CLOSE(v[2], 3.0, e);
}

BOOST_AUTO_TEST_CASE(testLength) {
    Vector3 v(1.0, 2.0, 3.0);

    const double e = 0.001;
    BOOST_CHECK_CLOSE(v.length(), 3.741657, e);
}

BOOST_AUTO_TEST_CASE(testNormalize) {
    Vector3 v(1.0, 2.0, 3.0);

    const double e = 0.001;
    BOOST_CHECK_CLOSE(v.normalize().length(), 1.0, e);
}

BOOST_AUTO_TEST_CASE(testIndex) {
    Vector3 v(1.0, 2.0, 3.0);

    const double e = 0.001;
    BOOST_CHECK_CLOSE(v[0], 1.0, e);
    BOOST_CHECK_CLOSE(v[1], 2.0, e);
    BOOST_CHECK_CLOSE(v[2], 3.0, e);
}

BOOST_AUTO_TEST_CASE(testCross) {
    Vector3 u(1.0, 0.0, 0.0);
    Vector3 v(0.0, 1.0, 0.0);
    Vector3 expected(0.0, 0.0, 1.0);
    Vector3 result = u.cross(v);
    BOOST_CHECK(expected == result);
}

BOOST_AUTO_TEST_CASE(testAngle) {
    Vector3 u(1.0, 0.0, 0.0);
    Vector3 v(0.0, 1.0, 1.0);
    double result = u.angle(v);
    const double e = 0.001;
    BOOST_CHECK_CLOSE(1.5708, result, e);
    result = v.angle(u);
    BOOST_CHECK_CLOSE(1.5708, result, e);
}

BOOST_AUTO_TEST_SUITE_END()

