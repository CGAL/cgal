#include <boost/test/unit_test.hpp>

#include "kernel/Vector2.h"

using kernel::Vector2;

BOOST_AUTO_TEST_SUITE(Vector2Test)

BOOST_AUTO_TEST_CASE(testConstructor) {
    Vector2 v(1.0, 2.0);

    const double e = 0.001;
    BOOST_CHECK_CLOSE(v[0], 1.0, e);
    BOOST_CHECK_CLOSE(v[1], 2.0, e);
}

BOOST_AUTO_TEST_CASE(testLength) {
    Vector2 v(1.0, 2.0);

    const double e = 0.001;
    BOOST_CHECK_CLOSE(v.length(), 2.236067, e);
}

BOOST_AUTO_TEST_CASE(testNormalize) {
    Vector2 v(1.0, 2.0);

    const double e = 0.001;
    BOOST_CHECK_CLOSE(v.normalize().length(), 1.0, e);
}

BOOST_AUTO_TEST_CASE(testIndex) {
    Vector2 v(1.0, 2.0);

    const double e = 0.001;
    BOOST_CHECK_CLOSE(v[0], 1.0, e);
    BOOST_CHECK_CLOSE(v[1], 2.0, e);
}

BOOST_AUTO_TEST_SUITE_END()

