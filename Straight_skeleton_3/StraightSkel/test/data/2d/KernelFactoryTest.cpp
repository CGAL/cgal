#include <boost/test/unit_test.hpp>

#include "data/2d/KernelFactory.h"
#include "data/2d/ptrs.h"

using namespace data::_2d;

BOOST_AUTO_TEST_SUITE(KernelFactoryTest)

BOOST_AUTO_TEST_CASE(testAll) {
    Point2SPtr p = KernelFactory::createPoint2(0.0, 0.0);
    Point2SPtr q = KernelFactory::createPoint2(2.0, 2.0);
    Vector2SPtr v = KernelFactory::createVector2(1.0, 1.0);
    Line2SPtr l1 = KernelFactory::createLine2(p, q);
    Line2SPtr l2 = KernelFactory::createLine2(p, v);
    BOOST_CHECK(*l1 == *l2);
}

BOOST_AUTO_TEST_SUITE_END()

