#include <boost/test/unit_test.hpp>

#include "algo/2d/KernelWrapper.h"
#include "data/2d/ptrs.h"
#include "data/2d/KernelFactory.h"

BOOST_AUTO_TEST_SUITE(KernelWrapperTest)

using algo::_2d::KernelWrapper;
using data::_2d::KernelFactory;
using data::_2d::Point2SPtr;
using data::_2d::Line2SPtr;

BOOST_AUTO_TEST_CASE(testIntersection) {
    Point2SPtr p = KernelFactory::createPoint2(1.0, 0.0);
    Point2SPtr q = KernelFactory::createPoint2(3.0, 2.0);
    Line2SPtr line1 = KernelFactory::createLine2(p, q);
    p = KernelFactory::createPoint2(1.0, 2.0);
    q = KernelFactory::createPoint2(3.0, 0.0);
    Line2SPtr line2 = KernelFactory::createLine2(p, q);
    Point2SPtr expected = KernelFactory::createPoint2(2.0, 1.0);
    Point2SPtr result = KernelWrapper::intersection(line1, line2);
    BOOST_CHECK(*expected == *result);
}

BOOST_AUTO_TEST_CASE(testDistance) {
    Point2SPtr p = KernelFactory::createPoint2(1.0, 1.0);
    Point2SPtr q = KernelFactory::createPoint2(3.0, 1.0);
    Line2SPtr l = KernelFactory::createLine2(p, q);
    q = KernelFactory::createPoint2(2.0, 2.0);
    double result = KernelWrapper::distance(l, q);
    BOOST_CHECK_EQUAL(1.0, result);
}

BOOST_AUTO_TEST_CASE(testOffsetLine) {
    Point2SPtr p = KernelFactory::createPoint2(1.0, 1.0);
    Point2SPtr q = KernelFactory::createPoint2(3.0, 1.0);
    Line2SPtr line = KernelFactory::createLine2(p, q);
    p = KernelFactory::createPoint2(1.0, 2.0);
    q = KernelFactory::createPoint2(3.0, 2.0);
    Line2SPtr expected = KernelFactory::createLine2(p, q);
    Line2SPtr result = KernelWrapper::offsetLine(line, 1.0);
    BOOST_CHECK(*expected == *result);

    p = KernelFactory::createPoint2(1.0, 1.0);
    q = KernelFactory::createPoint2(1.0, 3.0);
    line = KernelFactory::createLine2(p, q);
    p = KernelFactory::createPoint2(0.0, 1.0);
    q = KernelFactory::createPoint2(0.0, 3.0);
    expected = KernelFactory::createLine2(p, q);
    result = KernelWrapper::offsetLine(line, 1.0);
    BOOST_CHECK(*expected == *result);
}

BOOST_AUTO_TEST_SUITE_END()
