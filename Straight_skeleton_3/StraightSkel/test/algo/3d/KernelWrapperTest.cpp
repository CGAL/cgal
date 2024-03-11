#include <boost/test/unit_test.hpp>

#include "algo/3d/KernelWrapper.h"
#include "data/3d/ptrs.h"
#include "data/3d/KernelFactory.h"

BOOST_AUTO_TEST_SUITE(KernelWrapperTest)

using algo::_3d::KernelWrapper;
using data::_3d::KernelFactory;
using data::_3d::Point3SPtr;
using data::_3d::Plane3SPtr;

BOOST_AUTO_TEST_CASE(testOffsetPlane) {
    Point3SPtr p = KernelFactory::createPoint3(1.0, 1.0, 1.0);
    Point3SPtr q = KernelFactory::createPoint3(2.0, 1.0, 1.0);
    Point3SPtr r = KernelFactory::createPoint3(1.0, 2.0, 1.0);
    Plane3SPtr plane = KernelFactory::createPlane3(p, q, r);
    Plane3SPtr result = KernelWrapper::offsetPlane(plane, 1.0);
    Point3SPtr p_exp = KernelFactory::createPoint3(1.0, 1.0, 2.0);
    Point3SPtr q_exp = KernelFactory::createPoint3(2.0, 1.0, 2.0);
    Point3SPtr r_exp = KernelFactory::createPoint3(1.0, 2.0, 2.0);
    Plane3SPtr expected = KernelFactory::createPlane3(p_exp, q_exp, r_exp);
    BOOST_CHECK(*expected == *result);
}

BOOST_AUTO_TEST_SUITE_END()
