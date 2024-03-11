#include <boost/test/unit_test.hpp>

#include "algo/3d/PolyhedronBuilder.h"
#include "data/3d/ptrs.h"
#include "data/3d/KernelFactory.h"

BOOST_AUTO_TEST_SUITE(PolyhedronBuilderTest)

using algo::_3d::PolyhedronBuilder;
using data::_3d::KernelFactory;
using data::_3d::Point3SPtr;
using data::_3d::PolyhedronSPtr;

BOOST_AUTO_TEST_CASE(testMakeTetrahedron) {
    Point3SPtr p1 = KernelFactory::createPoint3(-10.0, -10.0, -10.0);
    Point3SPtr p2 = KernelFactory::createPoint3(10.0, 10.0, -10.0);
    Point3SPtr p3 = KernelFactory::createPoint3(10.0, -10.0, 10.0);
    Point3SPtr p4 = KernelFactory::createPoint3(-10.0, 10.0, 10.0);
    PolyhedronSPtr polyhedron =
            PolyhedronBuilder::makeTetrahedron(p1, p2, p3, p4);
    BOOST_CHECK_EQUAL(4, polyhedron->vertices().size());
    BOOST_CHECK_EQUAL(6, polyhedron->edges().size());
    BOOST_CHECK_EQUAL(4, polyhedron->facets().size());
    BOOST_CHECK(polyhedron->isConsistent());
    polyhedron = PolyhedronBuilder::makeTetrahedron(p1, p3, p2, p4);
    BOOST_CHECK(polyhedron->isConsistent());
    polyhedron = PolyhedronBuilder::makeTetrahedron(p4, p3, p2, p1);
    BOOST_CHECK(polyhedron->isConsistent());
}

BOOST_AUTO_TEST_SUITE_END()
