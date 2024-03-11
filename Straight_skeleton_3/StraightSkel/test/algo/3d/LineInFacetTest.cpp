#include <boost/test/unit_test.hpp>

#include "algo/3d/LineInFacet.h"
#include "data/3d/ptrs.h"
#include "data/3d/KernelFactory.h"
#include "data/3d/Vertex.h"
#include "data/3d/Facet.h"

BOOST_AUTO_TEST_SUITE(LineInFacetTest)

using namespace data::_3d;
using algo::_3d::IsLineInFacet;

BOOST_AUTO_TEST_CASE(testIsLineInFacet) {
    unsigned int num_vertices = 4;
    Point3SPtr points[num_vertices];
    points[0] = KernelFactory::createPoint3(-3.0, -3.0, 1.0);
    points[1] = KernelFactory::createPoint3(3.0, -3.0, 1.0);
    points[2] = KernelFactory::createPoint3(3.0, 3.0, 1.0);
    points[3] = KernelFactory::createPoint3(-3.0, 3.0, 1.0);
    VertexSPtr vertices[num_vertices];
    for (unsigned int i = 0; i < num_vertices; i++) {
        vertices[i] = Vertex::create(points[i]);
    }
    FacetSPtr facet = Facet::create(num_vertices, vertices);
    Point3SPtr p = KernelFactory::createPoint3(0.0, 0.0, 0.0);
    Point3SPtr q = KernelFactory::createPoint3(1.0, 0.5, 2.0);
    Line3SPtr line_inside = KernelFactory::createLine3(p, q);
    bool result = IsLineInFacet(facet, line_inside);
    BOOST_CHECK_EQUAL(true, result);
    q = KernelFactory::createPoint3(1.0, 0.0, 0.1);
    Line3SPtr line_outside = KernelFactory::createLine3(p, q);
    result = IsLineInFacet(facet, line_outside);
    BOOST_CHECK_EQUAL(false, result);
}

BOOST_AUTO_TEST_SUITE_END()
