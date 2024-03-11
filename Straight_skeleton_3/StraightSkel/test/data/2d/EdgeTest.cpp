#include <boost/test/unit_test.hpp>

#include "data/2d/Edge.h"
#include "data/2d/ptrs.h"
#include "data/2d/KernelFactory.h"
#include "data/2d/Vertex.h"

using namespace data::_2d;

BOOST_AUTO_TEST_SUITE(EdgeTest)

BOOST_AUTO_TEST_CASE(testConstructor) {
    Point2SPtr p1 = KernelFactory::createPoint2(1.0, 2.0);
    Point2SPtr p2 = KernelFactory::createPoint2(2.0, 3.0);
    VertexSPtr v1 = Vertex::create(p1);
    VertexSPtr v2 = Vertex::create(p2);
    EdgeSPtr e = Edge::create(v1, v2);
    BOOST_CHECK_EQUAL(v1, e->getVertexSrc());
    BOOST_CHECK_EQUAL(v2, e->getVertexDst());
    BOOST_CHECK_EQUAL(v1->getEdgeOut(), e);
    BOOST_CHECK_EQUAL(v2->getEdgeIn(), e);
}

BOOST_AUTO_TEST_SUITE_END()

