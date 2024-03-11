#include <boost/test/unit_test.hpp>

#include "data/3d/Edge.h"
#include "data/3d/ptrs.h"
#include "data/3d/KernelFactory.h"
#include "data/3d/Vertex.h"

using namespace data::_3d;

BOOST_AUTO_TEST_SUITE(EdgeTest)

BOOST_AUTO_TEST_CASE(testConstructor) {
    Point3SPtr p = KernelFactory::createPoint3(-1.0, -1.0, -1.0);
    Point3SPtr q = KernelFactory::createPoint3(1.0, 1.0, -1.0);
    VertexSPtr src = Vertex::create(p);
    VertexSPtr dst = Vertex::create(q);
    EdgeSPtr result = Edge::create(src, dst);
    BOOST_CHECK_EQUAL(src, result->getVertexSrc());
    BOOST_CHECK_EQUAL(dst, result->getVertexDst());
    BOOST_CHECK(src != dst);
    EdgeSPtr expected = src->edges().front().lock();
    BOOST_CHECK(result == expected);
    expected = dst->edges().front().lock();
    BOOST_CHECK(result == expected);
}

BOOST_AUTO_TEST_SUITE_END()
