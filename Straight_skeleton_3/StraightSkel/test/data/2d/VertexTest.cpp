#include <boost/test/unit_test.hpp>

#include "data/2d/Vertex.h"
#include "data/2d/ptrs.h"
#include "data/2d/KernelFactory.h"

using namespace data::_2d;

BOOST_AUTO_TEST_SUITE(VertexTest)

BOOST_AUTO_TEST_CASE(testConstructor) {
    Point2SPtr p = KernelFactory::createPoint2(1.0, 2.0);
    VertexSPtr v = Vertex::create(p);
    BOOST_CHECK_EQUAL(p, v->getPoint());
    // compare raw pointers
    BOOST_CHECK_EQUAL(p.get(), v->getPoint().get());
}

BOOST_AUTO_TEST_SUITE_END()
