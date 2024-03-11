#include <boost/test/unit_test.hpp>

#include "db/2d/VertexDAO.h"
#include "data/2d/ptrs.h"
#include "db/2d/ptrs.h"
#include "db/2d/DAOFactory.h"

using namespace db::_2d;

BOOST_AUTO_TEST_SUITE(VertexDAOTest)

BOOST_AUTO_TEST_CASE(testAll) {
    PolygonSPtr polygon = Polygon::create();
    Point2SPtr point = KernelFactory::createPoint2(1.0, 2.0);
    VertexSPtr vertex = Vertex::create(point);
    VertexDAOSPtr dao_vertex = DAOFactory::getVertexDAO();
    int returned = dao_vertex->insert(vertex);
    BOOST_CHECK_EQUAL(-1, returned);
    polygon->addVertex(vertex);
    BOOST_CHECK(dao_vertex->insert(vertex));
    VertexSPtr result = dao_vertex->find(polygon->getID(), vertex->getID());
    BOOST_CHECK_EQUAL(vertex->getID(), result->getID());
    const double e = 0.001;
    BOOST_CHECK_CLOSE(vertex->getX(), result->getX(), e);
    BOOST_CHECK_CLOSE(vertex->getY(), result->getY(), e);
    BOOST_CHECK(dao_vertex->del(vertex));
    VertexSPtr nothing = dao_vertex->find(polygon->getID(), result->getID());
    BOOST_CHECK_EQUAL(nothing, VertexSPtr());
}

BOOST_AUTO_TEST_SUITE_END()

