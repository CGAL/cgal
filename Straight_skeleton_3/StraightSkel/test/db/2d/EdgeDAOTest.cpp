#include <boost/test/unit_test.hpp>

#include "db/2d/EdgeDAO.h"
#include "data/2d/ptrs.h"
#include "db/2d/ptrs.h"
#include "db/2d/DAOFactory.h"

using namespace db::_2d;

BOOST_AUTO_TEST_SUITE(EdgeDAOTest)

BOOST_AUTO_TEST_CASE(testAll) {
    PolygonSPtr polygon = Polygon::create();
    Point2SPtr p = KernelFactory::createPoint2(1.0, 2.0);
    Point2SPtr q = KernelFactory::createPoint2(3.0, 4.0);
    VertexSPtr src = Vertex::create(p);
    VertexSPtr dst = Vertex::create(q);
    EdgeSPtr edge = Edge::create(src, dst);
    EdgeDAOSPtr dao_edge = DAOFactory::getEdgeDAO();
    int returned = dao_edge->insert(edge);
    BOOST_CHECK_EQUAL(-1, returned);
    polygon->addEdge(edge);
    BOOST_CHECK(dao_edge->insert(edge));
    EdgeSPtr result = dao_edge->find(polygon->getID(), edge->getID());
    BOOST_CHECK_EQUAL(edge->getID(), result->getID());
    BOOST_CHECK_EQUAL(edge->getVertexSrc()->getID(),
            result->getVertexSrc()->getID());
    BOOST_CHECK_EQUAL(edge->getVertexDst()->getID(),
            result->getVertexDst()->getID());
    BOOST_CHECK(dao_edge->del(edge));
    VertexDAOSPtr vertex_dao = DAOFactory::getVertexDAO();
    vertex_dao->del(edge->getVertexSrc());
    vertex_dao->del(edge->getVertexDst());
}

BOOST_AUTO_TEST_SUITE_END()

