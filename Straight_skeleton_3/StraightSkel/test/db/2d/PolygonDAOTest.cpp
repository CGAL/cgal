#include <boost/test/unit_test.hpp>

#include "db/2d/PolygonDAO.h"
#include "data/2d/ptrs.h"
#include "db/2d/ptrs.h"
#include "db/2d/DAOFactory.h"

using namespace db::_2d;

BOOST_AUTO_TEST_SUITE(PolygonDAOTest)

BOOST_AUTO_TEST_CASE(testAll) {
    unsigned int i = 0;
    const unsigned int num_points = 3;
    Point2SPtr p[num_points];
    p[0] = KernelFactory::createPoint2(1.0, 1.0);
    p[1] = KernelFactory::createPoint2(3.0, 1.0);
    p[2] = KernelFactory::createPoint2(2.0, 3.0);
    VertexSPtr v[num_points];
    for (i = 0; i < num_points; i++) {
        v[i] = Vertex::create(p[i]);
    }
    EdgeSPtr e[num_points];
    for (i = 0; i < num_points; i++) {
        e[i] = Edge::create(v[(i+1)%num_points], v[(i+2)%num_points]);
    }
    PolygonSPtr polygon = Polygon::create();
    for (i = 0; i < num_points; i++) {
        polygon->addVertex(v[i]);
    }
    for (i = 0; i < num_points; i++) {
        polygon->addEdge(e[i]);
    }

    PolygonDAOSPtr polygon_dao = DAOFactory::getPolygonDAO();
    polygon_dao->insert(polygon);
    PolygonSPtr result = polygon_dao->find(polygon->getID());
    BOOST_CHECK(result->isConsistent());
    BOOST_CHECK_EQUAL(polygon->getID(), result->getID());
    BOOST_CHECK_EQUAL(polygon->edges().size(), result->edges().size());
    BOOST_CHECK_EQUAL(polygon->vertices().size(), result->vertices().size());
    polygon_dao->del(polygon);
}

BOOST_AUTO_TEST_SUITE_END()

