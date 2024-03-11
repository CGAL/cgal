#include <boost/test/unit_test.hpp>

#include "db/3d/PolyhedronDAO.h"
#include "data/3d/ptrs.h"
#include "db/3d/ptrs.h"
#include "db/3d/DAOFactory.h"

using namespace db::_3d;

BOOST_AUTO_TEST_SUITE(PolyhedronDAOTest)

BOOST_AUTO_TEST_CASE(testAll) {
    const unsigned int num_vertices = 4;
    const unsigned int num_edges = 6;
    const unsigned int num_facets = 4;
    Point3SPtr points[num_vertices];
    points[0] = KernelFactory::createPoint3(-1.0, -1.0, -1.0);
    points[1] = KernelFactory::createPoint3(1.0, 1.0, -1.0);
    points[2] = KernelFactory::createPoint3(1.0, -1.0, 1.0);
    points[3] = KernelFactory::createPoint3(-1.0, 1.0, 1.0);
    VertexSPtr vertices[num_vertices];
    for (unsigned int i = 0; i < num_vertices; i++) {
        vertices[i] = Vertex::create(points[i]);
    }
    EdgeSPtr edges[num_edges];
    edges[0] = Edge::create(vertices[0], vertices[1]);
    edges[1] = Edge::create(vertices[1], vertices[2]);
    edges[2] = Edge::create(vertices[2], vertices[0]);
    edges[3] = Edge::create(vertices[2], vertices[3]);
    edges[4] = Edge::create(vertices[3], vertices[0]);
    edges[5] = Edge::create(vertices[3], vertices[1]);
    FacetSPtr facets[num_facets];
    facets[0] = Facet::create(3, edges);
    facets[1] = Facet::create(3, &(edges[2]));
    EdgeSPtr edges2[] = {edges[4], edges[5], edges[0]};
    facets[2] = Facet::create(3, edges2);
    EdgeSPtr edges3[] = {edges[5], edges[3], edges[1]};
    facets[3] = Facet::create(3, edges3);
    PolyhedronSPtr polyhedron = Polyhedron::create(num_facets, facets);
    BOOST_CHECK(polyhedron->isConsistent());
    PolyhedronDAOSPtr dao_polyhedron = DAOFactory::getPolyhedronDAO();
    dao_polyhedron->insert(polyhedron);
    PolyhedronSPtr result = dao_polyhedron->find(polyhedron->getID());
    BOOST_CHECK(result->isConsistent());
    BOOST_CHECK_EQUAL(polyhedron->vertices().size(), result->vertices().size());
    BOOST_CHECK_EQUAL(polyhedron->edges().size(), result->edges().size());
    BOOST_CHECK_EQUAL(polyhedron->facets().size(), result->facets().size());
    const double e = 0.001;
    list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    list<VertexSPtr>::iterator it_vr = result->vertices().begin();
    while (it_v != polyhedron->vertices().end() && it_vr != result->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        VertexSPtr vertex_r = *it_vr++;
        BOOST_CHECK_CLOSE(vertex->getX(), vertex_r->getX(), e);
        BOOST_CHECK_CLOSE(vertex->getY(), vertex_r->getY(), e);
        BOOST_CHECK_CLOSE(vertex->getZ(), vertex_r->getZ(), e);
    }
    dao_polyhedron->del(polyhedron);
}

BOOST_AUTO_TEST_SUITE_END()

