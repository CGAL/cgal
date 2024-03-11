#include <boost/test/unit_test.hpp>

#include "data/3d/Polyhedron.h"
#include "data/3d/ptrs.h"
#include "data/3d/KernelFactory.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "debug.h"


using namespace data::_3d;

BOOST_AUTO_TEST_SUITE(PolyhedronTest)

BOOST_AUTO_TEST_CASE(testConstructor) {
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
}

BOOST_AUTO_TEST_SUITE_END()
