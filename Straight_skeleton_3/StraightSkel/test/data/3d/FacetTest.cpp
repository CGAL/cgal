#include <boost/test/unit_test.hpp>

#include "data/3d/Facet.h"
#include "data/3d/ptrs.h"
#include "data/3d/KernelFactory.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"

using namespace data::_3d;

BOOST_AUTO_TEST_SUITE(FacetTest)

BOOST_AUTO_TEST_CASE(testConstructor) {
    const unsigned int num_vertices = 3;
    const unsigned int num_edges = 3;
    Point3SPtr points[num_vertices];
    points[0] = KernelFactory::createPoint3(-1.0, -1.0, -1.0);
    points[1] = KernelFactory::createPoint3(1.0, 1.0, -1.0);
    points[2] = KernelFactory::createPoint3(1.0, -1.0, 1.0);
    VertexSPtr vertices[num_vertices];
    for (unsigned int i = 0; i < num_vertices; i++) {
        vertices[i] = Vertex::create(points[i]);
    }
    EdgeSPtr edges[num_edges];
    edges[0] = Edge::create(vertices[0], vertices[1]);
    edges[1] = Edge::create(vertices[1], vertices[2]);
    edges[2] = Edge::create(vertices[2], vertices[0]);
    FacetSPtr result = Facet::create(num_edges, edges);

    BOOST_CHECK_EQUAL(num_vertices, result->vertices().size());
    BOOST_CHECK_EQUAL(num_edges, result->edges().size());
    for (unsigned int i = 0; i < num_vertices; i++) {
        BOOST_CHECK_EQUAL(1, vertices[i]->facets().size());
        BOOST_CHECK_EQUAL(2, vertices[i]->edges().size());
        FacetSPtr expected = vertices[i]->facets().front().lock();
        BOOST_CHECK(result == expected);
    }
    for (unsigned int i = 0; i < num_edges; i++) {
        BOOST_CHECK(result == edges[i]->getFacetL());
        BOOST_CHECK(!edges[i]->getFacetR());
    }
    unsigned int i = 0;
    list<EdgeSPtr>::iterator it_e = result->edges().begin();
    while(it_e != result->edges().end()) {
        EdgeSPtr edge = *it_e++;
        BOOST_CHECK_EQUAL(edges[i], edge);
        i++;
    }
}

BOOST_AUTO_TEST_CASE(testToPolygon) {
    const unsigned int num_vertices = 3;
    const unsigned int num_edges = 3;
    Point3SPtr points[num_vertices];
    points[0] = KernelFactory::createPoint3(1.0, 2.0, 1.0);
    points[1] = KernelFactory::createPoint3(6.0, 3.0, 1.0);
    points[2] = KernelFactory::createPoint3(4.0, 5.0, 1.0);
    VertexSPtr vertices[num_vertices];
    for (unsigned int i = 0; i < num_vertices; i++) {
        vertices[i] = Vertex::create(points[i]);
    }
    EdgeSPtr edges[num_edges];
    edges[0] = Edge::create(vertices[0], vertices[1]);
    edges[1] = Edge::create(vertices[1], vertices[2]);
    edges[2] = Edge::create(vertices[2], vertices[0]);
    FacetSPtr facet = Facet::create(num_edges, edges);
    data::_2d::PolygonSPtr result = facet->toPolygon();
    BOOST_CHECK_EQUAL(facet->vertices().size(), result->vertices().size());
    BOOST_CHECK_EQUAL(facet->edges().size(), result->edges().size());
    const double e = 0.001;
    list<VertexSPtr>::iterator it_v = facet->vertices().begin();
    list<data::_2d::VertexSPtr>::iterator it_v2 = result->vertices().begin();
    while(it_v != facet->vertices().end() && it_v2 != result->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        data::_2d::VertexSPtr vertex2 = *it_v2++;
        BOOST_CHECK_CLOSE(vertex->getX(), vertex2->getX(), e);
        BOOST_CHECK_CLOSE(vertex->getY(), vertex2->getY(), e);
    }

    points[0] = KernelFactory::createPoint3(1.0, 2.0, 1.0);
    points[1] = KernelFactory::createPoint3(6.0, 3.0, 6.0);
    points[2] = KernelFactory::createPoint3(4.0, 5.0, 4.0);
    for (unsigned int i = 0; i < num_vertices; i++) {
        vertices[i]->setPoint(points[i]);
    }
    facet->initPlane();
    result = facet->toPolygon();
    list<EdgeSPtr>::iterator it_e = facet->edges().begin();
    list<data::_2d::EdgeSPtr>::iterator it_e2 = result->edges().begin();
    while(it_e != facet->edges().end() && it_e2 != result->edges().end()) {
        EdgeSPtr edge = *it_e++;
        data::_2d::EdgeSPtr edge2 = *it_e2++;
        Vector3SPtr v3 = KernelFactory::createVector3(
                *(edge->getVertexDst()->getPoint()) - *(edge->getVertexSrc()->getPoint()));
        double length3 = 0.0;
        for (unsigned int i = 0; i < 3; i++) {
            length3 += (*v3)[i] * (*v3)[i];
        }
        data::_2d::Vector2SPtr v2 = data::_2d::KernelFactory::createVector2(
                *(edge2->getVertexDst()->getPoint()) - *(edge2->getVertexSrc()->getPoint()));
        double length2 = 0.0;
        for (unsigned int i = 0; i < 2; i++) {
            length2 += (*v2)[i] * (*v2)[i];
        }
        BOOST_CHECK_CLOSE(length3, length2, e);
    }
}

BOOST_AUTO_TEST_SUITE_END()
