#include <boost/test/unit_test.hpp>

#include "data/2d/Polygon.h"
#include <list>
#include "data/2d/ptrs.h"
#include "data/2d/KernelFactory.h"
#include "data/2d/Vertex.h"
#include "data/2d/Edge.h"

using namespace data::_2d;

BOOST_AUTO_TEST_SUITE(PolygonTest)

BOOST_AUTO_TEST_CASE(testConstructor) {
    unsigned int i = 0;
    const unsigned int num_points = 3;
    Point2SPtr p[num_points];
    p[0] = KernelFactory::createPoint2(0.0, 0.0);
    p[1] = KernelFactory::createPoint2(1.0, 0.0);
    p[2] = KernelFactory::createPoint2(1.0, 1.0);
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
    BOOST_CHECK(polygon->isConsistent());
    i = 0;
    for (list<VertexSPtr>::iterator it = polygon->vertices().begin();
            it != polygon->vertices().end(); it++) {
        BOOST_CHECK_EQUAL(*it, v[i++]);
    }
    BOOST_CHECK_EQUAL(num_points, i);
    i = 0;
    for (list<EdgeSPtr>::iterator it = polygon->edges().begin();
            it != polygon->edges().end(); it++) {
        BOOST_CHECK_EQUAL(*it, e[i++]);
    }
    BOOST_CHECK_EQUAL(num_points, i);
}

BOOST_AUTO_TEST_CASE(testList) {
    VertexSPtr v[3];
    v[0] = Vertex::create(KernelFactory::createPoint2(1.0, 2.0));
    v[1] = Vertex::create(KernelFactory::createPoint2(3.0, 4.0));
    v[2] = Vertex::create(KernelFactory::createPoint2(5.0, 6.0));
    PolygonSPtr polygon = Polygon::create();
    for (unsigned int i = 0; i < 3; i++) {
        polygon->addVertex(v[i]);
    }
    polygon->removeVertex(v[0]);
    polygon->removeVertex(v[2]);
    // BOOST_CHECK_EQUAL does not work here (Problem with iterators?)
    BOOST_CHECK(polygon->vertices().begin() == v[1]->getListIt());
    BOOST_CHECK_EQUAL(*(polygon->vertices().begin()), *(v[1]->getListIt()));
    BOOST_CHECK_EQUAL(*(polygon->vertices().begin()), v[1]);
}

BOOST_AUTO_TEST_CASE(testSortEdges) {
    VertexSPtr v[4];
    v[0] = Vertex::create(KernelFactory::createPoint2(1.0, 1.0));
    v[1] = Vertex::create(KernelFactory::createPoint2(2.0, 1.0));
    v[2] = Vertex::create(KernelFactory::createPoint2(2.0, 2.0));
    v[3] = Vertex::create(KernelFactory::createPoint2(1.0, 2.0));
    EdgeSPtr e[4];
    e[0] = Edge::create(v[0], v[1]);
    e[1] = Edge::create(v[1], v[2]);
    e[2] = Edge::create(v[2], v[3]);
    e[3] = Edge::create(v[3], v[0]);
    PolygonSPtr polygon = Polygon::create();
    polygon->addEdge(e[0]);
    polygon->addEdge(e[2]);
    polygon->addEdge(e[3]);
    polygon->addEdge(e[1]);
    BOOST_CHECK(polygon->isConsistent());
    polygon->sortEdges();
    list<EdgeSPtr>::iterator it_e = polygon->edges().begin();
    for (unsigned int i = 0; i < 4; i++) {
        EdgeSPtr edge = *it_e;
        BOOST_CHECK_EQUAL(e[i], edge);
        it_e++;
    }
}

BOOST_AUTO_TEST_SUITE_END()

