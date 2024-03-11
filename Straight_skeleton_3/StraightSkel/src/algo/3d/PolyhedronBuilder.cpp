/**
 * @file   algo/3d/PolyhedronBuilder.cpp
 * @author Gernot Walzl
 * @date   2012-04-03
 */

#include "algo/3d/PolyhedronBuilder.h"

#include "algo/3d/KernelWrapper.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/Facet.h"
#include "data/3d/Edge.h"
#include "data/3d/Vertex.h"

namespace algo { namespace _3d {

PolyhedronBuilder::PolyhedronBuilder() {
    // intentionally does nothing
}

PolyhedronBuilder::~PolyhedronBuilder() {
    // intentionally does nothing
}

PolyhedronSPtr PolyhedronBuilder::makeTetrahedron(
        Point3SPtr p1, Point3SPtr p2, Point3SPtr p3, Point3SPtr p4) {
    const unsigned int num_vertices = 4;
    const unsigned int num_edges = 6;
    const unsigned int num_facets = 4;
    VertexSPtr vertices[num_vertices];
    Plane3SPtr plane = KernelFactory::createPlane3(p1, p2, p3);
    if (KernelWrapper::side(plane, p4) < 0) {  // everything is fine.
        vertices[0] = Vertex::create(p1);
        vertices[1] = Vertex::create(p2);
        vertices[2] = Vertex::create(p3);
        vertices[3] = Vertex::create(p4);
    } else {
        vertices[0] = Vertex::create(p1);
        vertices[2] = Vertex::create(p2);
        vertices[1] = Vertex::create(p3);
        vertices[3] = Vertex::create(p4);
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
    PolyhedronSPtr result = Polyhedron::create(num_facets, facets);
    return result;
}

} }
