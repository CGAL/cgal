/**
 * @file   db/2d/AbstractFile.cpp
 * @author Gernot Walzl
 * @date   2013-12-18
 */

#include "db/2d/AbstractFile.h"

#include "data/2d/KernelFactory.h"
#include "data/2d/Vertex.h"
#include "data/2d/Edge.h"
#include "data/2d/Polygon.h"
#include "debug.h"
#include <cmath>

namespace db { namespace _2d {

AbstractFile::AbstractFile() {
    // intentionally does nothing
}

AbstractFile::~AbstractFile() {
    // intentionally does nothing
}

bool AbstractFile::hasCollinearEdges(VertexSPtr vertex, double epsilon) {
    bool result = false;
    EdgeSPtr edge_in = vertex->getEdgeIn();
    EdgeSPtr edge_out = vertex->getEdgeOut();
    if (edge_in && edge_out) {
        Vector2SPtr normal_l = KernelFactory::createVector2(edge_in->line());
        Vector2SPtr normal_r = KernelFactory::createVector2(edge_out->line());
        double length_l = 0.0;
        double length_r = 0.0;
        for (unsigned int i = 0; i < 2; i++) {
            length_l += (*normal_l)[i] * (*normal_l)[i];
            length_r += (*normal_r)[i] * (*normal_r)[i];
        }
        length_l = sqrt(length_l);
        length_r = sqrt(length_r);
        double diff = 0.0;
        double length_diff = 0.0;
        for (unsigned int i = 0; i < 2; i++) {
            diff = ((*normal_l)[i]/length_l) - ((*normal_r)[i]/length_r);
            length_diff += diff*diff;
        }
        length_diff = sqrt(length_diff);
        result = (length_diff < epsilon);
    }
    return result;
}

int AbstractFile::mergeCollinearEdges(PolygonSPtr polygon, double epsilon) {
    int result = 0;
    std::list<VertexSPtr> vertices_toremove;
    std::list<VertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        if (hasCollinearEdges(vertex, epsilon)) {
            vertices_toremove.push_back(vertex);
        }
    }
    if (vertices_toremove.size() > 0) {
        DEBUG_PRINT("Adjacent edges of the following vertices are detected to be collinear and will be merged.");
    }
    it_v = vertices_toremove.begin();
    while (it_v != vertices_toremove.end()) {
        VertexSPtr vertex = *it_v++;
        DEBUG_VAR(vertex->toString());
        EdgeSPtr edge_in = vertex->getEdgeIn();
        EdgeSPtr edge_out = vertex->getEdgeOut();
        VertexSPtr vertex_dst = vertex->next();
        edge_in->setVertexDst(vertex_dst);
        vertex_dst->setEdgeIn(edge_in);
        vertex->setEdgeIn(EdgeSPtr());
        polygon->removeEdge(edge_out);
        polygon->removeVertex(vertex);
        result++;
    }
    return result;
}

} }
