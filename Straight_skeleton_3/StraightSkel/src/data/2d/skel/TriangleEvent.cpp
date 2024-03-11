/**
 * @file   data/2d/skel/TriangleEvent.cpp
 * @author Gernot Walzl
 * @date   2013-05-06
 */

#include "data/2d/skel/TriangleEvent.h"

#include "data/2d/Vertex.h"
#include "data/2d/Edge.h"
#include "data/2d/skel/SkelVertexData.h"
#include "data/2d/skel/SkelEdgeData.h"

namespace data { namespace _2d { namespace skel {

TriangleEvent::TriangleEvent() {
    type_ = AbstractEvent::TRIANGLE_EVENT;
}

TriangleEvent::~TriangleEvent() {
    // intentionally does nothing
}

TriangleEventSPtr TriangleEvent::create() {
    TriangleEventSPtr result = TriangleEventSPtr(new TriangleEvent());
    return result;
}

void TriangleEvent::getVertices(VertexSPtr out[3]) const {
    for (unsigned int i = 0; i < 3; i++) {
        out[i] = VertexSPtr();
    }
    out[0] = edge_->getVertexSrc();
    out[1] = edge_->getVertexDst();
    out[2] = edge_->next()->getVertexDst();
}

void TriangleEvent::getEdges(EdgeSPtr out[3]) const {
    for (unsigned int i = 0; i < 3; i++) {
        out[i] = EdgeSPtr();
    }
    out[0] = edge_;
    out[1] = edge_->next();
    out[2] = edge_->prev();
}

void TriangleEvent::setHighlight(bool highlight) {
    VertexSPtr vertices[3];
    getVertices(vertices);
    for (unsigned int i = 0; i < 3; i++) {
        if (!vertices[i]->hasData()) {
            SkelVertexData::create(vertices[i]);
        }
        vertices[i]->getData()->setHighlight(highlight);
    }
    EdgeSPtr edges[3];
    getEdges(edges);
    for (unsigned int i = 0; i < 3; i++) {
        if (!edges[i]->hasData()) {
            SkelEdgeData::create(edges[i]);
        }
        edges[i]->getData()->setHighlight(highlight);
    }
}

} } }
