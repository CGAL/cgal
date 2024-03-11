/**
 * @file   data/3d/CircularEdgeData.cpp
 * @author Gernot Walzl
 * @date   2012-11-30
 */

#include "data/3d/CircularEdgeData.h"

#include "debug.h"
#include "data/3d/CircularEdge.h"

namespace data { namespace _3d {

CircularEdgeData::CircularEdgeData() {
    highlight_ = false;
}

CircularEdgeData::~CircularEdgeData() {
    // intentionally does nothing
}

CircularEdgeSPtr CircularEdgeData::getEdge() const {
    if (this->edge_.expired())
        return CircularEdgeSPtr();
    else
        return CircularEdgeSPtr(this->edge_);
}

void CircularEdgeData::setEdge(CircularEdgeSPtr edge) {
    this->edge_ = edge;
}

bool CircularEdgeData::isHighlight() const {
    return highlight_;
}

void CircularEdgeData::setHighlight(bool highlight) {
    highlight_ = highlight;
}

} }
