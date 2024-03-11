/**
 * @file   data/2d/EdgeData.cpp
 * @author Gernot Walzl
 * @date   2011-11-22
 */

#include "data/2d/EdgeData.h"

#include "data/2d/Edge.h"
#include "debug.h"

namespace data { namespace _2d {

EdgeData::EdgeData() {
    highlight_ = false;
}

EdgeData::~EdgeData() {
    // intentionally does nothing
}

EdgeSPtr EdgeData::getEdge() const {
    DEBUG_WPTR(edge_);
    if (this->edge_.expired())
        return EdgeSPtr();
    else
        return EdgeSPtr(this->edge_);
}

void EdgeData::setEdge(EdgeSPtr edge) {
    this->edge_ = edge;
}

bool EdgeData::isHighlight() const {
    return highlight_;
}

void EdgeData::setHighlight(bool highlight) {
    highlight_ = highlight;
}

} }
