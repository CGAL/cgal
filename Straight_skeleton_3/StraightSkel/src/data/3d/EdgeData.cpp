/**
 * @file   data/3d/EdgeData.cpp
 * @author Gernot Walzl
 * @date   2011-11-22
 */

#include "data/3d/EdgeData.h"

#include "data/3d/Edge.h"

namespace data { namespace _3d {

EdgeData::EdgeData() {
    highlight_ = false;
}

EdgeData::~EdgeData() {
    // intentionally does nothing
}

EdgeDataSPtr EdgeData::create(EdgeSPtr edge) {
    EdgeDataSPtr result = EdgeDataSPtr(new EdgeData());
    result->setEdge(edge);
    edge->setData(result);
    return result;
}

EdgeSPtr EdgeData::getEdge() const {
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
