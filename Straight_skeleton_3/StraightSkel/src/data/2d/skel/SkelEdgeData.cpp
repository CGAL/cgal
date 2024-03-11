/**
 * @file   data/2d/skel/SkelEdgeData.cpp
 * @author Gernot Walzl
 * @date   2012-05-12
 */

#include "data/2d/skel/SkelEdgeData.h"

#include "data/2d/Edge.h"
#include "debug.h"

namespace data { namespace _2d { namespace skel {

SkelEdgeData::SkelEdgeData() {
    speed_ = 1.0;
}

SkelEdgeData::~SkelEdgeData() {
    // intentionally does nothing
}

SkelEdgeDataSPtr SkelEdgeData::create(EdgeSPtr edge) {
    SkelEdgeDataSPtr result = SkelEdgeDataSPtr(new SkelEdgeData());
    result->setEdge(edge);
    result->setEdgeOrigin(edge);
    edge->setData(result);
    return result;
}

EdgeSPtr SkelEdgeData::getOffsetEdge() const {
    DEBUG_WPTR(offset_edge_);
    if (this->offset_edge_.expired())
        return EdgeSPtr();
    else
        return EdgeSPtr(this->offset_edge_);
}

void SkelEdgeData::setOffsetEdge(EdgeSPtr offset_edge) {
    this->offset_edge_ = offset_edge;
}

EdgeSPtr SkelEdgeData::getEdgeOrigin() const {
    DEBUG_WPTR(edge_origin_);
    if (this->edge_origin_.expired())
        return EdgeSPtr();
    else
        return EdgeSPtr(this->edge_origin_);
}

void SkelEdgeData::setEdgeOrigin(EdgeSPtr edge_origin) {
    this->edge_origin_ = edge_origin;
}

double SkelEdgeData::getSpeed() const {
    return speed_;
}

void SkelEdgeData::setSpeed(double speed) {
    speed_ = speed;
}

} } }
