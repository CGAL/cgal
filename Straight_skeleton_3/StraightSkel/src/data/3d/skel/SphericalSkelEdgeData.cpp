/**
 * @file   data/3d/skel/SphericalSkelEdgeData.cpp
 * @author Gernot Walzl
 * @date   2012-11-30
 */

#include "data/3d/skel/SphericalSkelEdgeData.h"

#include "debug.h"
#include "data/3d/CircularEdge.h"

namespace data { namespace _3d { namespace skel {

SphericalSkelEdgeData::SphericalSkelEdgeData() {
    speed_ = 1.0;
}

SphericalSkelEdgeData::~SphericalSkelEdgeData() {
    // intentionally does nothing
}

SphericalSkelEdgeDataSPtr SphericalSkelEdgeData::create(CircularEdgeSPtr edge) {
    SphericalSkelEdgeDataSPtr result = SphericalSkelEdgeDataSPtr(new SphericalSkelEdgeData());
    result->setEdge(edge);
    edge->setData(result);
    return result;
}

CircularEdgeSPtr SphericalSkelEdgeData::getOffsetEdge() const {
    DEBUG_WPTR(offset_edge_);
    if (this->offset_edge_.expired())
        return CircularEdgeSPtr();
    else
        return CircularEdgeSPtr(this->offset_edge_);
}

void SphericalSkelEdgeData::setOffsetEdge(CircularEdgeSPtr offset_edge) {
    this->offset_edge_ = offset_edge;
}

double SphericalSkelEdgeData::getSpeed() const {
    return speed_;
}

void SphericalSkelEdgeData::setSpeed(double speed) {
    speed_ = speed;
}

FacetSPtr SphericalSkelEdgeData::getFacetOrigin() const {
    DEBUG_WPTR(facet_origin_);
    if (this->facet_origin_.expired())
        return FacetSPtr();
    else
        return FacetSPtr(this->facet_origin_);
}

void SphericalSkelEdgeData::setFacetOrigin(FacetSPtr facet_origin) {
    this->facet_origin_ = facet_origin;
}

Line3SPtr SphericalSkelEdgeData::getRotationAxis() const {
    //DEBUG_SPTR(rotation_axis_);
    return this->rotation_axis_;
}

void SphericalSkelEdgeData::setRotationAxis(Line3SPtr rotation_axis) {
    this->rotation_axis_ = rotation_axis;
}

} } }
