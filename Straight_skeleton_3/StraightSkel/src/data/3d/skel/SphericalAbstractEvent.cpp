/**
 * @file   data/3d/skel/SphericalAbstractEvent.cpp
 * @author Gernot Walzl
 * @date   2012-11-28
 */

#include "data/3d/skel/SphericalAbstractEvent.h"

#include "debug.h"
#include "util/StringFactory.h"
#include <sstream>

namespace data { namespace _3d { namespace skel {

SphericalAbstractEvent::SphericalAbstractEvent() {
    // intentionally does nothing
}

SphericalAbstractEvent::~SphericalAbstractEvent() {
    polygon_result_.reset();
    skel_.reset();
}

SphericalPolygonSPtr SphericalAbstractEvent::getPolygonResult() const {
    DEBUG_SPTR(polygon_result_);
    return polygon_result_;
}

void SphericalAbstractEvent::setPolygonResult(SphericalPolygonSPtr polygon) {
    this->polygon_result_ = polygon;
}

SphericalSkeletonSPtr SphericalAbstractEvent::getSkel() const {
    DEBUG_WPTR(skel_);
    if (this->skel_.expired())
        return SphericalSkeletonSPtr();
    else
        return SphericalSkeletonSPtr(this->skel_);
}

void SphericalAbstractEvent::setSkel(SphericalSkeletonSPtr skel) {
    this->skel_ = skel;
}

std::list<SphericalAbstractEventSPtr>::iterator SphericalAbstractEvent::getListIt() const {
    return this->list_it_;
}

void SphericalAbstractEvent::setListIt(std::list<SphericalAbstractEventSPtr>::iterator list_it) {
    this->list_it_ = list_it;
}

void SphericalAbstractEvent::setHighlight(bool highlight) {
    // to be implemented by inherited class
}

int SphericalAbstractEvent::getType() const {
    return this->type_;
}

std::string SphericalAbstractEvent::toString() const {
    std::stringstream sstr;
    switch (getType()) {
        case CONST_OFFSET_EVENT:
            sstr << "ConstOffsetEvent";
            break;
        case EDGE_EVENT:
            sstr << "EdgeEvent";
            break;
        case SPLIT_EVENT:
            sstr << "SplitEvent";
            break;
        case TRIANGLE_EVENT:
            sstr << "TriangleEvent";
            break;
        case DBL_EDGE_EVENT:
            sstr << "DblEdgeEvent";
            break;
        case LEAVE_EVENT:
            sstr << "LeaveEvent";
            break;
        case RETURN_EVENT:
            sstr << "ReturnEvent";
            break;
        case DBL_LEAVE_EVENT:
            sstr << "DblLeaveEvent";
            break;
        case DBL_RETURN_EVENT:
            sstr << "DblReturnEvent";
            break;
        case VERTEX_EVENT:
            sstr << "VertexEvent";
            break;
        case EDGE_MERGE_EVENT:
            sstr << "EdgeMergeEvent";
            break;
        case INVERSION_EVENT:
            sstr << "InversionEvent";
            break;
        default:
            sstr << "AbstractEvent";
    }
    sstr << "(offset=" << util::StringFactory::fromDouble(getOffset()) << ")";
    return sstr.str();
}

} } }
