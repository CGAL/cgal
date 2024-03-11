/**
 * @file   data/2d/skel/AbstractEvent.cpp
 * @author Gernot Walzl
 * @date   2012-02-02
 */

#include "data/2d/skel/AbstractEvent.h"

#include "debug.h"
#include "data/2d/skel/Node.h"
#include "util/StringFactory.h"
#include <sstream>

namespace data { namespace _2d { namespace skel {

AbstractEvent::AbstractEvent() {
    id_ = -1;
}

AbstractEvent::~AbstractEvent() {
    polygon_result_.reset();
    skel_.reset();
}

PolygonSPtr AbstractEvent::getPolygonResult() const {
    //DEBUG_SPTR(polygon_result_);
    return polygon_result_;
}

void AbstractEvent::setPolygonResult(PolygonSPtr polygon) {
    this->polygon_result_ = polygon;
}

StraightSkeletonSPtr AbstractEvent::getSkel() const {
    DEBUG_WPTR(skel_);
    if (this->skel_.expired())
        return StraightSkeletonSPtr();
    else
        return StraightSkeletonSPtr(this->skel_);
}

void AbstractEvent::setSkel(StraightSkeletonSPtr skel) {
    this->skel_ = skel;
}

std::list<AbstractEventSPtr>::iterator AbstractEvent::getListIt() const {
    return this->list_it_;
}

void AbstractEvent::setListIt(std::list<AbstractEventSPtr>::iterator list_it) {
    this->list_it_ = list_it;
}

int AbstractEvent::getID() const {
    return this->id_;
}

void AbstractEvent::setID(int id) {
    this->id_ = id;
}

void AbstractEvent::setHighlight(bool highlight) {
    // to be implemented by inherited class
}

int AbstractEvent::getType() const {
    return this->type_;
}

std::string AbstractEvent::toString() const {
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
        default:
            sstr << "AbstractEvent";
    }
    sstr << "(offset=" << util::StringFactory::fromDouble(getOffset()) << ")";
    return sstr.str();
}

} } }
