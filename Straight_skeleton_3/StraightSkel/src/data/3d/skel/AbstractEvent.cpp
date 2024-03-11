/**
 * @file   data/3d/skel/AbstractEvent.cpp
 * @author Gernot Walzl
 * @date   2012-03-27
 */

#include "data/3d/skel/AbstractEvent.h"

#include "debug.h"
#include "util/StringFactory.h"
#include <sstream>

namespace data { namespace _3d { namespace skel {

AbstractEvent::AbstractEvent() {
    id_ = -1;
}

AbstractEvent::~AbstractEvent() {
    polyhedron_result_.reset();
    skel_.reset();
}

PolyhedronSPtr AbstractEvent::getPolyhedronResult() const {
    //DEBUG_SPTR(polyhedron_result_);
    return polyhedron_result_;
}

void AbstractEvent::setPolyhedronResult(PolyhedronSPtr polyhedron) {
    this->polyhedron_result_ = polyhedron;
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
        case SAVE_OFFSET_EVENT:
            sstr << "SaveOffsetEvent";
            break;
        case EDGE_EVENT:
            sstr << "EdgeEvent";
            break;
        case EDGE_MERGE_EVENT:
            sstr << "EdgeMergeEvent";
            break;
        case TRIANGLE_EVENT:
            sstr << "TriangleEvent";
            break;
        case DBL_EDGE_MERGE_EVENT:
            sstr << "DblEdgeMergeEvent";
            break;
        case DBL_TRIANGLE_EVENT:
            sstr << "DblTriangleEvent";
            break;
        case TETRAHEDRON_EVENT:
            sstr << "TetrahedronEvent";
            break;
        case VERTEX_EVENT:
            sstr << "VertexEvent";
            break;
        case FLIP_VERTEX_EVENT:
            sstr << "FlipVertexEvent";
            break;
        case SURFACE_EVENT:
            sstr << "SurfaceEvent";
            break;
        case POLYHEDRON_SPLIT_EVENT:
            sstr << "PolyhedronSplitEvent";
            break;
        case SPLIT_MERGE_EVENT:
            sstr << "SplitMergeEvent";
            break;
        case EDGE_SPLIT_EVENT:
            sstr << "EdgeSplitEvent";
            break;
        case PIERCE_EVENT:
            sstr << "PierceEvent";
            break;
        default:
            sstr << "AbstractEvent";
    }
    sstr << "(offset=" << util::StringFactory::fromDouble(getOffset()) << ")";
    return sstr.str();
}

} } }
