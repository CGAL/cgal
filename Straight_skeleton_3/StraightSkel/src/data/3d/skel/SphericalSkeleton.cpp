/**
 * @file   data/3d/skel/SphericalSkeleton.cpp
 * @author Gernot Walzl
 * @date   2012-11-28
 */

#include "data/3d/skel/SphericalSkeleton.h"

#include "debug.h"
#include "data/3d/KernelFactory.h"
#include "data/3d/skel/CircularNode.h"
#include "data/3d/skel/CircularArc.h"
#include "data/3d/skel/SphericalAbstractEvent.h"
#include "data/3d/skel/SphericalEdgeEvent.h"
#include "data/3d/skel/SphericalSplitEvent.h"
#include "data/3d/skel/SphericalTriangleEvent.h"
#include "util/StringFactory.h"
#include <sstream>

namespace data { namespace _3d { namespace skel {

SphericalSkeleton::SphericalSkeleton(Sphere3SPtr sphere) {
    sphere_ = sphere;
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    Vector3SPtr normal = KernelFactory::createVector3(0.0, 0.0, 1.0);
    rot_axes_plane_ = KernelFactory::createPlane3(p_center, normal);
}

SphericalSkeleton::~SphericalSkeleton() {
    sphere_.reset();
    rot_axes_plane_.reset();
}

SphericalSkeletonSPtr SphericalSkeleton::create(Sphere3SPtr sphere) {
    return SphericalSkeletonSPtr(new SphericalSkeleton(sphere));
}


Sphere3SPtr SphericalSkeleton::getSphere() const {
    return sphere_;
}

Plane3SPtr SphericalSkeleton::getRotAxesPlane() const {
    return rot_axes_plane_;
}

void SphericalSkeleton::setRotAxesPlane(Plane3SPtr rot_axes_plane) {
    rot_axes_plane_ = rot_axes_plane;
}


void SphericalSkeleton::addEvent(SphericalAbstractEventSPtr event) {
    std::list<SphericalAbstractEventSPtr>::iterator it = events_.insert(events_.end(), event);
    event->setSkel(shared_from_this());
    event->setListIt(it);
    CircularNodeSPtr node = CircularNodeSPtr();
    if (event->getType() == SphericalAbstractEvent::EDGE_EVENT) {
        SphericalEdgeEventSPtr edge_event =
                std::dynamic_pointer_cast<SphericalEdgeEvent>(event);
        node = edge_event->getNode();
    } else if (event->getType() == SphericalAbstractEvent::SPLIT_EVENT) {
        SphericalSplitEventSPtr split_event =
                std::dynamic_pointer_cast<SphericalSplitEvent>(event);
        node = split_event->getNode();
    } else if (event->getType() == SphericalAbstractEvent::TRIANGLE_EVENT) {
        SphericalTriangleEventSPtr triangle_event =
                std::dynamic_pointer_cast<SphericalTriangleEvent>(event);
        node = triangle_event->getNode();
    }
    if (node) {
        if (node->getSkel() != shared_from_this()) {
            addNode(node);
        }
    }
}

bool SphericalSkeleton::removeEvent(SphericalAbstractEventSPtr event) {
    bool result = false;
    if (event->getSkel() == shared_from_this()) {
        events_.erase(event->getListIt());
        event->setSkel(SphericalSkeletonSPtr());
        event->setListIt(std::list<SphericalAbstractEventSPtr>::iterator());
        CircularNodeSPtr node = CircularNodeSPtr();
        if (event->getType() == SphericalAbstractEvent::EDGE_EVENT) {
            SphericalEdgeEventSPtr edge_event =
                    std::dynamic_pointer_cast<SphericalEdgeEvent>(event);
            node = edge_event->getNode();
        } else if (event->getType() == SphericalAbstractEvent::SPLIT_EVENT) {
            SphericalSplitEventSPtr split_event =
                    std::dynamic_pointer_cast<SphericalSplitEvent>(event);
            node = split_event->getNode();
        } else if (event->getType() == SphericalAbstractEvent::TRIANGLE_EVENT) {
            SphericalTriangleEventSPtr triangle_event =
                    std::dynamic_pointer_cast<SphericalTriangleEvent>(event);
            node = triangle_event->getNode();
        }
        if (node) {
            result = removeNode(node);
        }
    }
    return result;
}

void SphericalSkeleton::addNode(CircularNodeSPtr node) {
    std::list<CircularNodeSPtr>::iterator it = nodes_.insert(nodes_.end(), node);
    node->setSkel(shared_from_this());
    node->setListIt(it);
}

bool SphericalSkeleton::removeNode(CircularNodeSPtr node) {
    bool result = false;
    if (node->getSkel() == shared_from_this()) {
        nodes_.erase(node->getListIt());
        node->setSkel(SphericalSkeletonSPtr());
        node->setListIt(std::list<CircularNodeSPtr>::iterator());
        result = true;
    }
    return result;
}

void SphericalSkeleton::addArc(CircularArcSPtr arc) {
    std::list<CircularArcSPtr>::iterator it = arcs_.insert(arcs_.end(), arc);
    arc->setSkel(shared_from_this());
    arc->setListIt(it);
}

bool SphericalSkeleton::removeArc(CircularArcSPtr arc) {
    bool result = false;
    if (arc->getSkel() == shared_from_this()) {
        arcs_.erase(arc->getListIt());
        arc->setSkel(SphericalSkeletonSPtr());
        arc->setListIt(std::list<CircularArcSPtr>::iterator());
        result = true;
    }
    return result;
}

SharedMutex& SphericalSkeleton::mutex() {
    return this->mutex_;
}

std::list<SphericalAbstractEventSPtr>& SphericalSkeleton::events() {
    return this->events_;
}

std::list<CircularNodeSPtr>& SphericalSkeleton::nodes() {
    return this->nodes_;
}

std::list<CircularArcSPtr>& SphericalSkeleton::arcs() {
    return this->arcs_;
}

bool SphericalSkeleton::isConsistent() const {
    bool result = true;
    std::list<CircularNodeSPtr>::const_iterator it_n = nodes_.begin();
    while (it_n != nodes_.end()) {
        CircularNodeSPtr node = *it_n++;
        if (node->getSkel() != shared_from_this()) {
            DEBUG_VAR(node->toString());
            result = false;
            break;
        }
        unsigned int num_arcs = 0;
        std::list<CircularArcWPtr>::const_iterator it_a = node->arcs().begin();
        while (it_a != node->arcs().end()) {
            CircularArcWPtr arc_wptr = *it_a++;
            if (arc_wptr.expired()) {
                DEBUG_VAR(node->toString());
            } else {
                CircularArcSPtr arc = CircularArcSPtr(arc_wptr);
                num_arcs++;
                if (node != arc->getNodeSrc() && node != arc->getNodeDst()) {
                    DEBUG_VAR(node->toString());
                    DEBUG_VAR(arc->toString());
                    result = false;
                    break;
                }
            }
        }
        if (!(num_arcs == 1 || num_arcs == 3)) {
            DEBUG_VAL("Warning: CircularNode does not have 1 or 3 CircularArcs.");
            DEBUG_VAR(node->toString());
            DEBUG_VAR(num_arcs);
        }
    }
    std::list<CircularArcSPtr>::const_iterator it_a = arcs_.begin();
    while (it_a != arcs_.end()) {
        CircularArcSPtr arc = *it_a++;
        CircularArcWPtr arc_wptr(arc);
        if (arc->getSkel() != shared_from_this()) {
            DEBUG_VAR(arc->toString());
            result = false;
            break;
        }
        std::list<CircularArcWPtr> warcs = arc->getNodeSrc()->arcs();
        if (warcs.end() == std::find(warcs.begin(), warcs.end(), arc_wptr)) {
            DEBUG_VAR(arc->toString());
            DEBUG_VAR(arc->getNodeSrc()->toString());
            result = false;
            break;
        }
        if (arc->hasNodeDst()) {
            warcs = arc->getNodeDst()->arcs();
            if (warcs.end() == std::find(warcs.begin(), warcs.end(), arc_wptr)) {
                DEBUG_VAR(arc->toString());
                DEBUG_VAR(arc->getNodeDst()->toString());
                result = false;
                break;
            }
        }
    }
    return result;
}

int SphericalSkeleton::countEvents(int type) const {
    int result = 0;
    std::list<SphericalAbstractEventSPtr>::const_iterator it_e = events_.begin();
    while (it_e != events_.end()) {
        SphericalAbstractEventSPtr event = *it_e++;
        if (event->getType() == type) {
            result += 1;
        }
    }
    return result;
}

double SphericalSkeleton::getRadius() const {
    double result = 0.0;
    if (sphere_) {
#ifdef USE_CGAL
        result = CGAL::sqrt(sphere_->squared_radius());
#else
        result = sphere_->getRadius();
#endif
    }
    return result;
}

std::string SphericalSkeleton::toString() const {
    std::stringstream sstr;
    sstr << "SphericalSkeleton(";
    sstr << util::StringFactory::fromPointer(this);
    sstr << "," << std::endl;
    sstr << "Nodes:  " << nodes_.size() << std::endl;
    sstr << "Arcs:   " << arcs_.size() << std::endl;
    sstr << "Events: " << events_.size() << std::endl;
    unsigned int num_events = 0;
    num_events = countEvents(SphericalAbstractEvent::CONST_OFFSET_EVENT);
    if (num_events > 0) {
        sstr << "  ConstOffsetEvents: " << num_events << std::endl;
    }
    num_events = countEvents(SphericalAbstractEvent::EDGE_EVENT);
    if (num_events > 0) {
        sstr << "  EdgeEvents:        " << num_events << std::endl;
    }
    num_events = countEvents(SphericalAbstractEvent::SPLIT_EVENT);
    if (num_events > 0) {
        sstr << "  SplitEvents:       " << num_events << std::endl;
    }
    num_events = countEvents(SphericalAbstractEvent::TRIANGLE_EVENT);
    if (num_events > 0) {
        sstr << "  TriangleEvents:    " << num_events << std::endl;
    }
    num_events = countEvents(SphericalAbstractEvent::DBL_EDGE_EVENT);
    if (num_events > 0) {
        sstr << "  DblEdgeEvents:     " << num_events << std::endl;
    }
    num_events = countEvents(SphericalAbstractEvent::LEAVE_EVENT);
    if (num_events > 0) {
        sstr << "  LeaveEvents:       " << num_events << std::endl;
    }
    num_events = countEvents(SphericalAbstractEvent::RETURN_EVENT);
    if (num_events > 0) {
        sstr << "  ReturnEvents:      " << num_events << std::endl;
    }
    num_events = countEvents(SphericalAbstractEvent::DBL_LEAVE_EVENT);
    if (num_events > 0) {
        sstr << "  DblLeaveEvents:     " << num_events << std::endl;
    }
    num_events = countEvents(SphericalAbstractEvent::DBL_RETURN_EVENT);
    if (num_events > 0) {
        sstr << "  DblReturnEvents:    " << num_events << std::endl;
    }
    num_events = countEvents(SphericalAbstractEvent::VERTEX_EVENT);
    if (num_events > 0) {
        sstr << "  VertexEvents:      " << num_events << std::endl;
    }
    num_events = countEvents(SphericalAbstractEvent::EDGE_MERGE_EVENT);
    if (num_events > 0) {
        sstr << "  EdgeMergeEvents:   " << num_events << std::endl;
    }
    num_events = countEvents(SphericalAbstractEvent::INVERSION_EVENT);
    if (num_events > 0) {
        sstr << "  InversionEvents:   " << num_events << std::endl;
    }
    sstr << ")" << std::endl;
    return sstr.str();
}

} } }

