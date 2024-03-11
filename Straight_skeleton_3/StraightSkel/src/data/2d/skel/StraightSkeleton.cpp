/**
 * @file   data/2d/skel/StraightSkeleton.cpp
 * @author Gernot Walzl
 * @date   2011-11-21
 */

#include "data/2d/skel/StraightSkeleton.h"

#include "debug.h"
#include "data/2d/skel/Node.h"
#include "data/2d/skel/Arc.h"
#include "data/2d/skel/AbstractEvent.h"
#include "data/2d/skel/EdgeEvent.h"
#include "data/2d/skel/SplitEvent.h"
#include "data/2d/skel/TriangleEvent.h"
#include "util/StringFactory.h"
#include <sstream>

namespace data { namespace _2d { namespace skel {

StraightSkeleton::StraightSkeleton() {
    id_ = -1;
}

StraightSkeleton::~StraightSkeleton() {
    events_.clear();
    arcs_.clear();
    nodes_.clear();
}

StraightSkeletonSPtr StraightSkeleton::create() {
    return StraightSkeletonSPtr(new StraightSkeleton());
}

void StraightSkeleton::addEvent(AbstractEventSPtr event) {
    std::list<AbstractEventSPtr>::iterator it = events_.insert(events_.end(), event);
    event->setSkel(shared_from_this());
    event->setListIt(it);
    NodeSPtr node;
    if (event->getType() == AbstractEvent::EDGE_EVENT) {
        EdgeEventSPtr edge_event = std::dynamic_pointer_cast<EdgeEvent>(event);
        node = edge_event->getNode();
    } else if (event->getType() == AbstractEvent::SPLIT_EVENT) {
        SplitEventSPtr split_event = std::dynamic_pointer_cast<SplitEvent>(event);
        node = split_event->getNode();
    } else if (event->getType() == AbstractEvent::TRIANGLE_EVENT) {
        TriangleEventSPtr triangle_event = std::dynamic_pointer_cast<TriangleEvent>(event);
        node = triangle_event->getNode();
    }
    if (node) {
        if (node->getSkel() != shared_from_this()) {
            this->addNode(node);
        }
    }
}

bool StraightSkeleton::removeEvent(AbstractEventSPtr event) {
    bool result = false;
    if (event->getSkel() == shared_from_this()) {
        events_.erase(event->getListIt());
        event->setSkel(StraightSkeletonSPtr());
        event->setListIt(std::list<AbstractEventSPtr>::iterator());
        NodeSPtr node;
        if (event->getType() == AbstractEvent::EDGE_EVENT) {
            EdgeEventSPtr edge_event = std::dynamic_pointer_cast<EdgeEvent>(event);
            node = edge_event->getNode();
        } else if (event->getType() == AbstractEvent::SPLIT_EVENT) {
            SplitEventSPtr split_event = std::dynamic_pointer_cast<SplitEvent>(event);
            node = split_event->getNode();
        } else if (event->getType() == AbstractEvent::TRIANGLE_EVENT) {
            TriangleEventSPtr triangle_event = std::dynamic_pointer_cast<TriangleEvent>(event);
            node = triangle_event->getNode();
        }
        if (node) {
            result = removeNode(node);
        }
    }
    return result;
}

void StraightSkeleton::addNode(NodeSPtr node) {
    std::list<NodeSPtr>::iterator it = nodes_.insert(nodes_.end(), node);
    node->setSkel(shared_from_this());
    node->setListIt(it);
}

bool StraightSkeleton::removeNode(NodeSPtr node) {
    bool result = false;
    if (node->getSkel() == shared_from_this()) {
        nodes_.erase(node->getListIt());
        node->setSkel(StraightSkeletonSPtr());
        node->setListIt(std::list<NodeSPtr>::iterator());
        result = true;
    }
    return result;
}

void StraightSkeleton::addArc(ArcSPtr arc) {
    std::list<ArcSPtr>::iterator it = arcs_.insert(arcs_.end(), arc);
    arc->setSkel(shared_from_this());
    arc->setListIt(it);
}

bool StraightSkeleton::removeArc(ArcSPtr arc) {
    bool result = false;
    if (arc->getSkel() == shared_from_this()) {
        arcs_.erase(arc->getListIt());
        arc->setSkel(StraightSkeletonSPtr());
        arc->setListIt(std::list<ArcSPtr>::iterator());
        result = true;
    }
    return result;
}

SharedMutex& StraightSkeleton::mutex() {
    return this->mutex_;
}

std::list<AbstractEventSPtr>& StraightSkeleton::events() {
    return this->events_;
}

std::list<NodeSPtr>& StraightSkeleton::nodes() {
    return this->nodes_;
}

std::list<ArcSPtr>& StraightSkeleton::arcs() {
    return this->arcs_;
}

PolygonSPtr StraightSkeleton::getPolygon() const {
    return this->polygon_;
}

void StraightSkeleton::setPolygon(PolygonSPtr polygon) {
    this->polygon_ = polygon;
}

int StraightSkeleton::getID() const {
    return this->id_;
}

void StraightSkeleton::setID(int id) {
    this->id_ = id;
}

void StraightSkeleton::resetAllIDs() {
    std::list<AbstractEventSPtr>::iterator it_e = events_.begin();
    while (it_e != events_.end()) {
        AbstractEventSPtr event = *it_e++;
        event->setID(-1);
    }
    std::list<ArcSPtr>::iterator it_a = arcs_.begin();
    while (it_a != arcs_.end()) {
        ArcSPtr arc = *it_a++;
        arc->setID(-1);
    }
    std::list<NodeSPtr>::iterator it_n = nodes_.begin();
    while (it_n != nodes_.end()) {
        NodeSPtr node = *it_n++;
        node->setID(-1);
    }
    setID(-1);
}

bool StraightSkeleton::isConsistent() const {
    bool result = true;

    std::list<NodeSPtr>::const_iterator it_n = nodes_.begin();
    while (it_n != nodes_.end()) {
        NodeSPtr node = *it_n++;
        if (node->getSkel() != shared_from_this()) {
            DEBUG_VAR(node->toString());
            result = false;
            break;
        }
        std::list<ArcWPtr>::const_iterator it_a = node->arcs().begin();
        while (it_a != node->arcs().end()) {
            ArcWPtr arc_wptr = *it_a++;
            if (arc_wptr.expired()) {
                DEBUG_VAR(node->toString());
            } else {
                ArcSPtr arc = ArcSPtr(arc_wptr);
                if (node != arc->getNodeSrc() && node != arc->getNodeDst()) {
                    DEBUG_VAR(node->toString());
                    DEBUG_VAR(arc->toString());
                    result = false;
                    break;
                }
            }
        }
    }

    std::list<ArcSPtr>::const_iterator it_a = arcs_.begin();
    while (it_a != arcs_.end()) {
        ArcSPtr arc = *it_a++;
        ArcWPtr arc_wptr(arc);
        if (arc->getSkel() != shared_from_this()) {
            DEBUG_VAR(arc->toString());
            result = false;
            break;
        }
        std::list<ArcWPtr> warcs = arc->getNodeSrc()->arcs();
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

int StraightSkeleton::countEvents(int type) const {
    int result = 0;
    std::list<AbstractEventSPtr>::const_iterator it_e = events_.begin();
    while (it_e != events_.end()) {
        AbstractEventSPtr event = *it_e++;
        if (event->getType() == type) {
            result += 1;
        }
    }
    return result;
}

std::string StraightSkeleton::getDescription() const {
    return this->description_;
}

void StraightSkeleton::setDescription(const std::string& desc) {
    this->description_ = desc;
}

void StraightSkeleton::appendDescription(const std::string& desc) {
    this->description_.append(desc);
}

std::string StraightSkeleton::toString() const {
    std::stringstream sstr;
    sstr << "StraightSkeleton(";
    if (id_ != -1) {
        sstr << "id=" << util::StringFactory::fromInteger(id_);
    } else {
        sstr << util::StringFactory::fromPointer(this);
    }
    sstr << "," << std::endl;
    sstr << "Nodes:  " << nodes_.size() << std::endl;
    sstr << "Arcs:   " << arcs_.size() << std::endl;
    sstr << "Events: " << events_.size() << std::endl;
    sstr << "  ConstOffsetEvents: " << countEvents(AbstractEvent::CONST_OFFSET_EVENT) << std::endl;
    sstr << "  EdgeEvents:        " << countEvents(AbstractEvent::EDGE_EVENT) << std::endl;
    sstr << "  SplitEvents:       " << countEvents(AbstractEvent::SPLIT_EVENT) << std::endl;
    sstr << "  TriangleEvents:    " << countEvents(AbstractEvent::TRIANGLE_EVENT) << std::endl;
    sstr << ")" << std::endl;
    return sstr.str();
}

} } }
