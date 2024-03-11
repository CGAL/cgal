/**
 * @file   data/3d/skel/StraightSkeleton.cpp
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#include "data/3d/skel/StraightSkeleton.h"

#include "debug.h"
#include "data/3d/skel/AbstractEvent.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/Arc.h"
#include "data/3d/skel/Sheet.h"
#include "util/StringFactory.h"
#include <sstream>

namespace data { namespace _3d { namespace skel {

StraightSkeleton::StraightSkeleton() {
    id_ = -1;
}

StraightSkeleton::~StraightSkeleton() {
    events_.clear();
    sheets_.clear();
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
}

bool StraightSkeleton::removeEvent(AbstractEventSPtr event) {
    bool result = false;
    if (event->getSkel() == shared_from_this()) {
        events_.erase(event->getListIt());
        event->setSkel(StraightSkeletonSPtr());
        event->setListIt(std::list<AbstractEventSPtr>::iterator());
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

void StraightSkeleton::addSheet(SheetSPtr sheet) {
    std::list<SheetSPtr>::iterator it = sheets_.insert(sheets_.end(), sheet);
    sheet->setSkel(shared_from_this());
    sheet->setListIt(it);
}

bool StraightSkeleton::removeSheet(SheetSPtr sheet) {
    bool result = false;
    if (sheet->getSkel() == shared_from_this()) {
        sheets_.erase(sheet->getListIt());
        sheet->setSkel(StraightSkeletonSPtr());
        sheet->setListIt(std::list<SheetSPtr>::iterator());
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

std::list<SheetSPtr>& StraightSkeleton::sheets() {
    return this->sheets_;
}

PolyhedronSPtr StraightSkeleton::getPolyhedron() const {
    return this->polyhedron_;
}

void StraightSkeleton::setPolyhedron(PolyhedronSPtr polyhedron) {
    this->polyhedron_ = polyhedron;
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
    std::list<SheetSPtr>::iterator it_s = sheets_.begin();
    while (it_s != sheets_.end()) {
        SheetSPtr sheet = *it_s++;
        sheet->setID(-1);
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

std::string StraightSkeleton::getConfig() const {
    return this->config_;
}

void StraightSkeleton::setConfig(const std::string& config) {
    this->config_ = config;
}

void StraightSkeleton::appendConfig(const std::string& config) {
    this->config_.append(config);
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
        std::list<SheetWPtr>::const_iterator it_s = node->sheets().begin();
        while (it_s != node->sheets().end()) {
            SheetWPtr sheet_wptr = *it_s++;
            if (sheet_wptr.expired()) {
                DEBUG_VAR(node->toString());
            } else {
                SheetSPtr sheet = SheetSPtr(sheet_wptr);
                std::list<NodeSPtr> nodes = sheet->nodes();
                if (nodes.end() ==
                        std::find(nodes.begin(), nodes.end(), node)) {
                    DEBUG_VAR(node->toString());
                    DEBUG_VAR(sheet->toString());
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
        unsigned int num_sheets = 0;
        std::list<SheetWPtr>::const_iterator it_s = arc->sheets().begin();
        while (it_s != arc->sheets().end()) {
            SheetWPtr sheet_wptr = *it_s++;
            if (sheet_wptr.expired()) {
                DEBUG_VAR(arc->toString());
            } else {
                SheetSPtr sheet = SheetSPtr(sheet_wptr);
                num_sheets++;
                std::list<ArcSPtr> arcs = sheet->arcs();
                if (arcs.end() ==
                        std::find(arcs.begin(), arcs.end(), arc)) {
                    DEBUG_VAR(arc->toString());
                    DEBUG_VAR(sheet->toString());
                    result = false;
                    break;
                }
            }
        }
        if (num_sheets != 3) {
            DEBUG_VAL("Warning: Arc does not have 3 sheets.");
            DEBUG_VAR(arc->toString());
            DEBUG_VAR(num_sheets);
        }
    }

    std::list<SheetSPtr>::const_iterator it_s = sheets_.begin();
    while (it_s != sheets_.end()) {
        SheetSPtr sheet = *it_s++;
        SheetWPtr sheet_wptr(sheet);
        if (sheet->getSkel() != shared_from_this()) {
            DEBUG_VAR(sheet->toString());
            result = false;
            break;
        }
        std::list<NodeSPtr>::const_iterator it_n = sheet->nodes().begin();
        while (it_n != sheet->nodes().end()) {
            NodeSPtr node = *it_n++;
            std::list<SheetWPtr> wsheets = node->sheets();
            if (wsheets.end() == std::find(
                    wsheets.begin(), wsheets.end(), sheet_wptr)) {
                DEBUG_VAR(sheet->toString());
                DEBUG_VAR(node->toString());
                result = false;
                break;
            }
        }
        std::list<ArcSPtr>::const_iterator it_a = sheet->arcs().begin();
        while (it_a != sheet->arcs().end()) {
            ArcSPtr arc = *it_a++;
            std::list<SheetWPtr> wsheets = arc->sheets();
            if (wsheets.end() == std::find(
                    wsheets.begin(), wsheets.end(), sheet_wptr)) {
                DEBUG_VAR(sheet->toString());
                DEBUG_VAR(arc->toString());
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
    sstr << "Sheets: " << sheets_.size() << std::endl;
    sstr << "Events: " << events_.size() << std::endl;
    sstr << "    ConstOffsetEvents:     " << countEvents(AbstractEvent::CONST_OFFSET_EVENT) << std::endl;
    sstr << "    SaveOffsetEvents:      " << countEvents(AbstractEvent::SAVE_OFFSET_EVENT) << std::endl;
    sstr << "  VanishEvents:" << std::endl;
    sstr << "    EdgeEvents:            " << countEvents(AbstractEvent::EDGE_EVENT) << std::endl;
    sstr << "    EdgeMergeEvents:       " << countEvents(AbstractEvent::EDGE_MERGE_EVENT) << std::endl;
    sstr << "    TriangleEvents:        " << countEvents(AbstractEvent::TRIANGLE_EVENT) << std::endl;
    sstr << "    DblEdgeMergeEvents:    " << countEvents(AbstractEvent::DBL_EDGE_MERGE_EVENT) << std::endl;
    sstr << "    DblTriangleEvents:     " << countEvents(AbstractEvent::DBL_TRIANGLE_EVENT) << std::endl;
    sstr << "    TetrahedronEvents:     " << countEvents(AbstractEvent::TETRAHEDRON_EVENT) << std::endl;
    sstr << "  ContactEvents:" << std::endl;
    sstr << "    VertexEvents:          " << countEvents(AbstractEvent::VERTEX_EVENT) << std::endl;
    sstr << "    FlipVertexEvents:      " << countEvents(AbstractEvent::FLIP_VERTEX_EVENT) << std::endl;
    sstr << "    SurfaceEvents:         " << countEvents(AbstractEvent::SURFACE_EVENT) << std::endl;
    sstr << "    PolyhedronSplitEvents: " << countEvents(AbstractEvent::POLYHEDRON_SPLIT_EVENT) << std::endl;
    sstr << "    SplitMergeEvents:      " << countEvents(AbstractEvent::SPLIT_MERGE_EVENT) << std::endl;
    sstr << "    EdgeSplitEvents:       " << countEvents(AbstractEvent::EDGE_SPLIT_EVENT) << std::endl;
    sstr << "    PierceEvents:          " << countEvents(AbstractEvent::PIERCE_EVENT) << std::endl;
    sstr << ")" << std::endl;
    return sstr.str();
}

} } }
