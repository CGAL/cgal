// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * @file   data/3d/skel/StraightSkeleton.h
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_SDS_STRAIGHT_SKELETON_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_SDS_STRAIGHT_SKELETON_H

#include <CGAL/Straight_skeleton_3/internal/SDS/Node.h>
#include <CGAL/Straight_skeleton_3/internal/SDS/Arc.h>
#include <CGAL/Straight_skeleton_3/internal/SDS/Sheet.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/events/Abstract_event.h>
#include <CGAL/Straight_skeleton_3/internal/weak_find.h>
#include <CGAL/Straight_skeleton_3/IO/String_factory.h>

#include <list>
#include <memory>
#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace SDS {

template <typename Traits>
class StraightSkeleton
  : public std::enable_shared_from_this<StraightSkeleton<Traits> >
{
  using Polyhedron = HDS::Polyhedron<Traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

public:
  using Node = SDS::Node<Traits>;
  using NodeWPtr = std::weak_ptr<Node>;
  using NodeSPtr = std::shared_ptr<Node>;
  using Arc = SDS::Arc<Traits>;
  using ArcWPtr = std::weak_ptr<Arc>;
  using ArcSPtr = std::shared_ptr<Arc>;
  using Sheet = SDS::Sheet<Traits>;
  using SheetWPtr = std::weak_ptr<Sheet>;
  using SheetSPtr = std::shared_ptr<Sheet>;

  using StraightSkeletonWPtr = std::weak_ptr<StraightSkeleton<Traits> >;
  using StraightSkeletonSPtr = std::shared_ptr<StraightSkeleton<Traits> >;

private:
  using AbstractEvent = algorithm::AbstractEvent<Traits>;
  using AbstractEventSPtr = std::shared_ptr<AbstractEvent>;

public:
  StraightSkeleton()
    : id_(-1)
  { }

  ~StraightSkeleton()
  {
    events_.clear();
    sheets_.clear();
    arcs_.clear();
    nodes_.clear();
  }

  static StraightSkeletonSPtr create()
  {
    return std::make_shared<StraightSkeleton>();
  }

  /**
  * also adds the node
  */
  void addEvent(const AbstractEventSPtr& event)
  {
    typename std::list<AbstractEventSPtr>::iterator it = events_.insert(events_.end(), event);
    event->setListIt(it);
  }

  bool removeEvent(const AbstractEventSPtr& event)
  {
    bool result = false;
    events_.erase(event->getListIt());
    event->setListIt(typename std::list<AbstractEventSPtr>::iterator());
    return result;
  }

  void addNode(NodeSPtr node)
  {
    typename std::list<NodeSPtr>::iterator it = nodes_.insert(nodes_.end(), node);
    node->setSkel(this->shared_from_this());
    node->setListIt(it);
  }

  bool removeNode(NodeSPtr node)
  {
    bool result = false;
    if (node->getSkel() == this->shared_from_this()) {
      nodes_.erase(node->getListIt());
      node->setSkel(StraightSkeletonSPtr());
      node->setListIt(std::list<NodeSPtr>::iterator());
      result = true;
    }
    return result;
  }

  void addArc(ArcSPtr arc)
  {
    typename std::list<ArcSPtr>::iterator it = arcs_.insert(arcs_.end(), arc);
    arc->setSkel(this->shared_from_this());
    arc->setListIt(it);
  }

  bool removeArc(ArcSPtr arc)
  {
    bool result = false;
    if (arc->getSkel() == this->shared_from_this()) {
      arcs_.erase(arc->getListIt());
      arc->setSkel(StraightSkeletonSPtr());
      arc->setListIt(typename std::list<ArcSPtr>::iterator());
      result = true;
    }
    return result;
  }

  void addSheet(SheetSPtr sheet)
  {
    typename std::list<SheetSPtr>::iterator it = sheets_.insert(sheets_.end(), sheet);
    sheet->setSkel(this->shared_from_this());
    sheet->setListIt(it);
  }

  bool removeSheet(SheetSPtr sheet)
  {
    bool result = false;
    if (sheet->getSkel() == this->shared_from_this()) {
      sheets_.erase(sheet->getListIt());
      sheet->setSkel(StraightSkeletonSPtr());
      sheet->setListIt(typename std::list<SheetSPtr>::iterator());
      result = true;
    }
    return result;
  }

  std::list<AbstractEventSPtr>& events()
  {
    return this->events_;
  }

  std::list<NodeSPtr>& nodes()
  {
    return this->nodes_;
  }

  std::list<ArcSPtr>& arcs()
  {
    return this->arcs_;
  }

  std::list<SheetSPtr>& sheets()
  {
    return this->sheets_;
  }

  PolyhedronSPtr getPolyhedron() const
  {
    return this->polyhedron_;
  }

  void setPolyhedron(PolyhedronSPtr polyhedron)
  {
    this->polyhedron_ = polyhedron;
  }

  int getID() const
  {
    return this->id_;
  }

  void setID(int id)
  {
    this->id_ = id;
  }

  void resetAllIDs()
  {
    typename std::list<AbstractEventSPtr>::iterator it_e = events_.begin();
    while (it_e != events_.end()) {
      AbstractEventSPtr event = *it_e++;
      event->setID(-1);
    }
    typename std::list<SheetSPtr>::iterator it_s = sheets_.begin();
    while (it_s != sheets_.end()) {
      SheetSPtr sheet = *it_s++;
      sheet->setID(-1);
    }
    typename std::list<ArcSPtr>::iterator it_a = arcs_.begin();
    while (it_a != arcs_.end()) {
      ArcSPtr arc = *it_a++;
      arc->setID(-1);
    }
    typename std::list<NodeSPtr>::iterator it_n = nodes_.begin();
    while (it_n != nodes_.end()) {
      NodeSPtr node = *it_n++;
      node->setID(-1);
    }
    setID(-1);
  }

  std::string getConfig() const
  {
    return this->config_;
  }

  void setConfig(const std::string& config)
  {
    this->config_ = config;
  }

  void appendConfig(const std::string& config)
  {
    this->config_.append(config);
  }

  bool isConsistent() const
  {
    bool result = true;

    typename std::list<NodeSPtr>::const_iterator it_n = nodes_.begin();
    while (it_n != nodes_.end()) {
      NodeSPtr node = *it_n++;
      if (node->getSkel() != this->shared_from_this()) {
        CGAL_SS3_SKEL_DS_TRACE(node->toString());
        result = false;
        break;
      }
      typename std::list<ArcWPtr>::const_iterator it_a = node->arcs().begin();
      while (it_a != node->arcs().end()) {
        ArcWPtr arc_wptr = *it_a++;
        if (ArcSPtr arc = arc_wptr.lock()) {
          if (node != arc->getNodeSrc() && node != arc->getNodeDst()) {
            CGAL_SS3_SKEL_DS_TRACE(node->toString());
            CGAL_SS3_SKEL_DS_TRACE(arc->toString());
            result = false;
            break;
          }
        } else {
          CGAL_SS3_SKEL_DS_TRACE(node->toString());
        }
      }
      typename std::list<SheetWPtr>::const_iterator it_s = node->sheets().begin();
      while (it_s != node->sheets().end()) {
        SheetWPtr sheet_wptr = *it_s++;
        if (SheetSPtr sheet = sheet_wptr.lock()) {
          std::list<NodeSPtr> nodes = sheet->nodes();
          if (nodes.end() == std::find(nodes.begin(), nodes.end(), node)) {
            CGAL_SS3_SKEL_DS_TRACE(node->toString());
            CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
            result = false;
            break;
          }
        } else {
          CGAL_SS3_SKEL_DS_TRACE(node->toString());
        }
      }
    }

    typename std::list<ArcSPtr>::const_iterator it_a = arcs_.begin();
    while (it_a != arcs_.end()) {
      ArcSPtr arc = *it_a++;
      ArcWPtr arc_wptr(arc);
      if (arc->getSkel() != this->shared_from_this()) {
        CGAL_SS3_SKEL_DS_TRACE(arc->toString());
        result = false;
        break;
      }
      std::list<ArcWPtr> warcs = arc->getNodeSrc()->arcs();
      if (warcs.end() == STL_Extension::internal::weak_find(warcs.begin(), warcs.end(), arc_wptr)) {
        CGAL_SS3_SKEL_DS_TRACE(arc->toString());
        CGAL_SS3_SKEL_DS_TRACE(arc->getNodeSrc()->toString());
        result = false;
        break;
      }
      if (arc->hasNodeDst()) {
        warcs = arc->getNodeDst()->arcs();
        if (warcs.end() == STL_Extension::internal::weak_find(warcs.begin(), warcs.end(), arc_wptr)) {
          CGAL_SS3_SKEL_DS_TRACE(arc->toString());
          CGAL_SS3_SKEL_DS_TRACE(arc->getNodeDst()->toString());
          result = false;
          break;
        }
      }
      unsigned int num_sheets = 0;
      typename std::list<SheetWPtr>::const_iterator it_s = arc->sheets().begin();
      while (it_s != arc->sheets().end()) {
        SheetWPtr sheet_wptr = *it_s++;
        if (SheetSPtr sheet = sheet_wptr.lock()) {
          ++num_sheets;
          std::list<ArcSPtr> arcs = sheet->arcs();
          if (arcs.end() == std::find(arcs.begin(), arcs.end(), arc)) {
            CGAL_SS3_SKEL_DS_TRACE(arc->toString());
            CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
            result = false;
            break;
          }
        } else {
          CGAL_SS3_SKEL_DS_TRACE(arc->toString());
        }
      }
      if (num_sheets != 3) {
        CGAL_SS3_SKEL_DS_TRACE("Warning: Arc does not have 3 sheets.");
        CGAL_SS3_SKEL_DS_TRACE(arc->toString());
        CGAL_SS3_SKEL_DS_TRACE(num_sheets);
      }
    }

    typename std::list<SheetSPtr>::const_iterator it_s = sheets_.begin();
    while (it_s != sheets_.end()) {
      SheetSPtr sheet = *it_s++;
      SheetWPtr sheet_wptr(sheet);
      if (sheet->getSkel() != this->shared_from_this()) {
        CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
        result = false;
        break;
      }
      typename std::list<NodeSPtr>::const_iterator it_n = sheet->nodes().begin();
      while (it_n != sheet->nodes().end()) {
        NodeSPtr node = *it_n++;
        std::list<SheetWPtr> wsheets = node->sheets();
        if (wsheets.end() == STL_Extension::internal::weak_find(wsheets.begin(), wsheets.end(), sheet_wptr)) {
          CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
          CGAL_SS3_SKEL_DS_TRACE(node->toString());
          result = false;
          break;
        }
      }
      typename std::list<ArcSPtr>::const_iterator it_a = sheet->arcs().begin();
      while (it_a != sheet->arcs().end()) {
        ArcSPtr arc = *it_a++;
        std::list<SheetWPtr> wsheets = arc->sheets();
        if (wsheets.end() == STL_Extension::internal::weak_find(wsheets.begin(), wsheets.end(), sheet_wptr)) {
          CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
          CGAL_SS3_SKEL_DS_TRACE(arc->toString());
          result = false;
          break;
        }
      }
    }

    return result;
  }

  int countEvents(int type) const
  {
    int result = 0;
    typename std::list<AbstractEventSPtr>::const_iterator it_e = events_.begin();
    while (it_e != events_.end()) {
      AbstractEventSPtr event = *it_e++;
      if (event->getType() == type) {
        result += 1;
      }
    }
    return result;
  }

  std::string toString() const
  {
    std::stringstream sstr;
    sstr << "StraightSkeleton(";
    if (id_ != -1) {
        sstr << "id=" << IO::StringFactory::fromInteger(id_);
    } else {
        sstr << IO::StringFactory::fromPointer(this);
    }
    sstr << "," << std::endl;
    sstr << "Nodes:  " << nodes_.size() << std::endl;
    sstr << "Arcs:   " << arcs_.size() << std::endl;
    sstr << "Sheets: " << sheets_.size() << std::endl;
    sstr << "Events: " << events_.size() << std::endl;
    sstr << "    ConstOffsetEvents:     " << countEvents(AbstractEvent::CONST_OFFSET_EVENT) << std::endl;
    sstr << "    SaveOffsetEvents:      " << countEvents(AbstractEvent::SAVE_OFFSET_EVENT) << std::endl;
    sstr << "  VanishEvents:" << std::endl;
    sstr << "    Generic VanishEvents:  " << countEvents(AbstractEvent::VANISH_EVENT) << std::endl;
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

protected:
  PolyhedronSPtr polyhedron_;

  std::list<AbstractEventSPtr> events_;
  std::list<NodeSPtr> nodes_;
  std::list<ArcSPtr> arcs_;
  std::list<SheetSPtr> sheets_;

  int id_;
  std::string config_;
};

} // namespace SDS
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_SDS_STRAIGHT_SKELETON_H */
