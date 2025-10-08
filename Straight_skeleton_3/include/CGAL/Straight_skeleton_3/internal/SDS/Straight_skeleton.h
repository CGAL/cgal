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
 * file   data/3d/skel/StraightSkeleton.h
 * author Gernot Walzl
 * date   2011-11-26
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
#include <unordered_set>

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
  ~StraightSkeleton()
  {
    polyhedron_.reset();
    sheets_.clear();
    arcs_.clear();
    nodes_.clear();
  }

  static StraightSkeletonSPtr create()
  {
    return std::make_shared<StraightSkeleton>();
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

  void mergeSheets(SheetSPtr sheet_into,
                   SheetSPtr sheet_from)
  {
    CGAL_precondition(sheet_into && sheet_from);

    sheet_into->merge(sheet_from);

    for (NodeSPtr node : sheet_from->nodes()) {
      node->removeSheet(sheet_from);
    }
    for (ArcSPtr arc : sheet_from->arcs()) {
      arc->removeSheet(sheet_from);
    }
    removeSheet(sheet_from);
  }

  PolyhedronSPtr getPolyhedron() const {
      return this->polyhedron_;
  }

  void setPolyhedron(PolyhedronSPtr polyhedron) {
      this->polyhedron_ = polyhedron;
  }

  std::list<SheetSPtr>& sheets()
  {
    return this->sheets_;
  }

  std::list<ArcSPtr>& arcs()
  {
    return this->arcs_;
  }

  std::list<NodeSPtr>& nodes()
  {
    return this->nodes_;
  }

  void resetAllIDs()
  {
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
  }

  bool isConsistent(bool is_partial = true) const
  {
    bool result = true;

    // Check for duplicate elements
    {
      std::unordered_set<NodeSPtr> node_set;
      for (const NodeSPtr& node : nodes_) {
        if (!node_set.insert(node).second) {
          CGAL_SS3_SKEL_DS_TRACE("Error: duplicate node in nodes_ list");
          CGAL_SS3_SKEL_DS_TRACE(node->toString());
          result = false;
          break;
        }
      }
      std::unordered_set<ArcSPtr> arc_set;
      for (const ArcSPtr& arc : arcs_) {
        if (!arc_set.insert(arc).second) {
          CGAL_SS3_SKEL_DS_TRACE("Error: duplicate arc in arcs_ list");
          CGAL_SS3_SKEL_DS_TRACE(arc->toString());
          result = false;
          break;
        }
      }
      std::unordered_set<SheetSPtr> sheet_set;
      for (const SheetSPtr& sheet : sheets_) {
        if (!sheet_set.insert(sheet).second) {
          CGAL_SS3_SKEL_DS_TRACE("Error: duplicate sheet in sheets_ list");
          CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
          result = false;
          break;
        }
      }
    }

    // Check uniqueness of Node IDs
    {
      std::unordered_set<int> node_ids;
      for (const NodeSPtr& node : nodes_) {
        int id = node->getID();
        if (!node_ids.insert(id).second) {
          CGAL_SS3_SKEL_DS_TRACE("Error: duplicate Node ID in nodes_ list: " << id);
          CGAL_SS3_SKEL_DS_TRACE(node->toString());
          result = false;
          break;
        }
      }
    }
    // Check uniqueness of Arc IDs
    {
      std::unordered_set<int> arc_ids;
      for (const ArcSPtr& arc : arcs_) {
        int id = arc->getID();
        if (!arc_ids.insert(id).second) {
          CGAL_SS3_SKEL_DS_TRACE("Error: duplicate Arc ID in arcs_ list: " << id);
          CGAL_SS3_SKEL_DS_TRACE(arc->toString());
          result = false;
          break;
        }
      }
    }
    // Check uniqueness of Sheet IDs
    {
      std::unordered_set<int> sheet_ids;
      for (const SheetSPtr& sheet : sheets_) {
        int id = sheet->getID();
        if (!sheet_ids.insert(id).second) {
          CGAL_SS3_SKEL_DS_TRACE("Error: duplicate Sheet ID in sheets_ list: " << id);
          CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
          result = false;
          break;
        }
      }
    }

    typename std::list<NodeSPtr>::const_iterator it_n = nodes_.begin();
    while (it_n != nodes_.end()) {
      NodeSPtr node = *it_n++;
      CGAL_SS3_DEBUG_SPTR(node);
      if (node->getSkel() != this->shared_from_this()) {
        CGAL_SS3_SKEL_DS_TRACE(node->toString());
        result = false;
        break;
      }
      typename std::list<ArcWPtr>::const_iterator it_a = node->arcs().begin();
      while (it_a != node->arcs().end()) {
        ArcWPtr arc_wptr = *it_a++;
        CGAL_SS3_DEBUG_WPTR(arc_wptr);
        ArcSPtr arc = arc_wptr.lock();
        if (arc_wptr.expired() || !arc) {
          CGAL_SS3_SKEL_DS_TRACE("Error: invalid ArcWPtr (expired or cannot lock) in node's arc list");
          CGAL_SS3_SKEL_DS_TRACE(node->toString());
          result = false;
          break;
        }

        if (std::find(arcs_.begin(), arcs_.end(), arc) == arcs_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: node references an arc not present in the skeleton's arc list");
          CGAL_SS3_SKEL_DS_TRACE(node->toString());
          CGAL_SS3_SKEL_DS_TRACE(arc->toString());
          result = false;
          break;
        }

        if (node != arc->getNodeSrc() && node != arc->getNodeDst()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: incident arc does not have node as endpoint");
          CGAL_SS3_SKEL_DS_TRACE(node->toString());
          CGAL_SS3_SKEL_DS_TRACE(arc->toString());
          result = false;
          break;
        }
      }
      typename std::list<SheetWPtr>::const_iterator it_s = node->sheets().begin();
      while (it_s != node->sheets().end()) {
        SheetWPtr sheet_wptr = *it_s++;
        CGAL_SS3_DEBUG_WPTR(sheet_wptr);
        SheetSPtr sheet = sheet_wptr.lock();
        if (sheet_wptr.expired() || !sheet) {
          CGAL_SS3_SKEL_DS_TRACE("Error: invalid SheetWPtr (expired or cannot lock) in node's sheet list");
          CGAL_SS3_SKEL_DS_TRACE(node->toString());
          result = false;
          break;
        }

        if (std::find(sheets_.begin(), sheets_.end(), sheet) == sheets_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: node references a sheet not present in the skeleton's sheet list");
          CGAL_SS3_SKEL_DS_TRACE(node->toString());
          CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
          result = false;
          break;
        }

        std::list<NodeSPtr> nodes = sheet->nodes();
        if (nodes.end() == std::find(nodes.begin(), nodes.end(), node)) {
          CGAL_SS3_SKEL_DS_TRACE("Error: incident sheet does not have node");
          CGAL_SS3_SKEL_DS_TRACE(node->toString());
          CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
          result = false;
          break;
        }
      }
    }

    typename std::list<ArcSPtr>::const_iterator it_a = arcs_.begin();
    while (it_a != arcs_.end()) {
      ArcSPtr arc = *it_a++;
      CGAL_SS3_DEBUG_SPTR(arc);
      ArcWPtr arc_wptr(arc);

      if (arc->getSkel() != this->shared_from_this()) {
        CGAL_SS3_SKEL_DS_TRACE(arc->toString());
        result = false;
        break;
      }

      if (!arc->getNodeSrc()) {
        CGAL_SS3_SKEL_DS_TRACE("Error: arc does not have source");
        result = false;
        break;
      }

      if (std::find(nodes_.begin(), nodes_.end(), arc->getNodeSrc()) == nodes_.end()) {
        CGAL_SS3_SKEL_DS_TRACE("Error: arc source is not a node present in the skeleton's node list");
        CGAL_SS3_SKEL_DS_TRACE(arc->getNodeSrc()->toString());
        result = false;
        break;
      }

      std::list<ArcWPtr> warcs = arc->getNodeSrc()->arcs();
      if (warcs.end() == STL_Extension::internal::weak_find(warcs.begin(), warcs.end(), arc_wptr)) {
        CGAL_SS3_SKEL_DS_TRACE("Error: arc's source node does not have arc");
        CGAL_SS3_SKEL_DS_TRACE(arc->toString());
        CGAL_SS3_SKEL_DS_TRACE(arc->getNodeSrc()->toString());
        result = false;
        break;
      }

      if (arc->hasNodeDst()) {
        if (std::find(nodes_.begin(), nodes_.end(), arc->getNodeDst()) == nodes_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: arc destination is not a node present in the skeleton's node list");
          CGAL_SS3_SKEL_DS_TRACE(arc->getNodeDst()->toString());
          result = false;
          break;
        }

        warcs = arc->getNodeDst()->arcs();
        if (warcs.end() == STL_Extension::internal::weak_find(warcs.begin(), warcs.end(), arc_wptr)) {
          CGAL_SS3_SKEL_DS_TRACE("Error: arc's target node does not have arc");
          CGAL_SS3_SKEL_DS_TRACE(arc->toString());
          CGAL_SS3_SKEL_DS_TRACE(arc->getNodeDst()->toString());
          result = false;
          break;
        }
      } else {
        if (is_partial) {
          if (!arc->getDirection()) {
            CGAL_SS3_SKEL_DS_TRACE("Error: arc does not have destination nor direction (partial skeleton)");
            result = false;
            break;
          }
        } else {
          CGAL_SS3_SKEL_DS_TRACE("Error: arc does not have destination (full skeleton)");
          result = false;
          break;
        }
      }
      unsigned int num_sheets = 0;
      typename std::list<SheetWPtr>::const_iterator it_s = arc->sheets().begin();
      while (it_s != arc->sheets().end()) {
        SheetWPtr sheet_wptr = *it_s++;
        CGAL_SS3_DEBUG_WPTR(sheet_wptr);
        SheetSPtr sheet = sheet_wptr.lock();
        if (sheet_wptr.expired() || !sheet) {
          CGAL_SS3_SKEL_DS_TRACE("Error: invalid SheetWPtr (expired or cannot lock) in arc's sheet list");
          CGAL_SS3_SKEL_DS_TRACE(arc->toString());
          result = false;
          break;
        }

        if (std::find(sheets_.begin(), sheets_.end(), sheet) == sheets_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: arc references a sheet that is not present in the skeleton's sheet list");
          CGAL_SS3_SKEL_DS_TRACE(arc->toString());
          CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
          result = false;
          break;
        }

        ++num_sheets;

        std::list<ArcSPtr> arcs = sheet->arcs();
        if (arcs.end() == std::find(arcs.begin(), arcs.end(), arc)) {
          CGAL_SS3_SKEL_DS_TRACE("Error: arc's sheet does not have arc");
          CGAL_SS3_SKEL_DS_TRACE(arc->toString());
          CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
          result = false;
          break;
        }
      }

      if (num_sheets != 3) {
        CGAL_SS3_SKEL_DS_TRACE("Error: Arc does not have 3 sheets.");
        CGAL_SS3_SKEL_DS_TRACE(arc->toString());
        CGAL_SS3_SKEL_DS_TRACE("num_sheets = " << num_sheets);
        result = false;
        break;
      }
    }

    typename std::list<SheetSPtr>::const_iterator it_s = sheets_.begin();
    while (it_s != sheets_.end()) {
      SheetSPtr sheet = *it_s++;
      CGAL_SS3_DEBUG_SPTR(sheet);
      SheetWPtr sheet_wptr(sheet);
      if (sheet->getSkel() != this->shared_from_this()) {
        CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
        result = false;
        break;
      }
      typename std::list<NodeSPtr>::const_iterator it_n = sheet->nodes().begin();
      while (it_n != sheet->nodes().end()) {
        NodeSPtr node = *it_n++;
        CGAL_SS3_DEBUG_SPTR(node);

        if (std::find(nodes_.begin(), nodes_.end(), node) == nodes_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: sheet references a node that is not present in the skeleton's node list");
          CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
          CGAL_SS3_SKEL_DS_TRACE(node->toString());
          result = false;
          break;
        }

        std::list<SheetWPtr> wsheets = node->sheets();
        if (wsheets.end() == STL_Extension::internal::weak_find(wsheets.begin(), wsheets.end(), sheet_wptr)) {
          CGAL_SS3_SKEL_DS_TRACE("Error: sheet's node does not have sheet");
          CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
          CGAL_SS3_SKEL_DS_TRACE(node->toString());
          result = false;
          break;
        }
      }
      typename std::list<ArcSPtr>::const_iterator it_a = sheet->arcs().begin();
      while (it_a != sheet->arcs().end()) {
        ArcSPtr arc = *it_a++;
        CGAL_SS3_DEBUG_SPTR(arc);

        if (std::find(arcs_.begin(), arcs_.end(), arc) == arcs_.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: sheet references an arc that is not present in the skeleton's arc list");
          CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
          CGAL_SS3_SKEL_DS_TRACE(arc->toString());
          result = false;
          break;
        }

        std::list<SheetWPtr> wsheets = arc->sheets();
        if (wsheets.end() == STL_Extension::internal::weak_find(wsheets.begin(), wsheets.end(), sheet_wptr)) {
          CGAL_SS3_SKEL_DS_TRACE("Error: sheet's arc does not have sheet");
          CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
          CGAL_SS3_SKEL_DS_TRACE(arc->toString());
          result = false;
          break;
        }

        std::vector<NodeSPtr> arc_nodes = { arc->getNodeSrc() };
        if (arc->hasNodeDst()) {
          arc_nodes.push_back(arc->getNodeDst());
        }
        for (NodeSPtr node : arc_nodes) {
          CGAL_SS3_DEBUG_SPTR(node);
          std::list<NodeSPtr> nodes = sheet->nodes();
          if (nodes.end() == std::find(nodes.begin(), nodes.end(), node)) {
            CGAL_SS3_SKEL_DS_TRACE("Error: Nodes of arc of sheet should be nodes of sheet");
            CGAL_SS3_SKEL_DS_TRACE(node->toString());
            CGAL_SS3_SKEL_DS_TRACE(arc->toString());
            CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
            result = false;
            break;
          }
        }
      }

      // Contour checks
      for (const ArcSPtr& contour : sheet->contours()) {
        // Check that contour has exactly one incident sheet (the current sheet)
        int sheet_count = 0;
        for (const SheetWPtr& contour_sheet_wptr : contour->sheets()) {
          SheetSPtr contour_sheet = contour_sheet_wptr.lock();
          if (contour_sheet == sheet) {
            ++sheet_count;
          }
        }
        if (sheet_count != 1) {
          CGAL_SS3_SKEL_DS_TRACE("Error: contour arc does not have exactly one incident sheet (should be the current sheet)");
          CGAL_SS3_SKEL_DS_TRACE(contour->toString());
          CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
          result = false;
          break;
        }

        // Check that contour has two valid, distinct incident nodes
        NodeSPtr node_src = contour->getNodeSrc();
        NodeSPtr node_dst = contour->getNodeDst();
        if (!node_src || !node_dst) {
          CGAL_SS3_SKEL_DS_TRACE("Error: contour does not have two valid incident nodes");
          CGAL_SS3_SKEL_DS_TRACE(contour->toString());
          result = false;
          break;
        }

        // Check that both nodes appear in the sheet's node list
        const auto& sheet_nodes = sheet->nodes();
        if (std::find(sheet_nodes.begin(), sheet_nodes.end(), node_src) == sheet_nodes.end() ||
            std::find(sheet_nodes.begin(), sheet_nodes.end(), node_dst) == sheet_nodes.end()) {
          CGAL_SS3_SKEL_DS_TRACE("Error: contour's nodes do not appear in sheet's node list");
          CGAL_SS3_SKEL_DS_TRACE(contour->toString());
          CGAL_SS3_SKEL_DS_TRACE(sheet->toString());
          result = false;
          break;
        }
      }
    }

    // @todo check contours

    return result;
  }

  std::string toString() const
  {
    std::stringstream sstr;
    sstr << "StraightSkeleton\n";
    sstr << "  Nodes:  " << nodes_.size() << std::endl;
    for (const NodeSPtr& node : nodes_) {
      sstr << "    " << node->toString() << "\n";
    }
    sstr << "  Arcs:   " << arcs_.size() << std::endl;
    for (const ArcSPtr& arc : arcs_) {
      sstr << "    " << arc->toString() << "\n";
    }
    sstr << "  Sheets: " << sheets_.size() << std::endl;
    for (const SheetSPtr& sheet : sheets_) {
      sstr << "    " << sheet->toString() << "\n";
    }
    return sstr.str();
  }

protected:
  PolyhedronSPtr polyhedron_;
  std::list<NodeSPtr> nodes_;
  std::list<ArcSPtr> arcs_;
  std::list<SheetSPtr> sheets_;
};

} // namespace SDS
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_SDS_STRAIGHT_SKELETON_H */
