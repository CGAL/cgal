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
 * file   data/3d/skel/Sheet.h
 * author Gernot Walzl
 * date   2012-03-27
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_SDS_SHEET_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_SDS_SHEET_H

#include <CGAL/Straight_skeleton_3/IO/String_factory.h>

#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>

#include <list>
#include <memory>
#include <sstream>
#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace SDS {

template <typename Traits>
class Node;

template <typename Traits>
class Arc;

template <typename Traits>
class StraightSkeleton;

template <typename Traits>
class Sheet
  : public std::enable_shared_from_this<Sheet<Traits> >
{
  using Plane_3 = typename Traits::Plane_3;

  using Plane3SPtr = std::shared_ptr<Plane_3>;

private:
  using NodeSPtr = std::shared_ptr<Node<Traits> >;
  using ArcSPtr = std::shared_ptr<Arc<Traits> >;
  using SheetSPtr = std::shared_ptr<Sheet<Traits> >;
  using StraightSkeletonWPtr = std::weak_ptr<StraightSkeleton<Traits> >;
  using StraightSkeletonSPtr = std::shared_ptr<StraightSkeleton<Traits> >;

  using Polyhedron = HDS::Polyhedron<Traits>;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

public:
  Sheet() {
    id_ = next_id_++;
  }

  ~Sheet()
  {
    facet_b_.reset();
    facet_f_.reset();
  }

  static SheetSPtr create()
  {
    return std::make_shared<Sheet>();
  }

  FacetSPtr getFacetB() const
  {
    CGAL_SS3_DEBUG_SPTR(facet_b_);
    return facet_b_;
  }

  void setFacetB(const FacetSPtr& facet_b)
  {
    CGAL_SS3_DEBUG_SPTR(facet_b);
    facet_b_ = facet_b;
  }

  FacetSPtr getFacetF() const
  {
    CGAL_SS3_DEBUG_SPTR(facet_f_);
    return facet_f_;
  }

  void setFacetF(const FacetSPtr& facet_f)
  {
    CGAL_SS3_DEBUG_SPTR(facet_f);
    facet_f_ = facet_f;
  }

  StraightSkeletonSPtr getSkel() const
  {
    return this->skel_.lock();
  }

  void setSkel(const StraightSkeletonSPtr& skel)
  {
    this->skel_ = skel;
  }

  typename std::list<SheetSPtr>::iterator getListIt() const
  {
    return this->list_it_;
  }

  void setListIt(typename std::list<SheetSPtr>::iterator list_it)
  {
    this->list_it_ = list_it;
  }

  int getID() const
  {
    return this->id_;
  }

  void setID(const int id)
  {
    this->id_ = id;
  }

  Plane3SPtr getPlane() const
  {
    CGAL_SS3_DEBUG_SPTR(this->plane_);
    return this->plane_;
  }

  void setPlane(const Plane3SPtr& plane)
  {
    CGAL_SS3_DEBUG_SPTR(plane);
    this->plane_ = plane;
  }

  void addNode(const NodeSPtr& node)
  {
    CGAL_SS3_DEBUG_SPTR(node);
    CGAL_precondition((std::find(nodes_.begin(), nodes_.end(), node) == nodes_.end()));
    /*std::list<NodeSPtr>::iterator it =*/ nodes_.insert(nodes_.end(), node);
    if (!node->containsSheet(this->shared_from_this())) {
      node->addSheet(this->shared_from_this());
    }
  }

  bool removeNode(const NodeSPtr& node)
  {
    CGAL_SS3_DEBUG_SPTR(node);
    bool result = false;
    nodes_.remove(node);
    result = node->removeSheet(this->shared_from_this());
    return result;
  }

  void addArc(const ArcSPtr& arc)
  {
    CGAL_SS3_DEBUG_SPTR(arc);
    /*std::list<ArcSPtr>::iterator it =*/ arcs_.insert(arcs_.end(), arc);
    arc->addSheet(this->shared_from_this());
  }

  bool removeArc(const ArcSPtr& arc)
  {
    CGAL_SS3_DEBUG_SPTR(arc);
    bool result = false;
    arcs_.remove(arc);
    result = arc->removeSheet(this->shared_from_this());
    return result;
  }

  void addContour(const ArcSPtr& contour)
  {
    CGAL_SS3_DEBUG_SPTR(contour);
    /*std::list<ArcSPtr>::iterator it =*/ contours_.insert(contours_.end(), contour);
    contour->addSheet(this->shared_from_this());
  }

  /**
  * merge 'sheet' into this sheet.
  */
  void merge(const SheetSPtr& sheet)
  {
    CGAL_SS3_DEBUG_SPTR(sheet);

    // @fixme something like this should be true once incident facets are fixed
    // CGAL_precondition(this->facet_b_ == sheet->getFacetB() && this->facet_f_ == sheet->getFacetF());

    // copy all arcs that do not yet exist in the sheet
    for (ArcSPtr arc : sheet->arcs()) {
      if (std::find(arcs_.begin(), arcs_.end(), arc) == arcs_.end()) {
        addArc(arc);
      }
    }

    // copy all contours that do not yet exist in the sheet
    for (ArcSPtr contour : sheet->contours()) {
      if (std::find(contours_.begin(), contours_.end(), contour) == contours_.end()) {
        addContour(contour);
      }
    }

    // copy all nodes that do not yet exist in the sheet
    for (NodeSPtr node : sheet->nodes()) {
      if (std::find(nodes_.begin(), nodes_.end(), node) == nodes_.end()) {
        addNode(node);
      }
    }
  }

  std::list<ArcSPtr>& arcs()
  {
    return this->arcs_;
  }

  std::list<ArcSPtr>& contours()
  {
    return this->contours_;
  }

  std::list<NodeSPtr>& nodes()
  {
    return this->nodes_;
  }

  std::string toString() const
  {
    std::stringstream result;
    result << "Sheet(";
    if (id_ != -1) {
      result << "id=" << IO::StringFactory::fromInteger(id_) << ", ";
    } else {
      result << IO::StringFactory::fromPointer(this) << ", ";
    }

    // node IDs
    result << "nodes={";
    bool first = true;
    for (const NodeSPtr& node : nodes_) {
      if (!first) result << ", ";
      result << IO::StringFactory::fromInteger(node->getID());
      first = false;
    }
    result << "}, ";

    // contour IDs
    result << "contours={";
    first = true;
    for (const ArcSPtr& contour : contours_) {
      if (!first) result << ", ";
      result << IO::StringFactory::fromInteger(contour->getID());
      first = false;
    }
    result << "}, ";

    // arc IDs
    result << "arcs={";
    first = true;
    for (const ArcSPtr& arc : arcs_) {
      if (!first) result << ", ";
      result << IO::StringFactory::fromInteger(arc->getID());
      first = false;
    }
    result << "}";

    result << ")";
    return result.str();
  }

protected:
  FacetSPtr facet_b_;
  FacetSPtr facet_f_;

  std::list<ArcSPtr> arcs_;
  // Squatting the arc type, but these are not really arcs: they represent edges of the original
  // polyhedron and are used to have a closed polygon (with holes) when drawing sheets.
  std::list<ArcSPtr> contours_;
  std::list<NodeSPtr> nodes_;
  StraightSkeletonWPtr skel_;

  typename std::list<SheetSPtr>::iterator list_it_;
  Plane3SPtr plane_;

  int id_;

private:
  static int next_id_;
};

template <typename Traits>
int Sheet<Traits>::next_id_ = 0;

} // namespace SDS
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_SDS_SHEET_H */
