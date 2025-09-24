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
  using FacetSPtr = typename Polyhedron::FacetSPtr;

public:
  Sheet() {
    id_ = -1;
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
    facet_b_ = facet_b;
  }

  FacetSPtr getFacetF() const
  {
    CGAL_SS3_DEBUG_SPTR(facet_f_);
    return facet_f_;
  }

  void setFacetF(const FacetSPtr& facet_f)
  {
    facet_f_ = facet_f;
  }

  StraightSkeletonSPtr getSkel() const
  {
    return this->skel_.lock();
  }

  void setSkel(StraightSkeletonSPtr skel)
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

  void setID(int id)
  {
    this->id_ = id;
  }

  Plane3SPtr getPlane() const
  {
    CGAL_SS3_DEBUG_SPTR(this->plane_);
    return this->plane_;
  }

  void setPlane(Plane3SPtr plane)
  {
    this->plane_ = plane;
  }

  void addNode(NodeSPtr node)
  {
    /*std::list<NodeSPtr>::iterator it =*/ nodes_.insert(nodes_.end(), node);
    if (!node->containsSheet(this->shared_from_this())) {
      node->addSheet(this->shared_from_this());
    }
  }

  bool removeNode(NodeSPtr node)
  {
    bool result = false;
    nodes_.remove(node);
    result = node->removeSheet(this->shared_from_this());
    return result;
  }

  void addArc(ArcSPtr arc)
  {
    /*std::list<ArcSPtr>::iterator it =*/ arcs_.insert(arcs_.end(), arc);
    arc->addSheet(this->shared_from_this());
  }

  bool removeArc(ArcSPtr arc)
  {
    bool result = false;
    arcs_.remove(arc);
    result = arc->removeSheet(this->shared_from_this());
    return result;
  }

  std::list<ArcSPtr>& arcs()
  {
    return this->arcs_;
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
        result << "id=" + IO::StringFactory::fromInteger(id_) << ", ";
    } else {
        result << IO::StringFactory::fromPointer(this) << ", ";
    }
    result << "nodes:" << nodes_.size() << ", ";
    result << "arcs:" << arcs_.size() << ")";
    return result.str();
  }

protected:
  FacetSPtr facet_b_;
  FacetSPtr facet_f_;

  std::list<ArcSPtr> arcs_;
  std::list<NodeSPtr> nodes_;
  StraightSkeletonWPtr skel_;

  typename std::list<SheetSPtr>::iterator list_it_;
  Plane3SPtr plane_;

  int id_;
};

} // namespace SDS
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_SDS_SHEET_H */
