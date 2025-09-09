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
 * @file   data/3d/skel/PolyhedronSplitEvent.h
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_POLYHEDRON_SPLIT_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_POLYHEDRON_SPLIT_EVENT_H

#include <CGAL/Straight_skeleton_3/internal/debug.h>
#include <CGAL/Straight_skeleton_3/IO/String_factory.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/events/Abstract_event.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/internal/SDS/Straight_skeleton.h>

#include <memory>
#include <string>
#include <sstream>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
class PolyhedronSplitEvent
  : public AbstractEvent<Traits>
{
  using Base = AbstractEvent<Traits>;
  using PolyhedronSplitEventSPtr = std::shared_ptr<PolyhedronSplitEvent<Traits> >;

private:
  using FT = typename Traits::FT;

  using Polyhedron = HDS::Polyhedron<Traits>;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;

private:
  using StraightSkeleton = SDS::StraightSkeleton<Traits>;
  using NodeSPtr = typename StraightSkeleton::NodeSPtr;

private:
  using EdgeFacetNeighborhood = algorithm::EdgeFacetNeighborhood<Traits>;

public:
  PolyhedronSplitEvent()
    : Base(Base::POLYHEDRON_SPLIT_EVENT)
  { }

  virtual ~PolyhedronSplitEvent()
  {
    node_.reset();
    edge1_.reset();
    edge2_.reset();
  }

  static PolyhedronSplitEventSPtr create()
  {
    return std::make_shared<PolyhedronSplitEvent>();
  }

  NodeSPtr getNode() const
  {
    CGAL_SS3_DEBUG_SPTR(node_);
    return node_;
  }

  void setNode(NodeSPtr node)
  {
    CGAL_SS3_DEBUG_SPTR(node);
    this->node_ = node;
  }

  const FT& getOffset() const
  {
    CGAL_SS3_DEBUG_SPTR(node_);
    return node_->getOffset();
  }

  EdgeSPtr getEdge1() const
  {
    CGAL_SS3_DEBUG_WPTR(edge1_);
    return edge1_.lock();
  }

  void setEdge1(EdgeSPtr edge1)
  {
    CGAL_SS3_DEBUG_SPTR(edge1);
    this->edge1_ = edge1;
    this->neighborhood1_ = EdgeFacetNeighborhood(edge1);
  }

  EdgeSPtr getEdge2() const
  {
    CGAL_SS3_DEBUG_WPTR(edge2_);
    return edge2_.lock();
  }

  void setEdge2(EdgeSPtr edge2)
  {
    CGAL_SS3_DEBUG_SPTR(edge2);
    this->edge2_ = edge2;
    this->neighborhood2_ = EdgeFacetNeighborhood(edge2);
  }

  bool isValid() const
  {
    return node_ && !edge1_.expired() && !edge2_.expired();
  }

  bool isObsolete() const
  {
    if (EdgeSPtr edge_1 = getEdge1()) {
      if (!neighborhood1_.checkNeighborhoodConsistency(edge_1)) {
        return true;
      }
    }
    if (EdgeSPtr edge_2 = getEdge2()) {
      if (!neighborhood2_.checkNeighborhoodConsistency(edge_2)) {
        return true;
      }
    }
    return false;
  }

  std::string toString() const
  {
    EdgeSPtr edge1 = getEdge1();
    EdgeSPtr edge2 = getEdge2();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "PolyhedronSplitEvent\n";
    sstr << "\t(ID=" << Base::getID() << ")\n";
    sstr << "\t(offset=" << IO::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(edgeA=" << edge1->toString() << ")"
         << "; edgeB=" << edge2->toString() << ")";
    return sstr.str();
  }

  bool operator==(const PolyhedronSplitEvent& other) const
  {
    return (node_->getOffset() == other.node_->getOffset()) &&
            (*(node_->getPoint()) == *(other.node_->getPoint())) &&
            ((edge1_.lock() == other.edge1_.lock() &&
              edge2_.lock() == other.edge2_.lock()) ||
            (edge1_.lock() == other.edge2_.lock() &&
              edge2_.lock() == other.edge1_.lock()));
  }

protected:
  NodeSPtr node_;
  EdgeWPtr edge1_;
  EdgeWPtr edge2_;

  EdgeFacetNeighborhood neighborhood1_;
  EdgeFacetNeighborhood neighborhood2_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL
#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_POLYHEDRON_SPLIT_EVENT_H */

