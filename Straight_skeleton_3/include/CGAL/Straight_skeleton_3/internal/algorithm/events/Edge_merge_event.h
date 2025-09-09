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
 * @file   data/3d/skel/EdgeMergeEvent.h
 * @author Gernot Walzl
 * @date   2012-09-14
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_EDGE_MERGE_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_EDGE_MERGE_EVENT_H

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
class EdgeMergeEvent
  : public AbstractEvent<Traits>
{
  using Base = AbstractEvent<Traits>;
  using EdgeMergeEventSPtr = std::shared_ptr<EdgeMergeEvent<Traits> >;

private:
  using FT = typename Traits::FT;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

private:
  using StraightSkeleton = SDS::StraightSkeleton<Traits>;
  using NodeSPtr = typename StraightSkeleton::NodeSPtr;

public:
  EdgeMergeEvent()
    : Base(Base::EDGE_MERGE_EVENT)
  { }

  virtual ~EdgeMergeEvent()
  {
    node_.reset();
    facet_.reset();
    edge1_.reset();
    edge2_.reset();
  }

  static EdgeMergeEventSPtr create()
  {
    return std::make_shared<EdgeMergeEvent>();
  }

  NodeSPtr getNode() const {
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

  FacetSPtr getFacet() const
  {
    CGAL_SS3_DEBUG_WPTR(facet_);
    return facet_.lock();
  }

  void setFacet(FacetSPtr facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    this->facet_ = facet;
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
  }

  EdgeSPtr getEdge2() const
  {
    CGAL_SS3_DEBUG_WPTR(edge2_);
    return edge2_.lock();
  }

  void setEdge2(EdgeSPtr edge2)
  {
    this->edge2_ = edge2;
  }

  bool isValid() const
  {
    return (node_ && !facet_.expired() && !edge1_.expired() && !edge2_.expired());
  }

  bool isObsolete() const
  {
    CGAL_warning_msg(false, "NYI");
    return false;
  }

  std::string toString() const
  {
    FacetSPtr facet = getFacet();
    EdgeSPtr edge1 = getEdge1();
    EdgeSPtr edge2 = getEdge2();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "EdgeMergeEvent\n";
    sstr << "\t(ID=" << Base::getID() << ")\n";
    sstr << "\t(offset=" << IO::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(facet=" << facet->getID() << ")\n";
    sstr << "\t(edgeA=" << edge1->getID() << "\n\t\t[" << edge1->getVertexSrc()->toString() << "\n\t\t "
                                                       << edge1->getVertexDst()->toString() << "])\n";
    sstr << "\t(edgeB=" << edge2->getID() << "\n\t\t[" << edge2->getVertexSrc()->toString() << "\n\t\t "
                                                       << edge2->getVertexDst()->toString() << "])";
    return sstr.str();
  }

  bool operator==(const EdgeMergeEvent& other) const
  {
    return (node_->getOffset() == other.node_->getOffset()) &&
            (*(node_->getPoint()) == *(other.node_->getPoint())) &&
            (facet_.lock() == other.facet_.lock()) &&
            ((edge1_.lock() == other.edge1_.lock() &&
              edge2_.lock() == other.edge2_.lock()) ||
            (edge1_.lock() == other.edge2_.lock() &&
              edge2_.lock() == other.edge1_.lock()));
  }

protected:
  NodeSPtr node_;
  FacetWPtr facet_;
  EdgeWPtr edge1_;
  EdgeWPtr edge2_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_EDGE_MERGE_EVENT_H */

