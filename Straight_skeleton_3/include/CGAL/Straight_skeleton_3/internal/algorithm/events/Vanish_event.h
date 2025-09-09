// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_VANISH_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_VANISH_EVENT_H

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

// @todo also add all the other boiler code as for other events (DAO, etc.)
template <typename Traits>
class VanishEvent
  : public AbstractEvent<Traits>
{
  using Base = AbstractEvent<Traits>;
  using VanishEventSPtr = std::shared_ptr<VanishEvent<Traits> >;

private:
  using FT = typename Traits::FT;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;

private:
  using StraightSkeleton = SDS::StraightSkeleton<Traits>;
  using NodeSPtr = typename StraightSkeleton::NodeSPtr;

private:
  using EdgeFacetNeighborhood = algorithm::EdgeFacetNeighborhood<Traits>;

public:
  VanishEvent()
    : Base(Base::VANISH_EVENT)
  { }

  virtual ~VanishEvent()
  {
    node_.reset();
    edge_.reset(); // @fixme is there still a point since edge_ is a weak pointer?
  }

  static VanishEventSPtr create()
  {
    return std::make_shared<VanishEvent>();
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

  EdgeSPtr getEdge() const
  {
    CGAL_SS3_DEBUG_WPTR(edge_);
    return edge_.lock();
  }

  void setEdge(EdgeSPtr edge)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    this->edge_ = edge;
    this->neighborhood_ = EdgeFacetNeighborhood(edge);
  }

  bool isValid() const
  {
    return node_ && !edge_.expired();
  }

  bool isObsolete() const
  {
    if (EdgeSPtr edge = getEdge()) {
      return ! neighborhood_.checkNeighborhoodConsistency(edge);
    }
    return false;
  }

  std::string toString() const
  {
    EdgeSPtr edge = getEdge();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "VanishEvent\n";
    sstr << "\t(ID=" << Base::getID() << ")\n";
    sstr << "\t(offset=" << IO::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(edgeA=" << edge->getID() << "\n\t\t[" << edge->getVertexSrc()->toString() << "\n\t\t "
                                                      << edge->getVertexDst()->toString() << "])";
    return sstr.str();
  }

  bool operator==(const VanishEvent& other) const
  {
    return (node_->getOffset() == other.node_->getOffset()) &&
            // && (edge_.lock() == other.edge_.lock()) // because of multiple reps...
            (*(node_->getPoint()) == *(other.node_->getPoint()));
  }

protected:
  NodeSPtr node_;
  EdgeWPtr edge_;

  EdgeFacetNeighborhood neighborhood_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_VANISH_EVENT_H */

