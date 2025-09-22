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
 * @file   data/3d/skel/SplitMergeEvent.h
 * @author Gernot Walzl
 * @date   2013-08-09
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_SPLIT_MERGE_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_SPLIT_MERGE_EVENT_H

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
class SplitMergeEvent
  : public AbstractEvent<Traits>
{
  using Base = AbstractEvent<Traits>;
  using SplitMergeEventSPtr = std::shared_ptr<SplitMergeEvent<Traits> >;

private:
  using FT = typename Traits::FT;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using VertexWPtr = typename Polyhedron::VertexWPtr;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

private:
  using StraightSkeleton = SDS::StraightSkeleton<Traits>;
  using NodeSPtr = typename StraightSkeleton::NodeSPtr;

private:
  using VertexFacetNeighborhood = algorithm::VertexFacetNeighborhood<Traits>;

public:
  SplitMergeEvent()
    : Base(Base::SPLIT_MERGE_EVENT)
  { }

  virtual ~SplitMergeEvent()
  {
    node_.reset();
  }

  static SplitMergeEventSPtr create()
  {
    return std::make_shared<SplitMergeEvent>();
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

  const FT& getTime() const
  {
    CGAL_SS3_DEBUG_SPTR(node_);
    return node_->getTime();
  }

  VertexSPtr getVertex1() const
  {
    CGAL_SS3_DEBUG_WPTR(vertex_1_);
    return vertex_1_.lock();
  }

  void setVertex1(const VertexSPtr& vertex_1)
  {
    CGAL_SS3_DEBUG_SPTR(vertex_1);
    this->vertex_1_ = vertex_1;
    this->neighborhood_1_ = VertexFacetNeighborhood(vertex_1);
  }

  VertexSPtr getVertex2() const
  {
    CGAL_SS3_DEBUG_WPTR(vertex_2_);
    return vertex_2_.lock();
  }

  void setVertex2(const VertexSPtr& vertex_2)
  {
    CGAL_SS3_DEBUG_SPTR(vertex_2);
    this->vertex_2_ = vertex_2;
    this->neighborhood_2_ = VertexFacetNeighborhood(vertex_2);
  }

  FacetSPtr getFacet1() const
  {
    CGAL_SS3_DEBUG_WPTR(facet_1_);
    return facet_1_.lock();
  }

  void setFacet1(const FacetSPtr& facet_1)
  {
    CGAL_SS3_DEBUG_SPTR(facet_1);
    this->facet_1_ = facet_1;
  }

  FacetSPtr getFacet2() const
  {
    CGAL_SS3_DEBUG_WPTR(facet_2_);
    return facet_2_.lock();
  }

  void setFacet2(const FacetSPtr& facet_2)
  {
    CGAL_SS3_DEBUG_SPTR(facet_2);
    this->facet_2_ = facet_2;
  }

  bool isValid() const
  {
    return (node_ &&
            !vertex_1_.expired() && !vertex_2_.expired() &&
            !facet_1_.expired() && !facet_2_.expired());
  }

  bool isObsolete() const
  {
    if (VertexSPtr vertex_1 = getVertex1()) {
      if (!neighborhood_1_.checkNeighborhoodConsistency(vertex_1)) {
        return true;
      }
    }
    if (VertexSPtr vertex_2 = getVertex2()) {
      if (!neighborhood_2_.checkNeighborhoodConsistency(vertex_2)) {
        return true;
      }
    }
    return false;
  }

  std::string toString() const
  {
    VertexSPtr vertex_1 = getVertex1();
    VertexSPtr vertex_2 = getVertex2();
    FacetSPtr facet_1 = getFacet1();
    FacetSPtr facet_2 = getFacet2();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "SplitMergeEvent\n";
    sstr << "\t(ID=" << Base::getID() << ")\n";
    sstr << "\t(offset=" << IO::StringFactory::fromDouble(CGAL::to_double(getTime())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(vertex1=" << vertex_1->toString() << ")\n";
    sstr << "\t(vertex2=" << vertex_2->toString() << ")\n";
    sstr << "\t(facet1=" << facet_1->getID() << ")\n";
    sstr << "\t(facet2=" << facet_2->getID() << ")";
    return sstr.str();
  }

  bool operator==(const SplitMergeEvent& other) const
  {
    return (node_->getTime() == other.node_->getTime()) &&
            (*(node_->getPoint()) == *(other.node_->getPoint())) &&
            ((facet_1_.lock() == other.facet_1_.lock() &&
              facet_2_.lock() == other.facet_2_.lock()) ||
            (facet_1_.lock() == other.facet_2_.lock() &&
              facet_2_.lock() == other.facet_1_.lock())) &&
            ((vertex_1_.lock() == other.vertex_1_.lock() &&
              vertex_2_.lock() == other.vertex_2_.lock()) ||
            (vertex_1_.lock() == other.vertex_2_.lock() &&
              vertex_2_.lock() == other.vertex_1_.lock()));
  }

protected:
  NodeSPtr node_;
  VertexWPtr vertex_1_;
  VertexWPtr vertex_2_;
  FacetWPtr facet_1_;
  FacetWPtr facet_2_;

  VertexFacetNeighborhood neighborhood_1_;
  VertexFacetNeighborhood neighborhood_2_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_SPLIT_MERGE_EVENT_H */
