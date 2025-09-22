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
 * @file   data/3d/skel/DblEdgeMergeEvent.h
 * @author Gernot Walzl
 * @date   2012-10-30
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_DBL_EDGE_MERGE_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_DBL_EDGE_MERGE_EVENT_H

#include <CGAL/Straight_skeleton_3/internal/debug.h>
#include <CGAL/Straight_skeleton_3/IO/String_factory.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/events/Abstract_event.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/internal/SDS/Straight_skeleton.h>

#include <array>
#include <memory>
#include <string>
#include <sstream>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
class DblEdgeMergeEvent
  : public AbstractEvent<Traits>
{
  using Base = AbstractEvent<Traits>;
  using DblEdgeMergeEventSPtr = std::shared_ptr<DblEdgeMergeEvent<Traits> >;

private:
  using FT = typename Traits::FT;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

private:
  using StraightSkeleton = SDS::StraightSkeleton<Traits>;
  using NodeSPtr = typename StraightSkeleton::NodeSPtr;

public:
  DblEdgeMergeEvent()
    : Base(Base::DBL_EDGE_MERGE_EVENT)
  { }

  virtual ~DblEdgeMergeEvent()
  {
    node_.reset();
  }

  static DblEdgeMergeEventSPtr create()
  {
    return std::make_shared<DblEdgeMergeEvent>();
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

  EdgeSPtr getEdge11() const
  {
    CGAL_SS3_DEBUG_WPTR(edge_11_);
    return edge_11_.lock();
  }

  void setEdge11(const EdgeSPtr& edge_11)
  {
    CGAL_SS3_DEBUG_SPTR(edge_11);
    this->edge_11_ = edge_11;
  }

  EdgeSPtr getEdge12() const
  {
    CGAL_SS3_DEBUG_WPTR(edge_12_);
    return edge_12_.lock();
  }

  void setEdge12(const EdgeSPtr& edge_12)
  {
    CGAL_SS3_DEBUG_SPTR(edge_12);
    this->edge_12_ = edge_12;
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

  EdgeSPtr getEdge21() const
  {
    CGAL_SS3_DEBUG_WPTR(edge_21_);
    return edge_21_.lock();
  }

  void setEdge21(const EdgeSPtr& edge_21)
  {
    CGAL_SS3_DEBUG_SPTR(edge_21);
    this->edge_21_ = edge_21;
  }

  EdgeSPtr getEdge22() const
  {
    CGAL_SS3_DEBUG_WPTR(edge_22_);
    return edge_22_.lock();
  }

  void setEdge22(const EdgeSPtr& edge_22)
  {
    CGAL_SS3_DEBUG_SPTR(edge_22);
    this->edge_22_ = edge_22;
  }

  void getVertices(VertexSPtr out[4]) const
  {
    EdgeSPtr edge_11 = getEdge11();
    EdgeSPtr edge_21 = getEdge21();
    EdgeSPtr edge_12 = getEdge12();
    EdgeSPtr edge_22 = getEdge22();
    FacetSPtr facet_1 = getFacet1();
    FacetSPtr facet_2 = getFacet2();

    for (unsigned int i = 0; i < 4; ++i) {
      out[i] = VertexSPtr();
    }
    out[0] = edge_11->dst(facet_1);
    out[1] = edge_21->dst(facet_2);
    out[2] = edge_12->src(facet_1);
    out[3] = edge_22->src(facet_2);
  }

  void getEdges(EdgeSPtr out[4]) const
  {
    EdgeSPtr edge_11 = getEdge11();
    EdgeSPtr edge_12 = getEdge12();
    FacetSPtr facet_1 = getFacet1();

    for (unsigned int i = 0; i < 4; ++i) {
      out[i] = EdgeSPtr();
    }
    out[0] = edge_11->next(facet_1);
    out[1] = out[0]->next(facet_1);
    FacetSPtr facet_other = edge_11->other(facet_1);
    out[2] = edge_12->next(facet_other);
    out[3] = out[2]->next(facet_other);
  }

  bool isValid() const
  {
    return (node_ &&
            !facet_1_.expired() && !edge_11_.expired() && !edge_12_.expired() &&
            !facet_2_.expired() && !edge_21_.expired() && !edge_22_.expired());
  }

  bool isObsolete() const
  {
    CGAL_assertion_msg(false, "NYI");
    return false;
  }

  std::string toString() const
  {
    VertexSPtr vertices[4];
    getVertices(vertices);

    EdgeSPtr edges[4];
    getEdges(edges);

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "DblEdgeMergeEvent\n";
    sstr << "\t(ID=" << Base::getID() << ")\n";
    sstr << "\t(offset=" << IO::StringFactory::fromDouble(CGAL::to_double(getTime())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")";
    for (unsigned int i = 0; i < 4; ++i) {
      sstr << "\n\t(" << vertices[i]->toString() << ")";
    }
    for (unsigned int i = 0; i < 4; ++i) {
      sstr << "\n\t(edge " << edges[i]->getID() << ")";
    }

    return sstr.str();
  }

  bool operator==(const DblEdgeMergeEvent& other) const
  {
    return (node_->getTime() == other.node_->getTime()) &&
            (*(node_->getPoint()) == *(other.node_->getPoint())) &&
            ((facet_1_.lock() == other.facet_1_.lock() &&
              edge_11_.lock() == other.edge_11_.lock() &&
              edge_12_.lock() == other.edge_12_.lock() &&
              facet_2_.lock() == other.facet_2_.lock() &&
              edge_21_.lock() == other.edge_21_.lock() &&
              edge_22_.lock() == other.edge_22_.lock()) ||
            (facet_1_.lock() == other.facet_2_.lock() &&
              edge_11_.lock() == other.edge_21_.lock() &&
              edge_12_.lock() == other.edge_22_.lock() &&
              facet_2_.lock() == other.facet_1_.lock() &&
              edge_21_.lock() == other.edge_11_.lock() &&
              edge_22_.lock() == other.edge_12_.lock()));
  }

protected:
  NodeSPtr node_;
  FacetWPtr facet_1_;
  EdgeWPtr edge_11_;
  EdgeWPtr edge_12_;
  FacetWPtr facet_2_;
  EdgeWPtr edge_21_;
  EdgeWPtr edge_22_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_DBL_EDGE_MERGE_EVENT_H */

