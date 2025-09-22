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
 * @file   data/3d/skel/DblTriangleEvent.h
 * @author Gernot Walzl
 * @date   2012-09-11
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_DBL_TRIANGLE_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_DBL_TRIANGLE_EVENT_H

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
class DblTriangleEvent
  : public AbstractEvent<Traits>
{
  using Base = AbstractEvent<Traits>;
  using DblTriangleEventSPtr = std::shared_ptr<DblTriangleEvent<Traits> >;

private:
  using FT = typename Traits::FT;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

private:
  using StraightSkeleton = SDS::StraightSkeleton<Traits>;
  using NodeSPtr = typename StraightSkeleton::NodeSPtr;

public:
  DblTriangleEvent()
    : Base(Base::DBL_TRIANGLE_EVENT)
  { }

  virtual ~DblTriangleEvent()
  {
    node_.reset();
  }

  static DblTriangleEventSPtr create()
  {
    return std::make_shared<DblTriangleEvent>();
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

  EdgeSPtr getEdge() const
  {
    CGAL_SS3_DEBUG_WPTR(edge_);
    return edge_.lock();
  }

  void setEdge(const EdgeSPtr& edge)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    this->edge_ = edge;
  }

  void getVertices(VertexSPtr out[4]) const
  {
    EdgeSPtr edge = getEdge();

    for (unsigned int i = 0; i < 4; ++i) {
      out[i] = VertexSPtr();
    }
    FacetSPtr facet_l = edge->getFacetL();
    FacetSPtr facet_r = edge->getFacetR();
    out[0] = edge->getVertexSrc();
    out[1] = edge->getVertexDst();
    out[2] = edge->next(facet_l)->dst(facet_l);
    out[3] = edge->next(facet_r)->dst(facet_r);
  }

  void getEdges(EdgeSPtr out[5]) const
  {
    EdgeSPtr edge = getEdge();

    for (unsigned int i = 0; i < 5; ++i) {
      out[i] = EdgeSPtr();
    }
    FacetSPtr facet_l = edge->getFacetL();
    FacetSPtr facet_r = edge->getFacetR();
    out[0] = edge;
    out[1] = edge->next(facet_l);
    out[2] = edge->prev(facet_l);
    out[3] = edge->next(facet_r);
    out[4] = edge->prev(facet_r);
  }

  bool isValid() const
  {
    return (node_ && !edge_.expired());
  }

  bool isObsolete() const
  {
    CGAL_warning_msg(false, "NYI");
    return false;
  }

  std::string toString() const
  {
    EdgeSPtr edge = getEdge();

    VertexSPtr vertices[4];
    getVertices(vertices);

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "DblTriangleEvent\n";
    sstr << "\t(ID=" << Base::getID() << ")\n";
    sstr << "\t(offset=" << IO::StringFactory::fromDouble(CGAL::to_double(getTime())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(edge=" << edge->getID() << "\n\t\t[" << edge->getVertexSrc()->toString() << "\n\t\t "
                                                     << edge->getVertexDst()->toString() << "])";
    return sstr.str();
  }

  bool operator==(const DblTriangleEvent& other) const {
    return (node_->getTime() == other.node_->getTime()) &&
            (*(node_->getPoint()) == *(other.node_->getPoint())) &&
            (edge_.lock() == other.edge_.lock());
  }

protected:
  NodeSPtr node_;
  EdgeWPtr edge_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_DBL_TRIANGLE_EVENT_H */

