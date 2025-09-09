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
 * @file   data/3d/skel/TetrahedronEvent.h
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_TETRAHEDRON_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_TETRAHEDRON_EVENT_H

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
class TetrahedronEvent
  : public AbstractEvent<Traits>
{
  using Base = AbstractEvent<Traits>;
  using TetrahedronEventSPtr = std::shared_ptr<TetrahedronEvent<Traits> >;

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

private:
  using EdgeFacetNeighborhood = algorithm::EdgeFacetNeighborhood<Traits>;

public:
  TetrahedronEvent()
    : Base(Base::TETRAHEDRON_EVENT)
  { }

  virtual ~TetrahedronEvent() {
    node_.reset();
  }

  static TetrahedronEventSPtr create()
  {
    return std::make_shared<TetrahedronEvent>();
  }

  NodeSPtr getNode() const
  {
    CGAL_SS3_DEBUG_SPTR(node_);
    return node_;
  }

  void setNode(NodeSPtr node)
  {
    this->node_ = node;
  }

  const FT& getOffset() const
  {
    CGAL_SS3_DEBUG_SPTR(node_);
    return node_->getOffset();
  }

  EdgeSPtr getEdgeBegin() const
  {
    CGAL_SS3_DEBUG_WPTR(edge_begin_);
    return edge_begin_.lock();
  }

  void setEdgeBegin(EdgeSPtr edge_begin)
  {
    this->edge_begin_ = edge_begin;
    this->neighborhood_ = EdgeFacetNeighborhood(edge_begin);
  }

  void getVertices(VertexSPtr out[4]) const
  {
    EdgeSPtr edge_begin = getEdgeBegin();

    for (unsigned int i = 0; i < 4; ++i) {
      out[i] = VertexSPtr();
    }
    out[0] = edge_begin->getVertexSrc();
    out[1] = edge_begin->getVertexDst();
    EdgeSPtr edge_l = edge_begin->next(edge_begin->getFacetL());
    out[2] = edge_l->dst(edge_begin->getFacetL());
    EdgeSPtr edge_r = edge_begin->next(edge_begin->getFacetR());
    out[3] = edge_r->dst(edge_begin->getFacetR());
  }

  void getEdges(EdgeSPtr out[6]) const
  {
    EdgeSPtr edge_begin = getEdgeBegin();

    for (unsigned int i = 0; i < 6; ++i) {
      out[i] = EdgeSPtr();
    }
    out[0] = edge_begin;
    out[1] = edge_begin->prev(edge_begin->getFacetL());
    out[2] = edge_begin->next(edge_begin->getFacetL());
    out[3] = edge_begin->prev(edge_begin->getFacetR());
    out[4] = edge_begin->next(edge_begin->getFacetR());
    FacetSPtr other = out[2]->other(edge_begin->getFacetL());
    out[5] = out[2]->prev(other);
  }

  void getFacets(FacetSPtr out[4]) const
  {
    EdgeSPtr edge_begin = getEdgeBegin();

    for (unsigned int i = 0; i < 4; ++i) {
      out[i] = FacetSPtr();
    }
    out[0] = edge_begin->getFacetL();
    out[1] = edge_begin->getFacetR();
    out[2] = out[0]->prev(edge_begin->getVertexDst());
    out[3] = out[1]->prev(edge_begin->getVertexSrc());
  }

  bool isValid() const
  {
    return (node_ && !edge_begin_.expired());
  }

  bool isObsolete() const
  {
    if (EdgeSPtr edge = getEdgeBegin()) {
      return ! neighborhood_.checkNeighborhoodConsistency(edge);
    }

    return false;
  }

  std::string toString() const
  {
    VertexSPtr vertices[4];
    getVertices(vertices);

    EdgeSPtr edges[6];
    getEdges(edges);

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "TetrahedronEvent\n";
    sstr << "\t(ID=" << Base::getID() << ")\n";
    sstr << "\t(offset=" << IO::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(vertices";
    for (int i=0; i<4; ++i)
      sstr << " " << vertices[i]->getID();
    sstr << ")\n";
    sstr << "\t(edges";
    for (int i=0; i<6; ++i)
      sstr << " " << edges[i]->getID();
    sstr << ")";
    return sstr.str();
  }

  bool operator==(const TetrahedronEvent& other) const
  {
    return (node_->getOffset() == other.node_->getOffset()) &&
            (*(node_->getPoint()) == *(other.node_->getPoint())) &&
            (edge_begin_.lock() == other.edge_begin_.lock());
  }

protected:
  NodeSPtr node_;
  EdgeWPtr edge_begin_;

  EdgeFacetNeighborhood neighborhood_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_TETRAHEDRON_EVENT_H */

