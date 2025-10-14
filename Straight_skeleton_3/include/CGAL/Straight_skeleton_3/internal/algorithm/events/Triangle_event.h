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
 * file   data/3d/skel/TriangleEvent.h
 * author Gernot Walzl
 * date   2012-04-23
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_TRIANGLE_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_TRIANGLE_EVENT_H

#include <CGAL/Straight_skeleton_3/internal/debug.h>
#include <CGAL/Straight_skeleton_3/IO/String_factory.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/events/Abstract_event.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>

#include <memory>
#include <string>
#include <sstream>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
class TriangleEvent
  : public AbstractEvent<Traits>
{
  using Base = AbstractEvent<Traits>;
  using TriangleEventSPtr = std::shared_ptr<TriangleEvent<Traits> >;

private:
  using Point_3 = typename Traits::Point_3;
  using Point3SPtr = std::shared_ptr<Point_3>;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

private:
  using EdgeFacetNeighborhood = algorithm::EdgeFacetNeighborhood<Traits>;

public:
  TriangleEvent()
    : Base(Base::TRIANGLE_EVENT)
  { }

  virtual ~TriangleEvent()
  { }

  static TriangleEventSPtr create()
  {
    return std::make_shared<TriangleEvent>();
  }

  Point3SPtr getPoint() const
  {
    CGAL_SS3_DEBUG_SPTR(point_);
    return point_;
  }

  void setPoint(const Point3SPtr& point)
  {
    this->point_ = point;
  }

  FacetSPtr getFacet() const
  {
    CGAL_SS3_DEBUG_WPTR(facet_);
    return facet_.lock();
  }

  void setFacet(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    this->facet_ = facet;
  }

  EdgeSPtr getEdgeBegin() const
  {
    CGAL_SS3_DEBUG_WPTR(edge_begin_);
    return edge_begin_.lock();
  }

  void setEdgeBegin(const EdgeSPtr& edge_begin) {
    CGAL_SS3_DEBUG_SPTR(edge_begin);
    this->edge_begin_ = edge_begin;
    this->neighborhood_ = EdgeFacetNeighborhood(edge_begin);
  }

  void getVertices(VertexSPtr out[3]) const
  {
    FacetSPtr facet = getFacet();
    EdgeSPtr edge_begin = getEdgeBegin();

    for (unsigned int i = 0; i < 3; ++i) {
      out[i] = VertexSPtr();
    }
    out[0] = edge_begin->src(facet);
    out[1] = edge_begin->dst(facet);
    out[2] = edge_begin->next(facet)->dst(facet);
  }

  void getEdges(EdgeSPtr out[3]) const
  {
    FacetSPtr facet = getFacet();
    EdgeSPtr edge_begin = getEdgeBegin();

    for (unsigned int i = 0; i < 3; ++i) {
      out[i] = EdgeSPtr();
    }
    out[0] = edge_begin;
    out[1] = edge_begin->next(facet);
    out[2] = edge_begin->prev(facet);
  }

  bool isValid() const
  {
    return (!facet_.expired() && !edge_begin_.expired());
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
    FacetSPtr facet = getFacet();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "TriangleEvent\n";
    sstr << "\t(ID=" << Base::getID() << ")\n";
    sstr << "\t(time=" << IO::StringFactory::fromDouble(CGAL::to_double(Base::getTime())) << ")\n";
    if (point_) {
      sstr << "\t(point=<" + IO::StringFactory::fromDouble(CGAL::to_double(point_->x())) + " "
                           + IO::StringFactory::fromDouble(CGAL::to_double(point_->y())) + " "
                           + IO::StringFactory::fromDouble(CGAL::to_double(point_->z())) + ">)";
    }
    sstr << "\t(facet=" << facet->getID() << ")";
    return sstr.str();
  }

  bool operator==(const TriangleEvent& other) const
  {
    return (Base::getTime() == other.getTime()) &&
            (!point_ || !other.point_ || *point_ == *(other.point_)) &&
            (facet_.lock() == other.facet_.lock()) &&
            (edge_begin_.lock() == other.edge_begin_.lock());
  }

protected:
  Point3SPtr point_;
  FacetWPtr facet_; // @todo shouldn't be needed, edge_begin_->getFacetL is enough
  EdgeWPtr edge_begin_;

  EdgeFacetNeighborhood neighborhood_; // this covers the four faces involved in the triangle event
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_TRIANGLE_EVENT_H */

