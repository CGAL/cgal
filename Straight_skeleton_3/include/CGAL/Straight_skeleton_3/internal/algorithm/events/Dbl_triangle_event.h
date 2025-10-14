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
 * file   data/3d/skel/DblTriangleEvent.h
 * author Gernot Walzl
 * date   2012-09-11
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
  using Point_3 = typename Traits::Point_3;
  using Point3SPtr = std::shared_ptr<Point_3>;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

public:
  DblTriangleEvent()
    : Base(Base::DBL_TRIANGLE_EVENT)
  { }

  virtual ~DblTriangleEvent()
  { }

  static DblTriangleEventSPtr create()
  {
    return std::make_shared<DblTriangleEvent>();
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
    return (!edge_.expired());
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
    sstr << "\t(time=" << IO::StringFactory::fromDouble(CGAL::to_double(Base::getTime())) << ")\n";
    if (point_) {
      sstr << "\t(point=<" + IO::StringFactory::fromDouble(CGAL::to_double(point_->x())) + " "
                           + IO::StringFactory::fromDouble(CGAL::to_double(point_->y())) + " "
                           + IO::StringFactory::fromDouble(CGAL::to_double(point_->z())) + ">)";
    }
    sstr << "\t(edge=" << edge->getID() << "\n\t\t[" << edge->getVertexSrc()->toString() << "\n\t\t "
                                                     << edge->getVertexDst()->toString() << "])";
    return sstr.str();
  }

  bool operator==(const DblTriangleEvent& other) const {
    return (Base::getTime() == other.getTime()) &&
            (!point_ || !other.point_ || *point_ == *(other.point_)) &&
            (edge_.lock() == other.edge_.lock());
  }

protected:
  Point3SPtr point_;
  EdgeWPtr edge_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_DBL_TRIANGLE_EVENT_H */

