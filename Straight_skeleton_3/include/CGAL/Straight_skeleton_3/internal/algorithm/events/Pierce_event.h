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
 * file   data/3d/skel/PierceEvent.h
 * author Gernot Walzl
 * date   2012-04-23
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_PIERCE_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_PIERCE_EVENT_H

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
class PierceEvent
  : public AbstractEvent<Traits>
{
  using Base = AbstractEvent<Traits>;
  using PierceEventSPtr = std::shared_ptr<PierceEvent<Traits> >;

private:
  using Point_3 = typename Traits::Point_3;
  using Point3SPtr = std::shared_ptr<Point_3>;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using VertexWPtr = typename Polyhedron::VertexWPtr;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

private:
  using VertexFacetNeighborhood = algorithm::VertexFacetNeighborhood<Traits>;

public:
  PierceEvent()
    : Base(Base::PIERCE_EVENT)
  { }

  virtual ~PierceEvent()
  { }

  static PierceEventSPtr create()
  {
    return std::make_shared<PierceEvent>();
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

  VertexSPtr getVertex() const
  {
    CGAL_SS3_DEBUG_WPTR(vertex_);
    return vertex_.lock();
  }

  void setVertex(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    this->vertex_ = vertex;
    this->neighborhood_ = VertexFacetNeighborhood(vertex);
  }

  bool isValid() const
  {
    return (!facet_.expired() && !vertex_.expired());
  }

  bool isObsolete() const
  {
    if (VertexSPtr vertex = getVertex()) {
      return ! neighborhood_.checkNeighborhoodConsistency(vertex);
    }
    return false;
  }

  std::string toString() const
  {
    FacetSPtr facet = getFacet();
    VertexSPtr vertex = getVertex();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "PierceEvent\n";
    sstr << "\t(ID=" << Base::getID() << ")\n";
    sstr << "\t(time=" << IO::StringFactory::fromDouble(CGAL::to_double(Base::getTime())) << ")\n";
    if (point_) {
      sstr << "\t(point=<" + IO::StringFactory::fromDouble(CGAL::to_double(point_->x())) + " "
                           + IO::StringFactory::fromDouble(CGAL::to_double(point_->y())) + " "
                           + IO::StringFactory::fromDouble(CGAL::to_double(point_->z())) + ">)";
    }
    sstr << "\t(vertex=" << vertex->toString() << ")\n";
    sstr << "\t(facet=" << facet->getID() << ")";
    return sstr.str();
  }

  bool operator==(const PierceEvent& other) const
  {
    return (Base::getTime() == other.getTime()) &&
            (!point_ || !other.point_ || *point_ == *(other.point_)) &&
            (facet_.lock() == other.facet_.lock()) &&
            (vertex_.lock() == other.vertex_.lock());
  }

protected:
  Point3SPtr point_;
  FacetWPtr facet_;
  VertexWPtr vertex_;

  VertexFacetNeighborhood neighborhood_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_PIERCE_EVENT_H */
