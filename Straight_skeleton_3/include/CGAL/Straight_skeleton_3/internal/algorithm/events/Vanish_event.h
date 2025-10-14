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

template <typename Traits>
class VanishEvent
  : public AbstractEvent<Traits>
{
  using Base = AbstractEvent<Traits>;
  using VanishEventSPtr = std::shared_ptr<VanishEvent<Traits> >;

private:
  using Point_3 = typename Traits::Point_3;
  using Point3SPtr = std::shared_ptr<Point_3>;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;

private:
  using EdgeFacetNeighborhood = algorithm::EdgeFacetNeighborhood<Traits>;

public:
  VanishEvent()
    : Base(Base::VANISH_EVENT)
  { }

  virtual ~VanishEvent()
  { }

  static VanishEventSPtr create()
  {
    return std::make_shared<VanishEvent>();
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
    this->neighborhood_ = EdgeFacetNeighborhood(edge);
  }

  bool isValid() const
  {
    return (!edge_.expired());
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
    sstr << "\t(time=" << IO::StringFactory::fromDouble(CGAL::to_double(Base::getTime())) << ")\n";
    if (point_) {
      sstr << "\t(point=<" + IO::StringFactory::fromDouble(CGAL::to_double(point_->x())) + " "
                           + IO::StringFactory::fromDouble(CGAL::to_double(point_->y())) + " "
                           + IO::StringFactory::fromDouble(CGAL::to_double(point_->z())) + ">)";
    }
    sstr << "\t(edgeA=" << edge->getID() << "\n\t\t[" << edge->getVertexSrc()->toString() << "\n\t\t "
                                                      << edge->getVertexDst()->toString() << "])";
    return sstr.str();
  }

  bool operator==(const VanishEvent& other) const
  {
    return (Base::getTime() == other.getTime()) &&
            // && (edge_.lock() == other.edge_.lock()) // because of multiple reps...
            (*(getPoint()) == *(other.getPoint()));
  }

protected:
  Point3SPtr point_;
  EdgeWPtr edge_;

  EdgeFacetNeighborhood neighborhood_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_VANISH_EVENT_H */

