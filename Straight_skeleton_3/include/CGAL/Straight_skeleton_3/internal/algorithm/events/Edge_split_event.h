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
 * file   data/3d/skel/EdgeSplitEvent.h
 * author Gernot Walzl
 * date   2012-04-23
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_EDGE_SPLIT_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_EDGE_SPLIT_EVENT_H

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
class EdgeSplitEvent
  : public AbstractEvent<Traits>
{
  using Base = AbstractEvent<Traits>;
  using EdgeSplitEventSPtr = std::shared_ptr<EdgeSplitEvent<Traits> >;

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
  EdgeSplitEvent()
    : Base(Base::EDGE_SPLIT_EVENT)
  { }

  virtual ~EdgeSplitEvent()
  { }

  static EdgeSplitEventSPtr create()
  {
    return std::make_shared<EdgeSplitEvent>();
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

  EdgeSPtr getEdge1() const
  {
    CGAL_SS3_DEBUG_WPTR(edge1_);
    return edge1_.lock();
  }

  void setEdge1(const EdgeSPtr& edge1)
  {
    CGAL_SS3_DEBUG_SPTR(edge1);
    this->edge1_ = edge1;
    this->neighborhood1_ = EdgeFacetNeighborhood(edge1);
  }

  EdgeSPtr getEdge2() const
  {
    CGAL_SS3_DEBUG_WPTR(edge2_);
    return edge2_.lock();
  }

  void setEdge2(const EdgeSPtr& edge2)
  {
    CGAL_SS3_DEBUG_SPTR(edge2);
    this->edge2_ = edge2;
    this->neighborhood2_ = EdgeFacetNeighborhood(edge2);
  }

  bool isValid() const
  {
    return (!edge1_.expired() && !edge2_.expired());
  }

  bool isObsolete() const
  {
    if (EdgeSPtr edge_1 = getEdge1()) {
      if (!neighborhood1_.checkNeighborhoodConsistency(edge_1)) {
        return true;
      }
    }
    if (EdgeSPtr edge_2 = getEdge2()) {
      if (!neighborhood2_.checkNeighborhoodConsistency(edge_2)) {
        return true;
      }
    }
    return false;
  }

  std::string toString() const
  {
    EdgeSPtr edge1 = getEdge1();
    EdgeSPtr edge2 = getEdge2();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "EdgeSplitEvent\n";
    sstr << "\t(ID=" << Base::getID() << ")\n";
    sstr << "\t(time=" << IO::StringFactory::fromDouble(CGAL::to_double(Base::getTime())) << ")\n";
    if (point_) {
      sstr << "\t(point=<" + IO::StringFactory::fromDouble(CGAL::to_double(point_->x())) + " "
                           + IO::StringFactory::fromDouble(CGAL::to_double(point_->y())) + " "
                           + IO::StringFactory::fromDouble(CGAL::to_double(point_->z())) + ">)";
    }
    sstr << "\t(edgeA=" << edge1->getID() << "[" << edge1->getVertexSrc()->getID() << "-"
                                                 << edge1->getVertexDst()->getID() << "]"
         << "; edgeB=" << edge2->getID() << "[" << edge2->getVertexSrc()->getID() << "-"
                                                << edge2->getVertexDst()->getID() << "])";
    return sstr.str();
  }

  bool operator==(const EdgeSplitEvent& other) const
  {
    return (Base::getTime() == other.getTime()) &&
            (!point_ || !other.point_ || *point_ == *(other.point_)) &&
            ((edge1_.lock() == other.edge1_.lock() &&
              edge2_.lock() == other.edge2_.lock()) ||
            (edge1_.lock() == other.edge2_.lock() &&
              edge2_.lock() == other.edge1_.lock()));
  }

protected:
  Point3SPtr point_;
  EdgeWPtr edge1_;
  EdgeWPtr edge2_;

  EdgeFacetNeighborhood neighborhood1_;
  EdgeFacetNeighborhood neighborhood2_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_EDGE_SPLIT_EVENT_H */

