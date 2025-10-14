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
 * file   data/3d/skel/EdgeMergeEvent.h
 * author Gernot Walzl
 * date   2012-09-14
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_EDGE_MERGE_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_EDGE_MERGE_EVENT_H

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
class EdgeMergeEvent
  : public AbstractEvent<Traits>
{
  using Base = AbstractEvent<Traits>;
  using EdgeMergeEventSPtr = std::shared_ptr<EdgeMergeEvent<Traits> >;

private:
  using Point_3 = typename Traits::Point_3;
  using Point3SPtr = std::shared_ptr<Point_3>;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

public:
  EdgeMergeEvent()
    : Base(Base::EDGE_MERGE_EVENT)
  { }

  virtual ~EdgeMergeEvent()
  { }

  static EdgeMergeEventSPtr create()
  {
    return std::make_shared<EdgeMergeEvent>();
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

  EdgeSPtr getEdge1() const
  {
    CGAL_SS3_DEBUG_WPTR(edge1_);
    return edge1_.lock();
  }

  void setEdge1(const EdgeSPtr& edge1)
  {
    CGAL_SS3_DEBUG_SPTR(edge1);
    this->edge1_ = edge1;
  }

  EdgeSPtr getEdge2() const
  {
    CGAL_SS3_DEBUG_WPTR(edge2_);
    return edge2_.lock();
  }

  void setEdge2(const EdgeSPtr& edge2)
  {
    this->edge2_ = edge2;
  }

  bool isValid() const
  {
    return (!facet_.expired() && !edge1_.expired() && !edge2_.expired());
  }

  bool isObsolete() const
  {
    CGAL_warning_msg(false, "NYI");
    return false;
  }

  std::string toString() const
  {
    FacetSPtr facet = getFacet();
    EdgeSPtr edge1 = getEdge1();
    EdgeSPtr edge2 = getEdge2();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "EdgeMergeEvent\n";
    sstr << "\t(ID=" << Base::getID() << ")\n";
    sstr << "\t(time=" << IO::StringFactory::fromDouble(CGAL::to_double(Base::getTime())) << ")\n";
    if (point_) {
      sstr << "\t(point=<" + IO::StringFactory::fromDouble(CGAL::to_double(point_->x())) + " "
                           + IO::StringFactory::fromDouble(CGAL::to_double(point_->y())) + " "
                           + IO::StringFactory::fromDouble(CGAL::to_double(point_->z())) + ">)";
    }
    sstr << "\t(facet=" << facet->getID() << ")\n";
    sstr << "\t(edgeA=" << edge1->getID() << "\n\t\t[" << edge1->getVertexSrc()->toString() << "\n\t\t "
                                                       << edge1->getVertexDst()->toString() << "])\n";
    sstr << "\t(edgeB=" << edge2->getID() << "\n\t\t[" << edge2->getVertexSrc()->toString() << "\n\t\t "
                                                       << edge2->getVertexDst()->toString() << "])";
    return sstr.str();
  }

  bool operator==(const EdgeMergeEvent& other) const
  {
    return (Base::getTime() == other.getTime()) &&
            (!point_ || !other.point_ || *point_ == *(other.point_)) &&
            (facet_.lock() == other.facet_.lock()) &&
            ((edge1_.lock() == other.edge1_.lock() &&
              edge2_.lock() == other.edge2_.lock()) ||
            (edge1_.lock() == other.edge2_.lock() &&
              edge2_.lock() == other.edge1_.lock()));
  }

protected:
  Point3SPtr point_;
  FacetWPtr facet_;
  EdgeWPtr edge1_;
  EdgeWPtr edge2_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_EDGE_MERGE_EVENT_H */

