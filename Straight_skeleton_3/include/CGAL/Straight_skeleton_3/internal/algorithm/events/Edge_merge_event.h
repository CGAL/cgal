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

#include <CGAL/license/Straight_skeleton_3.h>

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

template <typename GeomTraits>
class Edge_merge_event
  : public Abstract_event<GeomTraits>
{
  using Base = Abstract_event<GeomTraits>;
  using Edge_merge_event_sptr = std::shared_ptr<Edge_merge_event<GeomTraits> >;

private:
  using Point_3 = typename GeomTraits::Point_3;

private:
  using Polyhedron = HDS::Polyhedron<GeomTraits>;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

public:
  Edge_merge_event()
    : Base(Base::EDGE_MERGE_EVENT)
  { }

  virtual ~Edge_merge_event()
  { }

  static Edge_merge_event_sptr create()
  {
    return std::make_shared<Edge_merge_event>();
  }

  const Point_3& point() const
  {
    return point_;
  }

  void set_point(const Point_3& point)
  {
    this->point_ = point;
  }

  FacetSPtr get_facet() const
  {
    CGAL_SS3_DEBUG_WPTR(facet_);
    return facet_.lock();
  }

  void set_facet(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    this->facet_ = facet;
  }

  EdgeSPtr get_edge_1() const
  {
    CGAL_SS3_DEBUG_WPTR(edge1_);
    return edge1_.lock();
  }

  void set_edge_1(const EdgeSPtr& edge1)
  {
    CGAL_SS3_DEBUG_SPTR(edge1);
    this->edge1_ = edge1;
  }

  EdgeSPtr get_edge_2() const
  {
    CGAL_SS3_DEBUG_WPTR(edge2_);
    return edge2_.lock();
  }

  void set_edge_2(const EdgeSPtr& edge2)
  {
    this->edge2_ = edge2;
  }

  bool is_valid() const
  {
    return (!facet_.expired() && !edge1_.expired() && !edge2_.expired());
  }

  bool is_obsolete() const
  {
    CGAL_warning_msg(false, "NYI");
    return false;
  }

  std::string to_string() const
  {
    FacetSPtr facet = get_facet();
    EdgeSPtr edge1 = get_edge_1();
    EdgeSPtr edge2 = get_edge_2();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "Edge_merge_event\n";
    sstr << "\t(ID=" << Base::id() << ")\n";
    sstr << "\t(time=" << IO::String_factory::fromDouble(CGAL::to_double(Base::time())) << ")\n";
    sstr << "\t(point=<" + IO::String_factory::fromDouble(CGAL::to_double(point_.x())) + " "
                         + IO::String_factory::fromDouble(CGAL::to_double(point_.y())) + " "
                         + IO::String_factory::fromDouble(CGAL::to_double(point_.z())) + ">)";
    sstr << "\t(facet=" << facet->id() << ")\n";
    sstr << "\t(edgeA=" << edge1->id() << "\n\t\t[" << edge1->source()->to_string() << "\n\t\t "
                                                        << edge1->target()->to_string() << "])\n";
    sstr << "\t(edgeB=" << edge2->id() << "\n\t\t[" << edge2->source()->to_string() << "\n\t\t "
                                                        << edge2->target()->to_string() << "])";
    return sstr.str();
  }

  bool operator==(const Edge_merge_event& other) const
  {
    return (Base::time() == other.time()) &&
            (!point_ || !other.point_ || point_ == other.point_) &&
            (facet_.lock() == other.facet_.lock()) &&
            ((edge1_.lock() == other.edge1_.lock() &&
              edge2_.lock() == other.edge2_.lock()) ||
            (edge1_.lock() == other.edge2_.lock() &&
              edge2_.lock() == other.edge1_.lock()));
  }

protected:
  Point_3 point_;
  FacetWPtr facet_;
  EdgeWPtr edge1_;
  EdgeWPtr edge2_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_EDGE_MERGE_EVENT_H */

