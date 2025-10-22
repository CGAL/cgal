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
 * file   data/3d/skel/EdgeEvent.h
 * author Gernot Walzl
 * date   2012-04-23
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_EDGE_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_EDGE_EVENT_H

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
class Edge_event
  : public Abstract_event<Traits>
{
  using Base = Abstract_event<Traits>;
  using Edge_event_sptr = std::shared_ptr<Edge_event<Traits> >;

private:
  using Point_3 = typename Traits::Point_3;
  using Point3SPtr = std::shared_ptr<Point_3>;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;

private:
  using Edge_facet_neighborhood = algorithm::Edge_facet_neighborhood<Traits>;

public:
  Edge_event()
    : Base(Base::EDGE_EVENT)
  { }

  virtual ~Edge_event()
  { }

  static Edge_event_sptr create()
  {
    return std::make_shared<Edge_event>();
  }

  Point3SPtr point() const
  {
    CGAL_SS3_DEBUG_SPTR(point_);
    return point_;
  }

  void set_point(const Point3SPtr& point)
  {
    this->point_ = point;
  }

  EdgeSPtr get_edge() const
  {
    CGAL_SS3_DEBUG_WPTR(edge_);
    return edge_.lock();
  }

  void set_edge(const EdgeSPtr& e)
  {
    CGAL_SS3_DEBUG_SPTR(e);
    this->edge_ = e;
    this->neighborhood_ = Edge_facet_neighborhood(e);
  }

  bool is_valid() const
  {
    return (!edge_.expired());
  }

  bool is_obsolete() const
  {
    if (EdgeSPtr edge = get_edge()) {
      return ! neighborhood_.check_neighborhood_consistency(edge);
    }
    return false;
  }

  std::string to_string() const
  {
    EdgeSPtr edge = get_edge();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "Edge_event\n";
    sstr << "\t(ID=" << Base::get_ID() << ")\n";
    sstr << "\t(time=" << IO::String_factory::fromDouble(CGAL::to_double(Base::time())) << ")\n";
    if (point_) {
      sstr << "\t(point=<" + IO::String_factory::fromDouble(CGAL::to_double(point_->x())) + " "
                           + IO::String_factory::fromDouble(CGAL::to_double(point_->y())) + " "
                           + IO::String_factory::fromDouble(CGAL::to_double(point_->z())) + ">)";
    }
    sstr << "\t(edgeA=" << edge->get_ID() << "\n\t\t[" << edge->get_vertex_src()->to_string() << "\n\t\t "
                                                       << edge->get_vertex_dst()->to_string() << "])";
    return sstr.str();
  }

  bool operator==(const Edge_event& other) const
  {
    return (Base::time() == other.time()) &&
            (!point_ || !other.point_ || *point_ == *(other.point_)) &&
            (edge_.lock() == other.edge_.lock());
  }

protected:
  Point3SPtr point_;
  EdgeWPtr edge_;

  Edge_facet_neighborhood neighborhood_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_EDGE_EVENT_H */

