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

#include <CGAL/license/Straight_skeleton_3.h>

#include <CGAL/Straight_skeleton_3/internal/debug.h>
#include <CGAL/Straight_skeleton_3/IO/String_factory.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/events/Abstract_event.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/Straight_skeleton_3.h>

#include <memory>
#include <string>
#include <sstream>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename GeomTraits>
class Vanish_event
  : public Abstract_event<GeomTraits>
{
  using Base = Abstract_event<GeomTraits>;
  using Vanish_event_sptr = std::shared_ptr<Vanish_event<GeomTraits> >;

private:
  using Point_3 = typename GeomTraits::Point_3;

private:
  using Polyhedron = HDS::Polyhedron<GeomTraits>;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;

private:
  using Edge_facet_neighborhood = algorithm::Edge_facet_neighborhood<GeomTraits>;

public:
  Vanish_event()
    : Base(Base::VANISH_EVENT)
  { }

  virtual ~Vanish_event()
  { }

  static Vanish_event_sptr create()
  {
    return std::make_shared<Vanish_event>();
  }

  const Point_3& point() const
  {
    return point_;
  }

  void set_point(const Point_3& point)
  {
    this->point_ = point;
  }

  EdgeSPtr get_edge() const
  {
    CGAL_SS3_DEBUG_WPTR(edge_);
    return edge_.lock();
  }

  void set_edge(const EdgeSPtr& edge)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    this->edge_ = edge;
    this->neighborhood_ = Edge_facet_neighborhood(edge);
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
    sstr << "Vanish_event\n";
    sstr << "\t(ID=" << Base::id() << ")\n";
    sstr << "\t(time=" << IO::String_factory::fromDouble(CGAL::to_double(Base::time())) << ")\n";
    sstr << "\t(point=<" + IO::String_factory::fromDouble(CGAL::to_double(point_.x())) + " "
                         + IO::String_factory::fromDouble(CGAL::to_double(point_.y())) + " "
                         + IO::String_factory::fromDouble(CGAL::to_double(point_.z())) + ">)";
    sstr << "\t(edgeA=" << edge->id() << "\n\t\t[" << edge->source()->to_string() << "\n\t\t "
                                                       << edge->target()->to_string() << "])";
    return sstr.str();
  }

  bool operator==(const Vanish_event& other) const
  {
    return (Base::time() == other.time()) &&
            // && (edge_.lock() == other.edge_.lock()) // because of multiple reps...
            (point() == other.point());
  }

protected:
  Point_3 point_;
  EdgeWPtr edge_;

  Edge_facet_neighborhood neighborhood_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_VANISH_EVENT_H */

