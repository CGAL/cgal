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
#include <CGAL/Straight_skeleton_3/Straight_skeleton_3.h>

#include <CGAL/array.h>

#include <memory>
#include <string>
#include <sstream>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename GeomTraits>
class Dbl_triangle_event
  : public Abstract_event<GeomTraits>
{
  using Base = Abstract_event<GeomTraits>;
  using Dbl_triangle_event_sptr = std::shared_ptr<Dbl_triangle_event<GeomTraits> >;

private:
  using Point_3 = typename GeomTraits::Point_3;

private:
  using Polyhedron = HDS::Polyhedron<GeomTraits>;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

public:
  Dbl_triangle_event()
    : Base(Base::DBL_TRIANGLE_EVENT)
  { }

  virtual ~Dbl_triangle_event()
  { }

  static Dbl_triangle_event_sptr create()
  {
    return std::make_shared<Dbl_triangle_event>();
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
  }

  std::array<VertexSPtr, 4> get_vertices() const
  {
    EdgeSPtr edge = get_edge();
    FacetSPtr facet_l = edge->get_facet_L();
    FacetSPtr facet_r = edge->get_facet_R();
    return CGAL::make_array(edge->source(),
                            edge->target(),
                            edge->next(facet_l)->tgt(facet_l),
                            edge->next(facet_r)->tgt(facet_r));
  }

  std::array<EdgeSPtr, 5> get_edges() const
  {
    EdgeSPtr edge = get_edge();
    FacetSPtr facet_l = edge->get_facet_L();
    FacetSPtr facet_r = edge->get_facet_R();
    return CGAL::make_array(edge,
                            edge->next(facet_l),
                            edge->prev(facet_l),
                            edge->next(facet_r),
                            edge->prev(facet_r));
  }

  bool is_valid() const
  {
    return (!edge_.expired());
  }

  bool is_obsolete() const
  {
    CGAL_warning_msg(false, "NYI");
    return false;
  }

  std::string to_string() const
  {
    EdgeSPtr edge = get_edge();
    std::array<VertexSPtr, 4> vertices = get_vertices();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "Dbl_triangle_event\n";
    sstr << "\t(ID=" << Base::get_ID() << ")\n";
    sstr << "\t(time=" << IO::String_factory::fromDouble(CGAL::to_double(Base::time())) << ")\n";
    sstr << "\t(point=<" + IO::String_factory::fromDouble(CGAL::to_double(point_.x())) + " "
                         + IO::String_factory::fromDouble(CGAL::to_double(point_.y())) + " "
                         + IO::String_factory::fromDouble(CGAL::to_double(point_.z())) + ">)";
    sstr << "\t(edge=" << edge->get_ID() << "\n\t\t[" << edge->source()->to_string() << "\n\t\t "
                                                      << edge->target()->to_string() << "])";
    return sstr.str();
  }

  bool operator==(const Dbl_triangle_event& other) const {
    return (Base::time() == other.time()) &&
            (!point_ || !other.point_ || point_ == other.point_) &&
            (edge_.lock() == other.edge_.lock());
  }

protected:
  Point_3 point_;
  EdgeWPtr edge_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_DBL_TRIANGLE_EVENT_H */

