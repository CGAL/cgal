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
 * file   data/3d/skel/TetrahedronEvent.h
 * author Gernot Walzl
 * date   2012-04-23
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_TETRAHEDRON_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_TETRAHEDRON_EVENT_H

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

template <typename Traits>
class Tetrahedron_event
  : public Abstract_event<Traits>
{
  using Base = Abstract_event<Traits>;
  using Tetrahedron_event_sptr = std::shared_ptr<Tetrahedron_event<Traits> >;

private:
  using Point_3 = typename Traits::Point_3;
  using Point3SPtr = std::shared_ptr<Point_3>;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

private:
  using Edge_facet_neighborhood = algorithm::Edge_facet_neighborhood<Traits>;

public:
  Tetrahedron_event()
    : Base(Base::TETRAHEDRON_EVENT)
  { }

  virtual ~Tetrahedron_event()
  { }

  static Tetrahedron_event_sptr create()
  {
    return std::make_shared<Tetrahedron_event>();
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

  EdgeSPtr get_edge_begin() const
  {
    CGAL_SS3_DEBUG_WPTR(edge_begin_);
    return edge_begin_.lock();
  }

  void set_edge_begin(const EdgeSPtr& edge_begin)
  {
    CGAL_SS3_DEBUG_SPTR(edge_begin);
    this->edge_begin_ = edge_begin;
    this->neighborhood_ = Edge_facet_neighborhood(edge_begin);
  }

  std::array<VertexSPtr, 4> get_vertices() const
  {
    EdgeSPtr edge_begin = get_edge_begin();
    return CGAL::make_array(edge_begin->get_vertex_src(),
                            edge_begin->get_vertex_dst(),
                            edge_begin->next(edge_begin->get_facet_L())->dst(edge_begin->get_facet_L()),
                            edge_begin->next(edge_begin->get_facet_R())->dst(edge_begin->get_facet_R()));
  }

  std::array<EdgeSPtr, 6> get_edges() const
  {
    EdgeSPtr edge_begin = get_edge_begin();
    EdgeSPtr e2 = edge_begin->next(edge_begin->get_facet_L());
    return CGAL::make_array(edge_begin,
                            edge_begin->prev(edge_begin->get_facet_L()),
                            e2,
                            edge_begin->prev(edge_begin->get_facet_R()),
                            edge_begin->next(edge_begin->get_facet_R()),
                            e2->prev(e2->other(edge_begin->get_facet_L())));
  }

  std::array<FacetSPtr, 4> get_facets() const
  {
    EdgeSPtr edge_begin = get_edge_begin();
    return CGAL::make_array(edge_begin->get_facet_L(),
                            edge_begin->get_facet_R(),
                            edge_begin->get_facet_L()->prev(edge_begin->get_vertex_dst()),
                            edge_begin->get_facet_R()->prev(edge_begin->get_vertex_src()));
  }

  bool is_valid() const
  {
    return (!edge_begin_.expired());
  }

  bool is_obsolete() const
  {
    if (EdgeSPtr edge = get_edge_begin()) {
      return ! neighborhood_.check_neighborhood_consistency(edge);
    }

    return false;
  }

  std::string to_string() const
  {
    std::array<VertexSPtr, 4> vertices = get_vertices();
    std::array<EdgeSPtr, 6> edges = get_edges();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "Tetrahedron_event\n";
    sstr << "\t(ID=" << Base::get_ID() << ")\n";
    sstr << "\t(time=" << IO::String_factory::fromDouble(CGAL::to_double(Base::time())) << ")\n";
    if (point_) {
      sstr << "\t(point=<" + IO::String_factory::fromDouble(CGAL::to_double(point_->x())) + " "
                           + IO::String_factory::fromDouble(CGAL::to_double(point_->y())) + " "
                           + IO::String_factory::fromDouble(CGAL::to_double(point_->z())) + ">)";
    }
    sstr << "\t(vertices";
    for (int i=0; i<4; ++i)
      sstr << " " << vertices[i]->get_ID();
    sstr << ")\n";
    sstr << "\t(edges";
    for (int i=0; i<6; ++i)
      sstr << " " << edges[i]->get_ID();
    sstr << ")";
    return sstr.str();
  }

  bool operator==(const Tetrahedron_event& other) const
  {
    return (Base::time() == other.time()) &&
            (!point_ || !other.point_ || *point_ == *(other.point_)) &&
            (edge_begin_.lock() == other.edge_begin_.lock());
  }

protected:
  Point3SPtr point_;
  EdgeWPtr edge_begin_;

  Edge_facet_neighborhood neighborhood_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_TETRAHEDRON_EVENT_H */

