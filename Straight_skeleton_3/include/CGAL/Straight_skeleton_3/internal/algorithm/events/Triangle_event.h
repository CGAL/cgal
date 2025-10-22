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
 * file   data/3d/skel/TriangleEvent.h
 * author Gernot Walzl
 * date   2012-04-23
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_TRIANGLE_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_TRIANGLE_EVENT_H

#include <CGAL/Straight_skeleton_3/internal/debug.h>
#include <CGAL/Straight_skeleton_3/IO/String_factory.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/events/Abstract_event.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>

#include <CGAL/array.h>

#include <memory>
#include <string>
#include <sstream>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
class Triangle_event
  : public Abstract_event<Traits>
{
  using Base = Abstract_event<Traits>;
  using Triangle_event_sptr = std::shared_ptr<Triangle_event<Traits> >;

private:
  using Point_3 = typename Traits::Point_3;
  using Point3SPtr = std::shared_ptr<Point_3>;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

private:
  using Edge_facet_neighborhood = algorithm::Edge_facet_neighborhood<Traits>;

public:
  Triangle_event()
    : Base(Base::TRIANGLE_EVENT)
  { }

  virtual ~Triangle_event()
  { }

  static Triangle_event_sptr create()
  {
    return std::make_shared<Triangle_event>();
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

  EdgeSPtr get_edge_begin() const
  {
    CGAL_SS3_DEBUG_WPTR(edge_begin_);
    return edge_begin_.lock();
  }

  void set_edge_begin(const EdgeSPtr& edge_begin) {
    CGAL_SS3_DEBUG_SPTR(edge_begin);
    this->edge_begin_ = edge_begin;
    this->neighborhood_ = Edge_facet_neighborhood(edge_begin);
  }

  std::array<VertexSPtr, 3> get_vertices() const
  {
    EdgeSPtr edge_begin = get_edge_begin();
    FacetSPtr facet = get_facet();
    return CGAL::make_array(edge_begin->src(facet),
                            edge_begin->dst(facet),
                            edge_begin->next(facet)->dst(facet));
  }

  std::array<EdgeSPtr, 3> get_edges() const
  {
    FacetSPtr facet = get_facet();
    EdgeSPtr edge_begin = get_edge_begin();
    return CGAL::make_array(edge_begin,
                            edge_begin->next(facet),
                            edge_begin->prev(facet));
  }

  bool is_valid() const
  {
    return (!facet_.expired() && !edge_begin_.expired());
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
    FacetSPtr facet = get_facet();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "Triangle_event\n";
    sstr << "\t(ID=" << Base::get_ID() << ")\n";
    sstr << "\t(time=" << IO::String_factory::fromDouble(CGAL::to_double(Base::time())) << ")\n";
    if (point_) {
      sstr << "\t(point=<" + IO::String_factory::fromDouble(CGAL::to_double(point_->x())) + " "
                           + IO::String_factory::fromDouble(CGAL::to_double(point_->y())) + " "
                           + IO::String_factory::fromDouble(CGAL::to_double(point_->z())) + ">)";
    }
    sstr << "\t(facet=" << facet->get_ID() << ")";
    return sstr.str();
  }

  bool operator==(const Triangle_event& other) const
  {
    return (Base::time() == other.time()) &&
            (!point_ || !other.point_ || *point_ == *(other.point_)) &&
            (facet_.lock() == other.facet_.lock()) &&
            (edge_begin_.lock() == other.edge_begin_.lock());
  }

protected:
  Point3SPtr point_;
  FacetWPtr facet_; // @todo shouldn't be needed, edge_begin_->get_facet_L is enough
  EdgeWPtr edge_begin_;

  Edge_facet_neighborhood neighborhood_; // this covers the four faces involved in the triangle event
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_TRIANGLE_EVENT_H */

