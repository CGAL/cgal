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
 * file   data/3d/skel/FlipVertexEvent.h
 * author Gernot Walzl
 * date   2012-11-22
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_FLIP_VERTEX_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_FLIP_VERTEX_EVENT_H

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

template <typename Traits>
class Flip_vertex_event
  : public Abstract_event<Traits>
{
  using Base = Abstract_event<Traits>;
  using Flip_vertex_event_sptr = std::shared_ptr<Flip_vertex_event<Traits> >;

private:
  using Point_3 = typename Traits::Point_3;
  using Point3SPtr = std::shared_ptr<Point_3>;

  using Polyhedron = HDS::Polyhedron<Traits>;
  using VertexWPtr = typename Polyhedron::VertexWPtr;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

private:
  using Vertex_facet_neighborhood = algorithm::Vertex_facet_neighborhood<Traits>;

public:
  Flip_vertex_event()
    : Base(Base::FLIP_VERTEX_EVENT)
  { }

  virtual ~Flip_vertex_event()
  { }

  static Flip_vertex_event_sptr create()
  {
    return std::make_shared<Flip_vertex_event>();
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

  VertexSPtr get_vertex_1() const
  {
    CGAL_SS3_DEBUG_WPTR(vertex_1_);
    return vertex_1_.lock();
  }

  void set_vertex_1(const VertexSPtr& vertex_1)
  {
    CGAL_SS3_DEBUG_SPTR(vertex_1);
    this->vertex_1_ = vertex_1;
    this->neighborhood_1_ = Vertex_facet_neighborhood(vertex_1);
  }

  VertexSPtr get_vertex_2() const
  {
    CGAL_SS3_DEBUG_WPTR(vertex_2_);
    return vertex_2_.lock();
  }

  void set_vertex_2(const VertexSPtr& vertex_2)
  {
    CGAL_SS3_DEBUG_SPTR(vertex_2);
    this->vertex_2_ = vertex_2;
    this->neighborhood_2_ = Vertex_facet_neighborhood(vertex_2);
  }

  FacetSPtr get_facet_1() const
  {
    CGAL_SS3_DEBUG_WPTR(facet_1_);
    return facet_1_.lock();
  }

  void set_facet_1(const FacetSPtr& facet_1)
  {
    CGAL_SS3_DEBUG_SPTR(facet_1);
    this->facet_1_ = facet_1;
  }

  FacetSPtr get_facet_2() const
  {
    CGAL_SS3_DEBUG_WPTR(facet_2_);
    return facet_2_.lock();
  }

  void set_facet_2(const FacetSPtr& facet_2)
  {
    CGAL_SS3_DEBUG_SPTR(facet_2);
    this->facet_2_ = facet_2;
  }

  bool is_valid() const
  {
    return (!vertex_1_.expired() && !vertex_2_.expired() &&
            !facet_1_.expired() && !facet_2_.expired());
  }

  bool is_obsolete() const
  {
    if (VertexSPtr vertex_1 = get_vertex_1()) {
      if (!neighborhood_1_.check_neighborhood_consistency(vertex_1)) {
        return true;
      }
    }
    if (VertexSPtr vertex_2 = get_vertex_2()) {
      if (!neighborhood_2_.check_neighborhood_consistency(vertex_2)) {
        return true;
      }
    }
    return false;
  }

  std::string to_string() const
  {
    VertexSPtr vertex_1 = get_vertex_1();
    VertexSPtr vertex_2 = get_vertex_2();
    FacetSPtr facet_1 = get_facet_1();
    FacetSPtr facet_2 = get_facet_2();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "Flip_vertex_event\n";
    sstr << "\t(ID=" << Base::get_ID() << ")\n";
    sstr << "\t(time=" << IO::String_factory::fromDouble(CGAL::to_double(Base::time())) << ")\n";
    sstr << "\t(vertex A=" << vertex_1->to_string() << "; vertex B=" << vertex_2->to_string() << ")\n";
    sstr << "\t(facet A=" << facet_1->get_ID() << "; facet B=" << facet_2->get_ID() << ")";
    return sstr.str();
  }

  bool operator==(const Flip_vertex_event& other) const
  {
    return (Base::time() == other.time()) &&
            (!point_ || !other.point_ || *point_ == *(other.point_)) &&
            ((facet_1_.lock() == other.facet_1_.lock() &&
              facet_2_.lock() == other.facet_2_.lock()) ||
            (facet_1_.lock() == other.facet_2_.lock() &&
              facet_2_.lock() == other.facet_1_.lock())) &&
            ((vertex_1_.lock() == other.vertex_1_.lock() &&
              vertex_2_.lock() == other.vertex_2_.lock()) ||
            (vertex_1_.lock() == other.vertex_2_.lock() &&
              vertex_2_.lock() == other.vertex_1_.lock()));
  }

protected:
  Point3SPtr point_;
  VertexWPtr vertex_1_;
  VertexWPtr vertex_2_;
  FacetWPtr facet_1_;
  FacetWPtr facet_2_;

  Vertex_facet_neighborhood neighborhood_1_;
  Vertex_facet_neighborhood neighborhood_2_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_FLIP_VERTEX_EVENT_H */

