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
 * file   data/3d/skel/PierceEvent.h
 * author Gernot Walzl
 * date   2012-04-23
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_PIERCE_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_PIERCE_EVENT_H

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
class Pierce_event
  : public Abstract_event<GeomTraits>
{
  using Base = Abstract_event<GeomTraits>;
  using Pierce_event_sptr = std::shared_ptr<Pierce_event<GeomTraits> >;

private:
  using Point_3 = typename GeomTraits::Point_3;

private:
  using Polyhedron = HDS::Polyhedron<GeomTraits>;
  using VertexWPtr = typename Polyhedron::VertexWPtr;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

private:
  using Vertex_facet_neighborhood = algorithm::Vertex_facet_neighborhood<GeomTraits>;

public:
  Pierce_event()
    : Base(Base::PIERCE_EVENT)
  { }

  virtual ~Pierce_event()
  { }

  static Pierce_event_sptr create()
  {
    return std::make_shared<Pierce_event>();
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

  VertexSPtr get_vertex() const
  {
    CGAL_SS3_DEBUG_WPTR(vertex_);
    return vertex_.lock();
  }

  void set_vertex(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    this->vertex_ = vertex;
    this->neighborhood_ = Vertex_facet_neighborhood(vertex);
  }

  bool is_valid() const
  {
    return (!facet_.expired() && !vertex_.expired());
  }

  bool is_obsolete() const
  {
    if (VertexSPtr vertex = get_vertex()) {
      return ! neighborhood_.check_neighborhood_consistency(vertex);
    }
    return false;
  }

  std::string to_string() const
  {
    FacetSPtr facet = get_facet();
    VertexSPtr vertex = get_vertex();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "Pierce_event\n";
    sstr << "\t(ID=" << Base::id() << ")\n";
    sstr << "\t(time=" << IO::String_factory::fromDouble(CGAL::to_double(Base::time())) << ")\n";
    sstr << "\t(point=<" + IO::String_factory::fromDouble(CGAL::to_double(point_.x())) + " "
                         + IO::String_factory::fromDouble(CGAL::to_double(point_.y())) + " "
                         + IO::String_factory::fromDouble(CGAL::to_double(point_.z())) + ">)";
    sstr << "\t(vertex=" << vertex->to_string() << ")\n";
    sstr << "\t(facet=" << facet->id() << ")";
    return sstr.str();
  }

  bool operator==(const Pierce_event& other) const
  {
    return (Base::time() == other.time()) &&
            (!point_ || !other.point_ || point_ == other.point_) &&
            (facet_.lock() == other.facet_.lock()) &&
            (vertex_.lock() == other.vertex_.lock());
  }

protected:
  Point_3 point_;
  FacetWPtr facet_;
  VertexWPtr vertex_;

  Vertex_facet_neighborhood neighborhood_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_PIERCE_EVENT_H */
