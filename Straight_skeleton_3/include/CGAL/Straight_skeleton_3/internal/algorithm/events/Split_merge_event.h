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
 * file   data/3d/skel/SplitMergeEvent.h
 * author Gernot Walzl
 * date   2013-08-09
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_SPLIT_MERGE_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_SPLIT_MERGE_EVENT_H

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
class Split_merge_event
  : public Abstract_event<GeomTraits>
{
  using Base = Abstract_event<GeomTraits>;
  using Split_merge_event_sptr = std::shared_ptr<Split_merge_event<GeomTraits> >;

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
  Split_merge_event()
    : Base(Base::SPLIT_MERGE_EVENT)
  { }

  virtual ~Split_merge_event()
  { }

  static Split_merge_event_sptr create()
  {
    return std::make_shared<Split_merge_event>();
  }

  const Point_3& point() const
  {
    return point_;
  }

  void set_point(const Point_3& point)
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
    sstr << "Split_merge_event\n";
    sstr << "\t(ID=" << Base::id() << ")\n";
    sstr << "\t(time=" << IO::String_factory::fromDouble(CGAL::to_double(Base::time())) << ")\n";
    sstr << "\t(point=<" + IO::String_factory::fromDouble(CGAL::to_double(point_.x())) + " "
                         + IO::String_factory::fromDouble(CGAL::to_double(point_.y())) + " "
                         + IO::String_factory::fromDouble(CGAL::to_double(point_.z())) + ">)";
    sstr << "\t(vertex1=" << vertex_1->to_string() << ")\n";
    sstr << "\t(vertex2=" << vertex_2->to_string() << ")\n";
    sstr << "\t(facet1=" << facet_1->id() << ")\n";
    sstr << "\t(facet2=" << facet_2->id() << ")";
    return sstr.str();
  }

  bool operator==(const Split_merge_event& other) const
  {
    return (Base::time() == other.time()) &&
            (!point_ || !other.point_ || point_ == other.point_) &&
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
  Point_3 point_;
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

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_SPLIT_MERGE_EVENT_H */
