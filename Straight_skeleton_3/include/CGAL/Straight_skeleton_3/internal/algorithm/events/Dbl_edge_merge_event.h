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
 * file   data/3d/skel/DblEdgeMergeEvent.h
 * author Gernot Walzl
 * date   2012-10-30
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_DBL_EDGE_MERGE_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_DBL_EDGE_MERGE_EVENT_H

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
class Dbl_edge_merge_event
  : public Abstract_event<Traits>
{
  using Base = Abstract_event<Traits>;
  using Dbl_edge_merge_event_sptr = std::shared_ptr<Dbl_edge_merge_event<Traits> >;

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

public:
  Dbl_edge_merge_event()
    : Base(Base::DBL_EDGE_MERGE_EVENT)
  { }

  virtual ~Dbl_edge_merge_event()
  { }

  static Dbl_edge_merge_event_sptr create()
  {
    return std::make_shared<Dbl_edge_merge_event>();
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

  EdgeSPtr get_edge_11() const
  {
    CGAL_SS3_DEBUG_WPTR(edge_11_);
    return edge_11_.lock();
  }

  void set_edge_11(const EdgeSPtr& edge_11)
  {
    CGAL_SS3_DEBUG_SPTR(edge_11);
    this->edge_11_ = edge_11;
  }

  EdgeSPtr get_edge_12() const
  {
    CGAL_SS3_DEBUG_WPTR(edge_12_);
    return edge_12_.lock();
  }

  void set_edge_12(const EdgeSPtr& edge_12)
  {
    CGAL_SS3_DEBUG_SPTR(edge_12);
    this->edge_12_ = edge_12;
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

  EdgeSPtr get_edge_21() const
  {
    CGAL_SS3_DEBUG_WPTR(edge_21_);
    return edge_21_.lock();
  }

  void set_edge_21(const EdgeSPtr& edge_21)
  {
    CGAL_SS3_DEBUG_SPTR(edge_21);
    this->edge_21_ = edge_21;
  }

  EdgeSPtr get_edge_22() const
  {
    CGAL_SS3_DEBUG_WPTR(edge_22_);
    return edge_22_.lock();
  }

  void set_edge_22(const EdgeSPtr& edge_22)
  {
    CGAL_SS3_DEBUG_SPTR(edge_22);
    this->edge_22_ = edge_22;
  }

  std::array<VertexSPtr, 4> get_vertices() const
  {
    EdgeSPtr edge_11 = get_edge_11();
    EdgeSPtr edge_21 = get_edge_21();
    EdgeSPtr edge_12 = get_edge_12();
    EdgeSPtr edge_22 = get_edge_22();
    FacetSPtr facet_1 = get_facet_1();
    FacetSPtr facet_2 = get_facet_2();

    return CGAL::make_array(edge_11->dst(facet_1),
                            edge_21->dst(facet_2),
                            edge_12->src(facet_1),
                            edge_22->src(facet_2));
  }

  std::array<EdgeSPtr, 4> get_edges() const
  {
    EdgeSPtr edge_11 = get_edge_11();
    EdgeSPtr edge_12 = get_edge_12();
    FacetSPtr facet_1 = get_facet_1();
    FacetSPtr facet_other = edge_11->other(facet_1);

    return CGAL::make_array(edge_11->next(facet_1),
                            edge_11->next(facet_1)->next(facet_1),
                            edge_12->next(facet_other),
                            edge_12->next(facet_other)->next(facet_other));
  }

  bool is_valid() const
  {
    return (!facet_1_.expired() && !edge_11_.expired() && !edge_12_.expired() &&
            !facet_2_.expired() && !edge_21_.expired() && !edge_22_.expired());
  }

  bool is_obsolete() const
  {
    CGAL_assertion_msg(false, "NYI");
    return false;
  }

  std::string to_string() const
  {
    std::array<VertexSPtr, 4> vertices = get_vertices();
    std::array<EdgeSPtr, 4> edges = get_edges();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "Dbl_edge_merge_event\n";
    sstr << "\t(ID=" << Base::get_ID() << ")\n";
    sstr << "\t(time=" << IO::String_factory::fromDouble(CGAL::to_double(Base::time())) << ")\n";
    if (point_) {
      sstr << "\t(point=<" + IO::String_factory::fromDouble(CGAL::to_double(point_->x())) + " "
                           + IO::String_factory::fromDouble(CGAL::to_double(point_->y())) + " "
                           + IO::String_factory::fromDouble(CGAL::to_double(point_->z())) + ">)";
    }
    for (unsigned int i = 0; i < 4; ++i) {
      sstr << "\n\t(" << vertices[i]->to_string() << ")";
    }
    for (unsigned int i = 0; i < 4; ++i) {
      sstr << "\n\t(edge " << edges[i]->get_ID() << ")";
    }

    return sstr.str();
  }

  bool operator==(const Dbl_edge_merge_event& other) const
  {
    return (Base::time() == other.time()) &&
            (!point_ || !other.point_ || *point_ == *(other.point_)) &&
            ((facet_1_.lock() == other.facet_1_.lock() &&
              edge_11_.lock() == other.edge_11_.lock() &&
              edge_12_.lock() == other.edge_12_.lock() &&
              facet_2_.lock() == other.facet_2_.lock() &&
              edge_21_.lock() == other.edge_21_.lock() &&
              edge_22_.lock() == other.edge_22_.lock()) ||
            (facet_1_.lock() == other.facet_2_.lock() &&
              edge_11_.lock() == other.edge_21_.lock() &&
              edge_12_.lock() == other.edge_22_.lock() &&
              facet_2_.lock() == other.facet_1_.lock() &&
              edge_21_.lock() == other.edge_11_.lock() &&
              edge_22_.lock() == other.edge_12_.lock()));
  }

protected:
  Point3SPtr point_;
  FacetWPtr facet_1_;
  EdgeWPtr edge_11_;
  EdgeWPtr edge_12_;
  FacetWPtr facet_2_;
  EdgeWPtr edge_21_;
  EdgeWPtr edge_22_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_DBL_EDGE_MERGE_EVENT_H */

