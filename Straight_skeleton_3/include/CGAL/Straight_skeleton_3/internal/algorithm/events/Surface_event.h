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
 * file   data/3d/skel/SurfaceEvent.h
 * author Gernot Walzl
 * date   2012-09-10
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_SURFACE_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_SURFACE_EVENT_H

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
class Surface_event
  : public Abstract_event<Traits>
{
  using Base = Abstract_event<Traits>;

private:
  using Point_3 = typename Traits::Point_3;
  using Point3SPtr = std::shared_ptr<Point_3>;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;

private:
  using Surface_event_sptr = std::shared_ptr<Surface_event<Traits> >;

private:
  using Edge_facet_neighborhood = algorithm::Edge_facet_neighborhood<Traits>;

public:
  Surface_event()
    : Base(Base::SURFACE_EVENT)
  { }

  virtual ~Surface_event()
  { }

  Point3SPtr point() const
  {
    CGAL_SS3_DEBUG_SPTR(point_);
    return point_;
  }

  void set_point(const Point3SPtr& point)
  {
    this->point_ = point;
  }

  static Surface_event_sptr create()
  {
    return std::make_shared<Surface_event>();
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
    this->neighborhood1_ = Edge_facet_neighborhood(edge1);
  }

  EdgeSPtr get_edge_2() const
  {
    CGAL_SS3_DEBUG_WPTR(edge2_);
    return edge2_.lock();
  }

  void set_edge_2(const EdgeSPtr& edge2)
  {
    CGAL_SS3_DEBUG_SPTR(edge2);
    this->edge2_ = edge2;
    this->neighborhood2_ = Edge_facet_neighborhood(edge2);
  }

  bool is_valid() const
  {
    return (!edge1_.expired() && !edge2_.expired());
  }

  bool is_obsolete() const
  {
    if (EdgeSPtr edge_1 = get_edge_1()) {
      if (!neighborhood1_.check_neighborhood_consistency(edge_1)) {
        return true;
      }
    }
    if (EdgeSPtr edge_2 = get_edge_2()) {
      if (!neighborhood2_.check_neighborhood_consistency(edge_2)) {
        return true;
      }
    }
    return false;
  }

  std::string to_string() const
  {
    EdgeSPtr edge1 = get_edge_1();
    EdgeSPtr edge2 = get_edge_2();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "Surface_event\n";
    sstr << "\t(ID=" << Base::get_ID() << ")\n";
    sstr << "\t(time=" << IO::String_factory::fromDouble(CGAL::to_double(Base::time())) << ")\n";
    if (point_) {
      sstr << "\t(point=<" + IO::String_factory::fromDouble(CGAL::to_double(point_->x())) + " "
                           + IO::String_factory::fromDouble(CGAL::to_double(point_->y())) + " "
                           + IO::String_factory::fromDouble(CGAL::to_double(point_->z())) + ">)";
    }
    sstr << "\t(edgeA=" << edge1->get_ID() << "\n\t\t[" << edge1->get_vertex_src()->to_string() << "\n\t\t "
                                                       << edge1->get_vertex_dst()->to_string() << "]\n"
         << "\t edgeB=" << edge2->get_ID() << "\n\t\t[" << edge2->get_vertex_src()->to_string() << "\n\t\t "
                                                       << edge2->get_vertex_dst()->to_string() << "])";

    return sstr.str();
  }

bool operator==(const Surface_event& other) const
{
    return (Base::time() == other.time()) &&
           (!point_ || !other.point_ || *point_ == *(other.point_)) &&
           ((edge1_.lock() == other.edge1_.lock() &&
             edge2_.lock() == other.edge2_.lock()) ||
            (edge1_.lock() == other.edge2_.lock() &&
             edge2_.lock() == other.edge1_.lock()));
}


protected:
  Point3SPtr point_;
  EdgeWPtr edge1_;
  EdgeWPtr edge2_;

  Edge_facet_neighborhood neighborhood1_;
  Edge_facet_neighborhood neighborhood2_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL


#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_SURFACE_EVENT_H */

