// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_PERIODIC_3_MESH_3_TRIANGULATION_HELPERS_H
#define CGAL_PERIODIC_3_MESH_3_TRIANGULATION_HELPERS_H

#include <CGAL/license/Periodic_3_mesh_3.h>

#include <CGAL/Mesh_3/Triangulation_helpers.h>

#include <boost/mpl/identity.hpp>

namespace CGAL {
namespace Mesh_3 {

template <typename Tr, typename IsPeriodicTriangulation>
class Triangulation_helpers;

// specialize for periodic triangulations
template <typename Tr>
class Triangulation_helpers<Tr, CGAL::Tag_true /*Periodic_tag*/>
  : public Triangulation_helpers<Tr, CGAL::Tag_false> // inherit the non-periodic functions
{
  typedef typename Tr::Geom_traits              GT;
  typedef typename Tr::FT                       FT;
  typedef typename GT::Vector_3                 Vector_3;
  typedef typename GT::Triangle_3               Triangle_3;

  // If `Tr` is not a triangulation that has defined Bare_point,
  // use Point_3 as defined in the traits class.
  typedef typename boost::mpl::eval_if_c<
    CGAL::internal::Has_nested_type_Bare_point<Tr>::value,
    typename CGAL::internal::Bare_point_type<Tr>,
    boost::mpl::identity<typename GT::Point_3>
  >::type                                       Bare_point;

  typedef typename GT::Weighted_point_3         Weighted_point_3;

  typedef typename Tr::Offset                   Offset;

  typedef typename Tr::Vertex                   Vertex;
  typedef typename Tr::Vertex_handle            Vertex_handle;
  typedef typename Tr::Facet                    Facet;
  typedef typename Tr::Cell                     Cell;
  typedef typename Tr::Cell_handle              Cell_handle;

public:
  Bare_point canonicalize_point(const Tr& tr, const Bare_point& p) const
  {
    return P3T3::internal::construct_canonical_point(p, tr.geom_traits());
  }

  // @fixme it might be dangerous to call robust_canonicalize() without also changing
  // <p, offset> = construct_periodic_point(p) (lack of consistency in the result)
  Weighted_point_3 canonicalize_point(const Tr& tr, const Weighted_point_3& p) const
  {
    return P3T3::internal::construct_canonical_point(p, tr.geom_traits());
  }

  // Warning: This function finds which offset 'Oq' should be applied to 'q' such
  // that the distance between 'p' and '(q, Oq)' is minimal.
  //
  // \pre 'p' lives in the canonical instance.
  Bare_point get_closest_point(const Tr& tr,
                               const Bare_point& p,
                               const Bare_point& q) const
  {
    CGAL_precondition(p.x() < tr.domain().xmax());
    CGAL_precondition(p.y() < tr.domain().ymax());
    CGAL_precondition(p.z() < tr.domain().zmax());
    CGAL_precondition(p.x() >= tr.domain().xmin());
    CGAL_precondition(p.y() >= tr.domain().ymin());
    CGAL_precondition(p.z() >= tr.domain().zmin());

    typename GT::Compare_squared_distance_3 compare_sd =
      tr.geom_traits().compare_squared_distance_3_object();
    typename GT::Construct_point_3 cp =
      tr.geom_traits().construct_point_3_object();

    std::pair<Bare_point, Offset> pp_q = P3T3::internal::construct_periodic_point(q, tr.geom_traits());

    Offset min_off;
    Offset null_offset(0,0,0);

    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        for(int k = 0; k < 3; ++k)
        {
          const Offset o(i-1, j-1, k-1);

          if((i == 0 && j == 0 && k == 0) ||
             compare_sd(p, q, p, q,
                        null_offset, pp_q.second + o,
                        null_offset, pp_q.second + min_off) == SMALLER)
          {
            min_off = o;
          }
        }
      }
    }

    return cp(q, pp_q.second + min_off);
  }

  Weighted_point_3 get_closest_point(const Tr& tr,
                                     const Weighted_point_3& wp,
                                     const Weighted_point_3& wq) const
  {
    typename GT::Compute_weight_3 cw = tr.geom_traits().compute_weight_3_object();
    typename GT::Construct_point_3 cp = tr.geom_traits().construct_point_3_object();
    typename GT::Construct_weighted_point_3 cwp = tr.geom_traits().construct_weighted_point_3_object();

    return cwp(get_closest_point(tr, cp(wp), cp(wq)), cw(wq));
  }

  // returns the triangle corresponding to f, with a geometric shift
  // so that it is incident to ref_v's canonical position
  Triangle_3 get_incident_triangle(const Tr& tr,
                                   const Facet& f,
                                   const Vertex_handle ref_v) const
  {
    typename GT::Construct_point_3 cp = tr.geom_traits().construct_point_3_object();
    typename GT::Construct_translated_point_3 tp = tr.geom_traits().construct_translated_point_3_object();
    typename GT::Construct_vector_3 cv = tr.geom_traits().construct_vector_3_object();
    typename GT::Construct_triangle_3 ct = tr.geom_traits().construct_triangle_3_object();

    CGAL_precondition(f.first != Cell_handle() && f.first->has_vertex(ref_v));
    const int ref_v_pos = f.first->index(ref_v);
    const Bare_point& ref_p = cp(tr.point(ref_v));
    const Bare_point ref_p_in_f = cp(tr.point(f.first, ref_v_pos));
    Vector_3 move_to_canonical = cv(ref_p_in_f, ref_p);

    const int s = f.second;
    const Bare_point mp0 = tp(cp(tr.point(f.first, (s+1)%4)), move_to_canonical);
    const Bare_point mp1 = tp(cp(tr.point(f.first, (s+2)%4)), move_to_canonical);
    const Bare_point mp2 = tp(cp(tr.point(f.first, (s+3)%4)), move_to_canonical);
    const Triangle_3 t = ct(mp0, mp1, mp2);

    return t;
  }

  void set_point(const Tr& tr,
                 const Vertex_handle v,
                 const Vector_3& move,
                 const Weighted_point_3& new_position) const
  {
    // calling robust canonical here means we don't necessarily have
    // canonical(v + move) = new_position... @fixme
    return tr.set_point(v, move, canonicalize_point(tr, new_position));
  }

  FT compute_power_distance_to_power_sphere(const Tr& tr,
                                            const Cell_handle c,
                                            const int i) const
  {
    typename GT::Compute_power_distance_to_power_sphere_3 cr =
      tr.geom_traits().compute_power_distance_to_power_sphere_3_object();

    Offset o_nb = tr.neighbor_offset(c, i, c->neighbor(i));
    Offset o_vt = tr.get_offset(c->neighbor(i), c->neighbor(i)->index(c));

    const Weighted_point_3& wp0 = tr.point(c->vertex(0)); // need the canonical point
    const Weighted_point_3& wp1 = tr.point(c->vertex(1));
    const Weighted_point_3& wp2 = tr.point(c->vertex(2));
    const Weighted_point_3& wp3 = tr.point(c->vertex(3));
    const Weighted_point_3& wq = tr.point(c->neighbor(i)->vertex(c->neighbor(i)->index(c)));
    const Offset& op0 = tr.get_offset(c, 0);
    const Offset& op1 = tr.get_offset(c, 1);
    const Offset& op2 = tr.get_offset(c, 2);
    const Offset& op3 = tr.get_offset(c, 3);
    const Offset oq = o_vt - o_nb;

    return cr(wp0, wp1, wp2, wp3, wq, op0, op1, op2, op3, oq);
  }

  // The functions below are used in Mesh_3 and need a specific implementation
  // for the periodic case because we need to try with different offsets to get the result
  FT compute_power_distance_to_power_sphere(const Tr& tr,
                                            const Cell_handle c,
                                            const Vertex_handle v) const
  {
    // @fixme need to introduce Compare_power_distances_to_power_sphere_3(4 points, query)
    typename GT::Compute_power_distance_to_power_sphere_3 cr =
      tr.geom_traits().compute_power_distance_to_power_sphere_3_object();

    FT min_power_dist = std::numeric_limits<FT>::infinity();

    const Weighted_point_3& wp0 = tr.point(c->vertex(0)); // need the canonical point
    const Weighted_point_3& wp1 = tr.point(c->vertex(1));
    const Weighted_point_3& wp2 = tr.point(c->vertex(2));
    const Weighted_point_3& wp3 = tr.point(c->vertex(3));
    const Weighted_point_3& wq = tr.point(v);
    const Offset& op0 = tr.get_offset(c, 0);
    const Offset& op1 = tr.get_offset(c, 1);
    const Offset& op2 = tr.get_offset(c, 2);
    const Offset& op3 = tr.get_offset(c, 3);

    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        for(int k = 0; k < 3; ++k) {
          const Offset oq(i-1, j-1, k-1);

          FT power_dist = cr(wp0, wp1, wp2, wp3, wq, op0, op1, op2, op3, oq);

          if(power_dist < min_power_dist)
            min_power_dist = power_dist;
        }
      }
    }

    return min_power_dist;
  }

  // Warning: This is a periodic version that computes the smallest possible distance
  // between 'p' and 'q', for all possible combinations of offsets
  FT min_squared_distance(const Tr& tr,
                          const Bare_point& p,
                          const Bare_point& q) const
  {
    typename GT::Compare_squared_distance_3 compare_sd =
      tr.geom_traits().compare_squared_distance_3_object();
    typename GT::Compute_squared_distance_3 compute_sd =
      tr.geom_traits().compute_squared_distance_3_object();

    std::pair<Bare_point, Offset> pp_p = P3T3::internal::construct_periodic_point(p, tr.geom_traits());
    std::pair<Bare_point, Offset> pp_q = P3T3::internal::construct_periodic_point(q, tr.geom_traits());

    Offset min_off;

    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        for(int k = 0; k < 3; ++k)
        {
          const Offset o(i-1, j-1, k-1);

          if((i == 0 && j == 0 && k == 0) ||
              compare_sd(q, p, q, p,
                         pp_q.second, pp_p.second + o,
                         pp_q.second, pp_p.second + min_off) == SMALLER)
          {
            min_off = o;
          }
        }
      }
    }

    return compute_sd(q, p, pp_q.second, pp_p.second + min_off);
  }
};

} // namespace Mesh_3
} // namespace CGAL

#endif // CGAL_PERIODIC_3_MESH_3_TRIANGULATION_HELPERS_H
