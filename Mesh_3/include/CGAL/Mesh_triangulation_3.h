// Copyright (c) 2006-2009 INRIA Sophia-Antipolis (France).
// Copyright (c) 2011      GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau, Stephane Tayeb


#ifndef CGAL_MESH_TRIANGULATION_3_H
#define CGAL_MESH_TRIANGULATION_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Mesh_3/config.h>

#include <CGAL/Kernel_traits.h>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Robust_weighted_circumcenter_filtered_traits_3.h>

#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Compact_mesh_cell_base_3.h>
#include <CGAL/Mesh_3/io_signature.h>

namespace CGAL {

namespace details {

template<typename K>
struct Mesh_geom_traits_generator
{
private:
  typedef Robust_weighted_circumcenter_filtered_traits_3<K>   Geom_traits;

public:
  typedef Geom_traits                                         type;
  typedef type                                                Type;
};  // end struct Mesh_geom_traits_generator

} // end namespace details

template<class Gt_, class Tds_>
class Mesh_3_regular_triangulation_3_wrapper
  : public Regular_triangulation_3<Gt_, Tds_>
{
public:
  typedef Regular_triangulation_3<Gt_, Tds_>                  Base;

  typedef typename Base::Geom_traits                          Geom_traits;

  typedef typename Geom_traits::FT                            FT;
  typedef typename Base::Bare_point                           Bare_point;
  typedef typename Base::Weighted_point                       Weighted_point;
  typedef typename Base::Triangle                             Triangle;

  typedef typename Base::Vertex_handle                        Vertex_handle;
  typedef typename Base::Cell_handle                          Cell_handle;

  typedef typename Geom_traits::Vector_3                      Vector;

  using Base::geom_traits;
  using Base::point;

  static std::string io_signature() { return Get_io_signature<Base>()(); }

  // The undocumented, straightforward functions below are required for Mesh_3
  // because of Periodic_3_mesh_3: they are functions for which both triangulations
  // have fundamentally different implementations (usually, Periodic_3_mesh_3
  // does not know the offset of a point and must brute-force check for all
  // possibilities). To allow Periodic_3_mesh_3 to use Mesh_3's files,
  // each mesh triangulation implements its own version.

  const Bare_point& get_closest_point(const Bare_point& /*p*/, const Bare_point& q) const
  {
    return q;
  }

  const Triangle& get_closest_triangle(const Bare_point& /*p*/, const Triangle& t) const
  {
    return t;
  }

  void set_point(const Vertex_handle v,
                 const Vector& /*move*/,
                 const Weighted_point& new_position)
  {
    v->set_point(new_position);
  }

  FT compute_power_distance_to_power_sphere(const Cell_handle c,
                                            const Vertex_handle v) const
  {
    typedef typename Geom_traits::Compute_power_distance_to_power_sphere_3 Critical_radius;

    Critical_radius critical_radius =
        geom_traits().compute_power_distance_to_power_sphere_3_object();

    return critical_radius(point(c, 0), point(c, 1), point(c, 2), point(c, 3), point(v));
  }

  // \pre c->neighbor(i) is finite
  FT compute_power_distance_to_power_sphere(const Cell_handle c, const int i) const
  {
    Cell_handle nc = c->neighbor(i);
    CGAL_precondition(!this->is_infinite(nc));
    Vertex_handle v = nc->vertex(nc->index(c));

    return compute_power_distance_to_power_sphere(c, v);
  }

  typename Geom_traits::FT min_squared_distance(const Bare_point& p, const Bare_point& q) const
  {
    return geom_traits().compute_squared_distance_3_object()(p, q);
  }
};

// Struct Mesh_triangulation_3
//
template<class MD,
         class K_ = Default,
         class Concurrency_tag_ = Sequential_tag,
         class Vertex_base_ = Default,
         class Cell_base_   = Default>
struct Mesh_triangulation_3
{
private:
  using K = typename Default::Lazy_get<K_, Kernel_traits<MD> >::type;

  using Geom_traits = typename details::Mesh_geom_traits_generator<K>::type;

  using Indices_tuple = Mesh_3::internal::Indices_tuple_t<MD>;
  using Vertex_base = typename Default::Get<
    Vertex_base_,
    Mesh_vertex_generator_3<Geom_traits,
                            Indices_tuple,
                            typename MD::Index> >::type;
  using Cell_base = typename Default::Get<
    Cell_base_,
    Compact_mesh_cell_generator_3<Geom_traits,
                                  typename MD::Subdomain_index,
                                  typename MD::Surface_patch_index,
                                  typename MD::Index> >::type;
  using Concurrency_tag =
      typename Default::Get<Concurrency_tag_, Sequential_tag>::type;
  struct Tds : public Triangulation_data_structure_3<Vertex_base, Cell_base,
                                                     Concurrency_tag> {};
  using Triangulation =
      Mesh_3_regular_triangulation_3_wrapper<Geom_traits, Tds>;
  ;

public:
  using type = Triangulation;
  using Type = type;
};  // end struct Mesh_triangulation_3
}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_TRIANGULATION_3_H
