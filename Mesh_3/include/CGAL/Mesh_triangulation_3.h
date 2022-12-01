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
#include <CGAL/SMDS_3/io_signature.h>

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
  typedef typename Base::Facet                                Facet;
  typedef typename Base::Cell_handle                          Cell_handle;

  typedef typename Geom_traits::Vector_3                      Vector;

  using Base::geom_traits;
  using Base::point;
  using Base::triangle;

  static std::string io_signature() { return Get_io_signature<Base>()(); }

  // The undocumented functions below are required for Mesh_3 because of Periodic_3_mesh_3:
  // they are functions for which both triangulations
  // have fundamentally different implementations (for example, Periodic_3_mesh_3
  // might not know the offset of the closest point and must brute-force check for all
  // possibilities). To enable Periodic_3_mesh_3 to use Mesh_3's files,
  // each mesh triangulation implements its own version.

  const Bare_point& get_closest_point(const Bare_point& /*p*/, const Bare_point& q) const
  {
    return q;
  }

  Triangle get_incident_triangle(const Facet& f, const Vertex_handle) const
  {
    return triangle(f);
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

  // Duals
  template < class Gt, class Tds, class Lds >
  void
  Regular_triangulation_3<Gt,Tds,Lds>::
  dual_segment(Cell_handle c, int i, Bare_point& p, Bare_point&q) const
  {
    Cell_handle n = c->neighbor(i);
    CGAL_assertion(! is_infinite(c) && ! is_infinite(n));

    p = dual(c);
    q = dual(n);
  }

  template < class Gt, class Tds, class Lds >
  void
  Regular_triangulation_3<Gt,Tds,Lds>::
  dual_segment(const Facet& facet, Bare_point& p, Bare_point&q) const
  {
    return dual_segment(facet.first, facet.second, p, q);
  }

  template < class Gt, class Tds, class Lds >
  void
  Regular_triangulation_3<Gt,Tds,Lds>::
  dual_ray(Cell_handle c, int i, Ray& ray) const
  {
    Cell_handle n = c->neighbor(i);
    CGAL_precondition((!is_infinite(c) != !is_infinite(n))); // xor
    // either n or c is infinite
    int in;
    if(is_infinite(c))
    {
      in = n->index(c);
    }
    else
    {
      n = c;
      in = i;
    }

    // n now denotes a finite cell, either c or c->neighbor(i)
    int ind[3] = {(in+1)&3,(in+2)&3,(in+3)&3};
    if((in&1) == 1)
      std::swap(ind[0], ind[1]);

    const Weighted_point& p = n->vertex(ind[0])->point();
    const Weighted_point& q = n->vertex(ind[1])->point();
    const Weighted_point& r = n->vertex(ind[2])->point();
    const Bare_point& bp = construct_point(p);
    const Bare_point& bq = construct_point(q);
    const Bare_point& br = construct_point(r);

    Line l = construct_perpendicular_line(construct_plane(bp, bq, br),
                                          construct_weighted_circumcenter(p,q,r));

    ray = construct_ray(dual(n), l);
  }

  template < class Gt, class Tds, class Lds >
  void
  Regular_triangulation_3<Gt,Tds,Lds>::
  dual_ray(const Facet& facet, Ray& ray) const
  {
    return dual_ray(facet.first, facet.second, ray);
  }

  // Exact versions of dual_segment() and dual_ray() for Mesh_3.
  // These functions are really dirty: they assume that the point type is nice enough
  // such that EPECK can manipulate it (e.g. convert it to EPECK::Point_3) AND
  // that the result of these manipulations will make sense.
  template < class Gt, class Tds, class Lds >
  void
  Regular_triangulation_3<Gt,Tds,Lds>::
  dual_segment_exact(const Facet& facet, Bare_point& p, Bare_point&q) const
  {
    typedef typename Kernel_traits<Bare_point>::Kernel           K;
    typedef Exact_predicates_exact_constructions_kernel          EK;
    typedef Cartesian_converter<K, EK>                           To_exact;
    typedef Cartesian_converter<EK,K>                            Back_from_exact;

    typedef EK                                                   Exact_Rt;

    To_exact to_exact;
    Back_from_exact back_from_exact;
    Exact_Rt::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
        Exact_Rt().construct_weighted_circumcenter_3_object();

    Cell_handle c = facet.first;
    int i = facet.second;
    Cell_handle n = c->neighbor(i);

    const typename Exact_Rt::Weighted_point_3& cp = to_exact(c->vertex(0)->point());
    const typename Exact_Rt::Weighted_point_3& cq = to_exact(c->vertex(1)->point());
    const typename Exact_Rt::Weighted_point_3& cr = to_exact(c->vertex(2)->point());
    const typename Exact_Rt::Weighted_point_3& cs = to_exact(c->vertex(3)->point());

    const typename Exact_Rt::Weighted_point_3& np = to_exact(n->vertex(0)->point());
    const typename Exact_Rt::Weighted_point_3& nq = to_exact(n->vertex(1)->point());
    const typename Exact_Rt::Weighted_point_3& nr = to_exact(n->vertex(2)->point());
    const typename Exact_Rt::Weighted_point_3& ns = to_exact(n->vertex(3)->point());

    p = back_from_exact(exact_weighted_circumcenter(cp, cq, cr, cs));
    q = back_from_exact(exact_weighted_circumcenter(np, nq, nr, ns));
  }

  template < class Gt, class Tds, class Lds >
  void
  Regular_triangulation_3<Gt,Tds,Lds>::
  dual_ray_exact(const Facet& facet, Ray& ray) const
  {
    Cell_handle c = facet.first;
    int i = facet.second;
    Cell_handle n = c->neighbor(i);
    CGAL_precondition(!is_infinite(c) != !is_infinite(n)); // xor
    // either n or c is infinite
    int in;
    if(is_infinite(c))
    {
      in = n->index(c);
    }
    else
    {
      n = c;
      in = i;
    }

    // n now denotes a finite cell, either c or c->neighbor(i)
    int ind[3] = {(in+1)&3,(in+2)&3,(in+3)&3};
    if((in&1) == 1)
      std::swap(ind[0], ind[1]);

    // exact part
    typedef typename Kernel_traits<Bare_point>::Kernel           K;
    typedef Exact_predicates_exact_constructions_kernel          EK;
    typedef Cartesian_converter<K, EK>                           To_exact;
    typedef Cartesian_converter<EK,K>                            Back_from_exact;

    typedef EK                                                   Exact_Rt;

    To_exact to_exact;
    Back_from_exact back_from_exact;

    Exact_Rt::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
        Exact_Rt().construct_weighted_circumcenter_3_object();
    Exact_Rt::Construct_perpendicular_line_3 exact_perpendicular_line =
        Exact_Rt().construct_perpendicular_line_3_object();
    Exact_Rt::Construct_plane_3 exact_plane_3 = Exact_Rt().construct_plane_3_object();
    Exact_Rt::Construct_ray_3 exact_ray_3 = Exact_Rt().construct_ray_3_object();
    Exact_Rt::Construct_point_3 exact_point_3 = Exact_Rt().construct_point_3_object();

    const typename Exact_Rt::Weighted_point_3& p = to_exact(n->vertex(ind[0])->point());
    const typename Exact_Rt::Weighted_point_3& q = to_exact(n->vertex(ind[1])->point());
    const typename Exact_Rt::Weighted_point_3& r = to_exact(n->vertex(ind[2])->point());
    const typename Exact_Rt::Weighted_point_3& s = to_exact(n->vertex(in)->point());

    const typename Exact_Rt::Point_3& bp = exact_point_3(p);
    const typename Exact_Rt::Point_3& bq = exact_point_3(q);
    const typename Exact_Rt::Point_3& br = exact_point_3(r);

    typename Exact_Rt::Line_3 l = exact_perpendicular_line(
                                    exact_plane_3(bp, bq, br),
                                    exact_weighted_circumcenter(p, q, r));

    ray = back_from_exact(exact_ray_3(exact_weighted_circumcenter(p, q, r, s), l));
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

public:
  using type = Triangulation;
  using Type = type;
};  // end struct Mesh_triangulation_3
}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_TRIANGULATION_3_H
