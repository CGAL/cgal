// Copyright (c) 2010, 2017 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s)     : Mikhail Bogdanov
//                 Aymeric Pellé
//                 Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_3_MESH_TRIANGULATION_3_H
#define CGAL_PERIODIC_3_MESH_TRIANGULATION_3_H

#include <CGAL/license/Periodic_3_mesh_3.h>

#include <CGAL/Periodic_3_mesh_3/config.h>

// traits class
#include <CGAL/Kernel_traits.h>
#include <CGAL/Robust_weighted_circumcenter_filtered_traits_3.h>
#include <CGAL/internal/Robust_periodic_weighted_circumcenter_traits_3.h>

// periodic triangulations
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Periodic_3_regular_triangulation_3.h>

// vertex and cell bases
#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Mesh_cell_base_3.h>
#include <CGAL/Mesh_3/io_signature.h>

#include <CGAL/assertions.h>
#include <CGAL/array.h>
#include <CGAL/tags.h>

#include <iostream>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>

namespace CGAL {

/// This class currently provides an interface between the classe
/// `CGAL::Periodic_3_regular_triangulation_3` and the mesher `Mesh_3`.
/// As periodic triangulations are parallelized, a lot of these functions will
/// become obsolete.
template<class Gt, class Tds>
class Periodic_3_regular_triangulation_3_mesher_3
  : public Periodic_3_regular_triangulation_3<Gt, Tds>
{
public:
  typedef Sequential_tag                                      Concurrency_tag;
  typedef void                                                Lock_data_structure;

  void *get_lock_data_structure() const { return 0; }
  void set_lock_data_structure(void *) const { }

  typedef Periodic_3_regular_triangulation_3<Gt, Tds>         Base;

  typedef Gt                                                  Geom_traits;
  typedef Tds                                                 Triangulation_data_structure;

  typedef typename Base::Cell_iterator                        Finite_cells_iterator;
  typedef typename Base::Facet_iterator                       Finite_facets_iterator;
  typedef typename Base::Edge_iterator                        Finite_edges_iterator;
  typedef typename Base::Vertex_iterator                      Finite_vertices_iterator;

  typedef typename Base::Vertex_handle                        Vertex_handle;
  typedef typename Base::Edge                                 Edge;
  typedef typename Base::Facet                                Facet;
  typedef typename Base::Cell_handle                          Cell_handle;

  typedef typename Base::FT                                   FT;

  typedef typename Base::Bare_point                           Bare_point;
  typedef typename Base::Weighted_point                       Weighted_point;
  typedef typename Base::Periodic_bare_point                  Periodic_bare_point;
  typedef typename Base::Periodic_weighted_point              Periodic_weighted_point;

  typedef typename Base::Locate_type                          Locate_type;

  typedef typename Base::Segment                              Segment;
  typedef typename Base::Periodic_segment                     Periodic_segment;

  typedef typename Base::Triangle                             Triangle;
  typedef typename Base::Periodic_triangle                    Periodic_triangle;
  typedef typename Base::Tetrahedron                          Tetrahedron;
  typedef typename Base::Periodic_tetrahedron                 Periodic_tetrahedron;

  typedef typename Base::Offset                               Offset;
  typedef typename Base::Iso_cuboid                           Iso_cuboid;
  typedef typename Base::Conflict_tester                      Conflict_tester;
  typedef typename Base::Covering_sheets                      Covering_sheets;

  typedef typename Gt::Vector_3                               Vector_3;
  typedef typename Gt::Ray_3                                  Ray;

  using Base::construct_point;
  using Base::construct_weighted_point;
  using Base::construct_segment;
  using Base::construct_triangle;
  using Base::construct_tetrahedron;
  using Base::construct_periodic_point;
  using Base::construct_periodic_weighted_point;
  using Base::construct_periodic_segment;
  using Base::construct_periodic_triangle;
  using Base::construct_periodic_tetrahedron;
  using Base::domain;
  using Base::dual;
  using Base::get_offset;
  using Base::geom_traits;
  using Base::adjacent_vertices;
  using Base::incident_cells;
  using Base::incident_edges;
  using Base::incident_facets;
  using Base::insert_dummy_points;
  using Base::number_of_vertices;
  using Base::periodic_triangle;
  using Base::periodic_tetrahedron;
  using Base::point;
  using Base::set_point;
  using Base::tds;
  using Base::set_offsets;
#ifndef CGAL_NO_STRUCTURAL_FILTERING
  using Base::inexact_locate;
#endif

  static std::string io_signature() { return Get_io_signature<Base>()(); }

  /// Constructor
  Periodic_3_regular_triangulation_3_mesher_3(const Iso_cuboid& domain = Iso_cuboid(0,0,0,1,1,1),
                                              const Geom_traits& gt = Geom_traits())
    : Base(domain, gt)
  {
    insert_dummy_points();
  }

  /// Concurrency related
  template <typename Cell_handle>
  bool try_lock_cell(const Cell_handle &, int = 0) const { return true; }

  bool try_lock_and_get_incident_cells(Vertex_handle /* v */,
                                       std::vector<Cell_handle>& /* cells */) const
  {
    std::cerr << "ERROR: implement try_lock_and_get_incident_cells()"<< std::endl;
    return true;
  }

  /// Basic setters-getters
  // there are no infinite elements in a periodic triangulation
  bool is_infinite(const Vertex_handle) const { return false; }
  bool is_infinite(const Edge&) const { return false; }
  bool is_infinite(const Facet&) const { return false; }
  bool is_infinite(const Cell_handle) const { return false; }
  bool is_infinite(const Cell_handle, int) const { return false; }
  bool is_infinite(const Cell_handle c, int i, int j) const;

  Cell_handle infinite_cell() const
  {
    // there are no infinite cells in a periodic triangulation
    CGAL_assertion(false);
    return Cell_handle();
  }

  Vertex_handle infinite_vertex() const
  {
    // there is no infinite vertex in a periodic triangulation
    CGAL_assertion(false);
    return Vertex_handle();
  }

  int dimension() const
  {
    // there can be no degenerate dimensions in a periodic triangulation
    return (number_of_vertices() == 0) ? -2 : 3;
  }

  void set_domain(const Iso_cuboid& domain)
  {
    Base::set_domain(domain);
    insert_dummy_points();
  }

  Bare_point snap_to_domain_border(const Bare_point& p) const
  {
    const FT px = p.x();
    const FT py = p.y();
    const FT pz = p.z();
    const FT dxm = domain().xmin();
    const FT dym = domain().ymin();
    const FT dzm = domain().zmin();
    const FT dxM = domain().xmax();
    const FT dyM = domain().ymax();
    const FT dzM = domain().zmax();
    FT sx = px, sy = py, sz = pz;

    // simply comparing to FT::epsilon() is probably not completely satisfactory
    const FT eps = std::numeric_limits<FT>::epsilon();

    if(CGAL::abs(px - dxm) < eps) sx = domain().xmin();
    if(CGAL::abs(px - dxM) < eps) sx = domain().xmax();
    if(CGAL::abs(py - dym) < eps) sy = domain().ymin();
    if(CGAL::abs(py - dyM) < eps) sy = domain().ymax();
    if(CGAL::abs(pz - dzm) < eps) sz = domain().zmin();
    if(CGAL::abs(pz - dzM) < eps) sz = domain().zmax();

    std::cout << "snapped " << p << " to " << sx << " " << sy << " " << sz << std::endl;
    return geom_traits().construct_point_3_object()(sx, sy, sz);
  }

  Weighted_point snap_to_domain_border(const Weighted_point& p) const
  {
    const Bare_point snapped_p = snap_to_domain_border(
                                   geom_traits().construct_point_3_object()(p));
    return geom_traits().construct_weighted_point_3_object()(snapped_p, p.weight());
  }

  /// transform a bare point (living anywhere in space) into the canonical
  /// instance of the same bare point that lives inside the base domain
  Bare_point robust_canonicalize_point(const Bare_point& p) const
  {
    bool had_to_use_exact = false;
    Periodic_bare_point pbp = construct_periodic_point(p, had_to_use_exact);

    if(had_to_use_exact)
    {
      // the point is close to a border, snap it !
      Bare_point sp = snap_to_domain_border(p);

      // might have snapped to a 'max' of the domain, which is not in the domain
      // note: we could snap to 'min' all the time in 'snap_to_domain_border'
      // but this is clearer like that (and costs very little since we should
      // not have to use exact computations too often)
      return canonicalize_point(sp);
    }

    return construct_point(pbp);
  }

  Bare_point canonicalize_point(const Bare_point& p) const
  {
    if(p.x() >= domain().xmin() && p.x() < domain().xmax() &&
       p.y() >= domain().ymin() && p.y() < domain().ymax() &&
       p.z() >= domain().zmin() && p.z() < domain().zmax())
      return p;

    return robust_canonicalize_point(p);
  }

  /// transform a weighted point (living anywhere in space) into the canonical
  /// instance of the same weighted point that lives inside the base domain
  Weighted_point robust_canonicalize_point(const Weighted_point& p) const
  {
    bool had_to_use_exact = false;
    Periodic_weighted_point pwp = construct_periodic_weighted_point(p, had_to_use_exact);

    if(had_to_use_exact)
    {
      // the point is close to a border, snap it !
      Weighted_point sp = snap_to_domain_border(p);
      return canonicalize_point(sp);
    }

    return construct_weighted_point(pwp);
  }

  Weighted_point canonicalize_point(const Weighted_point& p) const
  {
    if(p.x() >= domain().xmin() && p.x() < domain().xmax() &&
       p.y() >= domain().ymin() && p.y() < domain().ymax() &&
       p.z() >= domain().zmin() && p.z() < domain().zmax())
      return p;

    return robust_canonicalize_point(p);
  }

  Triangle triangle(const Facet& f) const
  {
    Periodic_triangle ptri = periodic_triangle(f);
    return construct_triangle(ptri);
  }

  Tetrahedron tetrahedron(const Cell_handle c) const
  {
    Periodic_tetrahedron ptet = periodic_tetrahedron(c);
    return construct_tetrahedron(ptet);
  }

  /*!
  Copies all finite `Edge`s incident to `v` to the
  output iterator `edges`. Returns the resulting output iterator.

  Since there are no infinite edges in a periodic triangulation, this
  function is simply a wrapper around `incident_edges`

  \pre `t.dimension() > 0`, `v != Vertex_handle()`, `t.is_vertex(v)`.
  */
  template<class OutputIterator>
  OutputIterator
  finite_incident_edges(Vertex_handle v, OutputIterator edges) const
  {
    return incident_edges(v, edges);
  }

  /*!
  Copies the `Cell_handle`s of all finite cells incident to `v` to the output
  iterator `cells`.
  Returns the resulting output iterator.

  Since there are no infinite cells in a periodic triangulation, this
  function is simply a wrapper around `incident_cells`

  \pre `t.dimension() == 3`, `v != Vertex_handle()`, `t.is_vertex(v)`.
  */
  template<class OutputIterator>
  OutputIterator
  finite_incident_cells(Vertex_handle v, OutputIterator cells) const
  {
    return incident_cells(v, cells);
  }

  /*!
  Copies all finite `Facet`s incident to `v` to the output iterator
  `facets`.
  Returns the resulting output iterator.

  Since there are no infinite facets in a periodic triangulation, this
  function is simply a wrapper around `incident_facets`

  \pre `t.dimension() > 1`, `v != Vertex_handle()`, `t.is_vertex(v)`.
  */
  template<class OutputIterator>
  OutputIterator
  finite_incident_facets(Vertex_handle v, OutputIterator facets) const
  {
    return incident_facets(v, facets);
  }

  // Periodic triangulations cannot be used in parallel (yet), but the functions below
  // are required for compilation. Note that these functions could already
  // be implemented and moved to P3T3.h but to make it clear that P3M3 is
  // not available in parallel, they are put here and left empty.

  template <class OutputIterator>
  OutputIterator
  incident_edges_threadsafe(Vertex_handle /* v */, OutputIterator edges) const
  {
    CGAL_assertion(false); // not yet supported
    return edges;
  }

  template <class OutputIterator>
  OutputIterator
  incident_facets_threadsafe(Vertex_handle /* v */, OutputIterator facets) const
  {
    CGAL_assertion(false); // not yet supported
    return facets;
  }

  template <typename OutputIterator>
  void
  incident_cells_threadsafe(Vertex_handle /* v */, OutputIterator /* cells */) const
  {
    CGAL_assertion(false); // not yet supported
  }

  void clear_v_offsets() const
  {
    for (typename std::vector<Vertex_handle>::iterator voit = this->v_offsets.begin();
         voit != this->v_offsets.end() ; ++voit) {
      (*voit)->clear_offset();
    }
    this->v_offsets.clear();
  }

  FT compute_power_distance_to_power_sphere(const Cell_handle& c, const int i) const
  {
    typename Gt::Compute_power_distance_to_power_sphere_3 cr =
      geom_traits().compute_power_distance_to_power_sphere_3_object();

    Offset o_nb = this->neighbor_offset(c, i, c->neighbor(i));
    Offset o_vt = this->get_offset(c->neighbor(i), c->neighbor(i)->index(c));

    const Weighted_point& wp0 = this->point(c->vertex(0)); // need the canonical point
    const Weighted_point& wp1 = this->point(c->vertex(1));
    const Weighted_point& wp2 = this->point(c->vertex(2));
    const Weighted_point& wp3 = this->point(c->vertex(3));
    const Weighted_point& wq = this->point(c->neighbor(i)->vertex(c->neighbor(i)->index(c)));
    const Offset& op0 = this->get_offset(c, 0);
    const Offset& op1 = this->get_offset(c, 1);
    const Offset& op2 = this->get_offset(c, 2);
    const Offset& op3 = this->get_offset(c, 3);
    const Offset& oq = o_vt - o_nb;

    return cr(wp0, wp1, wp2, wp3, wq, op0, op1, op2, op3, oq);
  }

  // The functions below are needed by Mesh_3 but need a specific implementation
  // for the periodic case because we need to try with different offsets
  // to get a result
  FT compute_power_distance_to_power_sphere(const Cell_handle& c,
                                            const Vertex_handle v) const
  {
    typename Gt::Compute_power_distance_to_power_sphere_3 cr =
      geom_traits().compute_power_distance_to_power_sphere_3_object();

    FT min_power_dist = std::numeric_limits<FT>::infinity();

    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        for(int k = 0; k < 3; ++k) {
          const Weighted_point& wp0 = this->point(c->vertex(0)); // need the canonical point
          const Weighted_point& wp1 = this->point(c->vertex(1));
          const Weighted_point& wp2 = this->point(c->vertex(2));
          const Weighted_point& wp3 = this->point(c->vertex(3));
          const Weighted_point& wq = this->point(v);
          const Offset& op0 = this->get_offset(c, 0);
          const Offset& op1 = this->get_offset(c, 1);
          const Offset& op2 = this->get_offset(c, 2);
          const Offset& op3 = this->get_offset(c, 3);
          const Offset oq(i-1, j-1, k-1);

          FT power_dist = cr(wp0, wp1, wp2, wp3, wq, op0, op1, op2, op3, oq);

          if(power_dist < min_power_dist)
            min_power_dist = power_dist;
        }
      }
    }

    return min_power_dist;
  }

  // Return the tetrahedron made of 'f' + 'wp'
  // \pre there exists an offset such that 'f.first' is in conflict with 'wp'
  Tetrahedron tetrahedron(const Facet& f, const Weighted_point& wp) const
  {
    Weighted_point canonic_wp = canonicalize_point(wp);
    Conflict_tester tester(canonic_wp, this);

    const Cell_handle c = f.first;
    const int index = f.second;

    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        for(int k = 0; k < 3; ++k) {
          const Offset off(i-1, j-1, k-1);
          if(tester(c, off))
          {
            return construct_tetrahedron(
                       canonic_wp, this->point(c->vertex((index+1)&3)),
                       this->point(c->vertex((index+2)&3)), this->point(c->vertex((index+3)&3)),
                       off, this->get_offset(c, (index+1)&3),
                       this->get_offset(c, (index+2)&3), this->get_offset(c, (index+3)&3));
          }
        }
      }
    }

    CGAL_assertion(false);
    return Tetrahedron();
  }

  Bounded_side side_of_power_sphere(const Cell_handle& c, const Weighted_point& p,
                                    bool perturb = false) const
  {
    Weighted_point canonical_p = canonicalize_point(p);

//    std::cout << "Pside_of_power_sphere with at " << &*c << std::endl
//              << "                              " << &*(c->vertex(0)) << " : " << c->vertex(0)->point() << std::endl
//              << "                              " << &*(c->vertex(1)) << " : " << c->vertex(1)->point() << std::endl
//              << "                              " << &*(c->vertex(2)) << " : " << c->vertex(2)->point() << std::endl
//              << "                              " << &*(c->vertex(3)) << " : " << c->vertex(3)->point() << std::endl
//              << " Foreign: " << p << std::endl;

    Bounded_side bs = ON_UNBOUNDED_SIDE;
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        for(int k = 0; k < 3; ++k) {
          bs = Base::_side_of_power_sphere(c, canonical_p,
                                           Offset(i-1, j-1, k-1), perturb);

          if(bs != ON_UNBOUNDED_SIDE)
            return bs;
        }
      }
    }

    return bs;

    return Base::side_of_power_sphere(c, canonical_p, Offset(), perturb);
  }

  // Warning : This is a periodic version that computes the smallest possible
  // between 'p' and 'q', for all possible combinations of offsets
  FT min_squared_distance(const Bare_point& p, const Bare_point& q) const
  {
    typename Geom_traits::Compute_squared_distance_3 csd =
      geom_traits().compute_squared_distance_3_object();

    const Bare_point cp = canonicalize_point(p);
    const Bare_point cq = canonicalize_point(q);

//    std::cout << "minsqd: " << p << " // " << q << std::endl;
//    std::cout << "canon: " << cp << " // " << cq << std::endl;

    FT min_sq_dist = std::numeric_limits<FT>::infinity();

    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        for(int k = 0; k < 3; ++k) {
          FT sq_dist = csd(cq, construct_point(std::make_pair(cp, Offset(i-1, j-1, k-1))));

          if(sq_dist < min_sq_dist)
            min_sq_dist = sq_dist;
        }
      }
    }

//    std::cout << "minsqdt: " << min_sq_dist << std::endl;
    return min_sq_dist;
  }

  // Warning : This function finds which offset 'Oq' should be applied to 'q' so
  // that the distance between 'p' and '(q, Oq)' is minimal.
  //
  // \pre 'p' lives in the canonical domain.
  Bare_point get_closest_point(const Bare_point& p, const Bare_point& q) const
  {
    Bare_point rq;
    const Bare_point cq = canonicalize_point(q);
    FT min_sq_dist = std::numeric_limits<FT>::infinity();

    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        for(int k = 0; k < 3; ++k) {
          const Bare_point tcq = construct_point(std::make_pair(cq, Offset(i-1, j-1, k-1)));
          FT sq_dist = geom_traits().compute_squared_distance_3_object()(p, tcq);

          if(sq_dist < min_sq_dist)
          {
            rq = tcq;
            min_sq_dist = sq_dist;
          }
        }
      }
    }

    return rq;
  }

  // Warning: This is a periodic version that computes the smallest possible
  // distances between p and q, and between p and r FOR ALL POSSIBLE OFFSETS
  // before comparing these distances.
  //
  // It is used in facet encroachement checks in Periodic_3_mesh_3.
  bool greater_or_equal_power_distance(const Bare_point& p,
                                       const Weighted_point& q,
                                       const Weighted_point& r) const
  {
    typename Geom_traits::Compute_power_product_3 power_distance =
      geom_traits().compute_power_product_3_object();

    // canonicalize the points
    const Weighted_point cp =
      geom_traits().construct_weighted_point_3_object()(canonicalize_point(p));
    const Weighted_point cq = canonicalize_point(q);
    const Weighted_point cr = canonicalize_point(r);

    FT min_power_distance_to_q = std::numeric_limits<FT>::infinity();
    FT min_power_distance_to_r = std::numeric_limits<FT>::infinity();

    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        for(int k = 0; k < 3; ++k) {
          const Weighted_point cp_copy =
            construct_weighted_point(std::make_pair(cp, Offset(i-1, j-1, k-1)));

          const FT power_distance_to_q = power_distance(cp_copy, cq);
          if(power_distance_to_q < min_power_distance_to_q)
            min_power_distance_to_q = power_distance_to_q;

          const FT power_distance_to_r = power_distance(cp_copy, cr);
          if (power_distance_to_r < min_power_distance_to_r)
            min_power_distance_to_r = power_distance_to_r;
        }
      }
    }

    CGAL_postcondition(min_power_distance_to_r < 0.5 && min_power_distance_to_q < 0.5);
    return min_power_distance_to_q >= min_power_distance_to_r;
  }

  /// \name Locate functions
  ///
  /// Locate points within a periodic triangulation.
  ///
  /// These functions are temporarily here to interface between Mesh_3 and
  /// the periodic triangulations, until the latter are made parallel.
  ///
  /// \sa `CGAL::Regular_triangulation_3::locate`
  /// @{
  Vertex_handle nearest_power_vertex(const Bare_point& p, Cell_handle start) const
  {
    return Base::nearest_power_vertex(canonicalize_point(p), start);
  }

  /// Return the squared distance (note: _NOT_ the power distance) between the
  /// 'p' and the closest vertex for the power distance.
  std::pair<Vertex_handle, FT>
  nearest_power_vertex_with_sq_distance(const Bare_point& p, Cell_handle start) const
  {
    // The function below is almost a copy from 'nearest_power_vertex' in P3RT3.
    // Any change should be mirrored.
    CGAL_precondition(number_of_vertices() > 0);

    Bare_point canonical_p = canonicalize_point(p);

    Locate_type lt;
    int li, lj;

    typename Gt::Compute_squared_distance_3 csd =
      geom_traits().compute_squared_distance_3_object();
    typename Gt::Construct_point_3 cp =
      geom_traits().construct_point_3_object();
    typename Gt::Construct_weighted_point_3 cwp =
      geom_traits().construct_weighted_point_3_object();

    Offset query_offset;
    Cell_handle c = Base::locate(cwp(canonical_p), query_offset, lt, li, lj, start);

    // - start with the closest vertex from the located cell.
    // - repeatedly take the nearest of its incident vertices if any
    // - if not, we're done.

    // Take the opposite because periodic_locate() returns the offset such that
    // cell + offset contains 'p' but here we need to move 'p'
    query_offset = this->combine_offsets(Offset(), -query_offset);

    Vertex_handle nearest = Base::nearest_vertex_in_cell(c, canonical_p, query_offset);
    const Weighted_point& nearest_wp = this->point(nearest);
    Offset offset_of_nearest = Base::get_min_dist_offset(canonical_p, query_offset, nearest);
    FT min_sq_dist = csd(canonical_p, cp(nearest_wp), query_offset, offset_of_nearest);

    std::cout << "seeking nearest to " << canonical_p << std::endl;
    std::cout << "query offset: " << query_offset << std::endl;
    std::cout << "so point: " << this->construct_point(canonical_p, query_offset) << std::endl;
    std::cout << "nearest: " << nearest->point() << std::endl;
    std::cout << "initial min offset: " << offset_of_nearest << std::endl;
    std::cout << "so point: " << this->construct_point(nearest->point(), offset_of_nearest) << std::endl;
    std::cout << "giving a distance of: " << min_sq_dist << std::endl;

    std::vector<Vertex_handle> vs;
    vs.reserve(32);

    while(true)
    {
      Vertex_handle tmp = nearest;

      adjacent_vertices(nearest, std::back_inserter(vs));
      for(typename std::vector<Vertex_handle>::const_iterator vsit = vs.begin();
                                                              vsit != vs.end(); ++vsit)
      {
        // Can happen in 27-sheeted triangulations composed of few points
        if((*vsit)->point() == nearest->point())
          continue;

        const Offset min_dist_offset = this->get_min_dist_offset(canonical_p, query_offset, *vsit);
        if(this->compare_distance(canonical_p, (*vsit)->point(), tmp->point(),
                                  query_offset, min_dist_offset, offset_of_nearest) == SMALLER)
        {
          tmp = *vsit;
          offset_of_nearest = min_dist_offset;
          const Weighted_point& vswp = this->point(tmp);
          min_sq_dist = csd(canonical_p, cp(vswp), query_offset, min_dist_offset);
          std::cout << "new closest: " << tmp->point()
                    << " offset: " << offset_of_nearest << std::endl
                    << " sq dist: " << min_sq_dist << std::endl;
        }
      }

      if(tmp == nearest)
        break;

      vs.clear();
      nearest = tmp;
    }

    return std::make_pair(nearest, min_sq_dist);
  }

  Cell_handle locate(const Weighted_point& p,
                     Cell_handle start = Cell_handle(),
                     bool* CGAL_assertion_code(could_lock_zone) = NULL) const
  {
    CGAL_assertion(could_lock_zone == NULL);
    return Base::locate(canonicalize_point(p), start);
  }

  Cell_handle locate(const Weighted_point& p,
                     Vertex_handle hint,
                     bool* CGAL_assertion_code(could_lock_zone) = NULL) const
  {
    CGAL_assertion(could_lock_zone == NULL);
    // Compared to the non-periodic version in T3, the infinite cell cannot
    // be used as default hint, so `Cell_handle()` is used instead.
    return Base::locate(canonicalize_point(p),
                        hint == Vertex_handle() ? Cell_handle() : hint->cell());
  }

  Cell_handle locate(const Weighted_point& p,
                     Locate_type& l, int& i, int& j,
                     Cell_handle start = Cell_handle(),
                     bool* CGAL_assertion_code(could_lock_zone) = NULL) const
  {
    CGAL_assertion(could_lock_zone == NULL);
    return Base::locate(canonicalize_point(p), l, i, j, start);
  }

  Cell_handle locate(const Weighted_point& p,
                     Locate_type& l, int& i, int& j,
                     Vertex_handle hint,
                     bool* CGAL_assertion_code(could_lock_zone) = NULL) const
  {
    CGAL_assertion(could_lock_zone == NULL);
    return Base::locate(canonicalize_point(p), l, i, j,
                        hint == Vertex_handle() ? Cell_handle() : hint->cell());
  }
  /// @}

  /// \name Conflict functions
  /// Returns the vertices on the interior of the conflict hole.
  ///
  /// These functions are temporarily here to interface between Mesh_3 and
  /// the periodic triangulations, until the latter are made parallel.
  ///
  /// @{
  template <class OutputIterator>
  OutputIterator
  vertices_inside_conflict_zone(const Weighted_point& /* p */,
                                Cell_handle /* c */,
                                OutputIterator res) const
  {
    return res;
    CGAL_assertion(false); // not yet supported
  }

  //template < class Gt, class Tds >
  template <class OutputIteratorBoundaryFacets, class OutputIteratorCells,
            class OutputIteratorInternalFacets>
  Triple<OutputIteratorBoundaryFacets, OutputIteratorCells, OutputIteratorInternalFacets>
  find_conflicts(const Weighted_point& p,
                 Cell_handle c,
                 OutputIteratorBoundaryFacets bfit,
                 OutputIteratorCells cit,
                 OutputIteratorInternalFacets ifit,
                 bool* CGAL_assertion_code(could_lock_zone) = NULL,
                 const Facet* /* this_facet_must_be_in_the_cz */ = NULL,
                 bool* /* the_facet_is_in_its_cz */ = NULL) const
  {
    CGAL_triangulation_precondition(could_lock_zone == NULL);
    CGAL_triangulation_precondition(number_of_vertices() != 0);

    clear_v_offsets();

    Weighted_point canonic_p = canonicalize_point(p);

    Locate_type lt;
    int li, lj;
    c = locate(canonic_p, lt, li, lj, c);

    std::vector<Facet> facets;
    facets.reserve(64);
    std::vector<Cell_handle> cells;
    cells.reserve(32);

    Conflict_tester tester(canonic_p, this);
    Triple<typename std::back_insert_iterator<std::vector<Facet> >,
           typename std::back_insert_iterator<std::vector<Cell_handle> >,
           OutputIteratorInternalFacets> tit =
             Base::find_conflicts(c, tester,
                                  make_triple(std::back_inserter(facets),
                                              std::back_inserter(cells), ifit));
    ifit = tit.third;

    // Reset the conflict flag on the boundary.
    for(typename std::vector<Facet>::iterator fit=facets.begin();
        fit != facets.end(); ++fit) {
      fit->first->neighbor(fit->second)->tds_data().clear();
      *bfit++ = *fit;
    }

    // Reset the conflict flag in the conflict cells.
    for(typename std::vector<Cell_handle>::iterator ccit=cells.begin();
        ccit != cells.end(); ++ccit) {
      (*ccit)->tds_data().clear();
      *cit++ = *ccit;
    }

    return make_triple(bfit, cit, ifit);
  }

  template <class OutputIteratorBoundaryFacets, class OutputIteratorCells>
  std::pair<OutputIteratorBoundaryFacets, OutputIteratorCells>
  find_conflicts(const Weighted_point &p, Cell_handle c,
                 OutputIteratorBoundaryFacets bfit,
                 OutputIteratorCells cit,
                 bool* could_lock_zone = NULL) const
  {
    CGAL_assertion(could_lock_zone == NULL);

    Triple<OutputIteratorBoundaryFacets,
           OutputIteratorCells,
           Emptyset_iterator> t = find_conflicts(p, c, bfit, cit,
                                                 Emptyset_iterator(),
                                                 could_lock_zone);
    return std::make_pair(t.first, t.second);
  }
  /// @}

  void set_point(const Vertex_handle v,
                 const Vector_3& move,
                 const Weighted_point& new_position)
  {
    // calling robust canonical here means we don't necessarily have
    // canonical(v + move) = new_position... @fixme
    return Base::set_point(v, move, canonicalize_point(new_position));
  }

  /// \name Insert functions
  ///
  /// Insert points in the triangulation.
  ///
  /// These functions are temporarily here to interface between Mesh_3 and
  /// the periodic triangulations, until the latter are made parallel.
  ///
  /// @{
  template <class CellIt>
  Vertex_handle insert_in_hole(const Weighted_point& p,
                               CellIt cell_begin, CellIt cell_end,
                               Cell_handle begin, int i)
  {
    Vertex_handle v = tds().insert_in_hole(cell_begin, cell_end, begin, i);
    v->set_point(canonicalize_point(p));

    std::vector<Cell_handle> incident_cells_;
    incident_cells_.reserve(64);
    incident_cells(v, std::back_inserter(incident_cells_));

    // For all cells incident to the newly added vertex v: fetch their offsets from
    // the tester and reset them in the triangulation data structure.
    typename std::vector<Cell_handle>::iterator cit = incident_cells_.begin(),
                                                cend = incident_cells_.end();
    for (; cit != cend; cit++)
    {
      Offset off[4];
      for (int i=0 ; i<4 ; ++i)
        off[i] = (*cit)->vertex(i)->offset();

      set_offsets(*cit, off[0], off[1], off[2], off[3]);
    }

    clear_v_offsets();

    return v;
  }

  Vertex_handle insert(const Weighted_point& p,
                       Cell_handle start = Cell_handle(),
                       bool* CGAL_assertion_code(could_lock_zone) = NULL)
  {
    CGAL_assertion(could_lock_zone == NULL);
    return Base::insert(canonicalize_point(p), start);
  }

  Vertex_handle insert(const Weighted_point& p,
                       Vertex_handle hint,
                       bool* CGAL_assertion_code(could_lock_zone) = NULL)
  {
    CGAL_assertion(could_lock_zone == NULL);
    // compared to the non-periodic version in T3, the infinite cell cannot
    // be used; `Cell_handle()` is used instead
    return Base::insert(canonicalize_point(p),
                        hint == Vertex_handle() ? Cell_handle() : hint->cell());
  }

  Vertex_handle insert(const Weighted_point& p,
                       Locate_type lt, Cell_handle loc, int li, int lj,
                       bool* CGAL_assertion_code(could_lock_zone) = NULL)
  {
    CGAL_assertion(could_lock_zone == NULL);
    return Base::insert(canonicalize_point(p), lt, loc, li, lj);
  }
  /// @}

  /// Remove function
  void remove(Vertex_handle v,
              bool* CGAL_assertion_code(could_lock_zone) = NULL)
  {
    CGAL_assertion(could_lock_zone == NULL);
    return Base::remove(v);
  }

  /// Dual computations
  Object dual(const Facet & f) const
  {
    Segment s = construct_segment(Base::dual(f));
    return make_object(s);
  }

  void dual_segment(const Facet& facet, Bare_point& p, Bare_point& q) const {
    Segment s = construct_segment(Base::dual(facet));
    p = s.source();
    q = s.target();
    return;
  }

  void dual_segment_exact(const Facet& facet, Bare_point& p, Bare_point& q) const
  {
    // @fixme (below is not exact)
    return dual_segment(facet, p, q);
  }

  // dual rays are impossible in a periodic triangulation since there are no
  // infinite cells, but these functions are still required to compile Mesh_3
  void dual_ray(const Facet& /*f*/, Ray& /*ray*/) const { CGAL_assertion(false); }
  void dual_ray_exact(const Facet& /*facet*/, Ray& /*ray*/) const { CGAL_assertion(false); }
};

namespace details {

template<typename K>
struct Periodic_3_mesh_geom_traits_generator
{
private:
  typedef Robust_periodic_weighted_circumcenter_traits_3<
            Periodic_3_regular_triangulation_traits_3<
              Robust_weighted_circumcenter_filtered_traits_3<K> > > Geom_traits;

public:
  typedef Geom_traits type;
  typedef type Type;
};  // end struct Periodic_3_mesh_geom_traits_generator

}  // end namespace details

template<class MD,
         class K_ = Default,
         class Vertex_base_ = Default,
         class Cell_base_ = Default>
class Periodic_3_mesh_triangulation_3
{
  // default K
  typedef typename Default::Get<K_, typename Kernel_traits<MD>::Kernel>::type K;

  // traits
  typedef typename details::Periodic_3_mesh_geom_traits_generator<K>::type Geom_traits;

  // Periodic vertex and cell bases
  typedef Periodic_3_triangulation_ds_vertex_base_3<> VbDS;
  typedef Regular_triangulation_vertex_base_3<Geom_traits, VbDS> PVb;

  typedef Periodic_3_triangulation_ds_cell_base_3<> CbDS;
  typedef Regular_triangulation_cell_base_3<Geom_traits, CbDS> RCb;
  typedef Regular_triangulation_cell_base_with_weighted_circumcenter_3<Geom_traits, RCb> PCb;

  typedef Mesh_vertex_base_3<Geom_traits, MD, PVb> Default_Vb;
  typedef Mesh_cell_base_3<Geom_traits, MD, PCb> Default_Cb;

  // default Vb/Cb
  typedef typename Default::Get<Vertex_base_, Default_Vb>::type Vertex_base;
  typedef typename Default::Get<Cell_base_, Default_Cb>::type Cell_base;

  // Triangulation and tds
  typedef Triangulation_data_structure_3<Vertex_base, Cell_base> Tds;
  typedef Periodic_3_regular_triangulation_3_mesher_3<Geom_traits, Tds> Triangulation;

public:
  typedef Triangulation type;
  typedef type Type;
};  // end struct Mesh_triangulation_3

} // namespace CGAL

#endif // CGAL_PERIODIC_3_MESH_TRIANGULATION_3_H

