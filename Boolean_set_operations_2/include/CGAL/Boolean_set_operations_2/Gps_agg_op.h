// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Baruch Zukerman  <baruchzu@post.tau.ac.il>
//             Ophir Setter     <ophir.setter@cs.tau.ac.il>
//             Efi Fogel        <efifogel@gmail.com>

#ifndef CGAL_GPS_AGG_OP_H
#define CGAL_GPS_AGG_OP_H

#include <CGAL/license/Boolean_set_operations_2.h>

/*! \file Gps_agg_op.h
 *
 * The class Gps_agg_op is responsible for aggregated Boolean set operations
 * depending on a visitor template parameter.  It uses the surface-sweep
 * algorithm from the surface-sweep package to overlay all the polygon sets, and
 * then it uses a BFS that determines which of the faces are contained in the
 * result using the visitor.
 */

#include <CGAL/Boolean_set_operations_2/Gps_agg_meta_traits.h>
#include <CGAL/Boolean_set_operations_2/Gps_agg_op_surface_sweep_2.h>
#include <CGAL/Boolean_set_operations_2/Indexed_event.h>
#include <CGAL/Surface_sweep_2/Arr_construction_subcurve.h>

#include <CGAL/Boolean_set_operations_2/Gps_agg_op_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_scanner.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Arr_accessor.h>
#include <CGAL/iterator.h>

namespace CGAL {

template <typename Arrangement_, typename BfsVisitor, template <typename, typename> class SweepVisitor>
class Gps_agg_op {
  using Arrangement_2 = Arrangement_;
  using Bfs_visitor = BfsVisitor;

  using Geometry_traits_2 = typename Arrangement_2::Traits_adaptor_2;
  using Topology_traits = typename Arrangement_2::Topology_traits;

  using Arr = Arrangement_2;
  using Gt2 = Geometry_traits_2;
  using Tt = Topology_traits;

  using Curve_const_iterator = typename Gt2::Curve_const_iterator;
  using Mgt2 = Gps_agg_meta_traits<Arr>;
  using Curve_data = typename Mgt2::Curve_data;
  using Meta_X_monotone_curve_2 = typename Mgt2::X_monotone_curve_2;

  using Halfedge_handle = typename Arr::Halfedge_handle;
  using Halfedge_iterator = typename Arr::Halfedge_iterator;
  using Face_handle = typename Arr::Face_handle;
  using Edge_iterator = typename Arr::Edge_iterator;
  using Vertex_handle = typename Arr::Vertex_handle;
  using Allocator = typename Arr::Allocator;

  using Arr_entry = std::pair<Arr*, std::vector<Vertex_handle> *>;

  // We obtain a proper helper type from the topology traits of the arrangement.
  // However, the arrangement is parametrized with the Gt2 geometry traits,
  // while we need the Mgt2 geometry traits (which derives from Gt2).
  // Thus, we rebind the helper.
  // We cannot parameterized the arrangement with the Mgt2 geometry
  // traits to start with, because it extends the curve type with arrangement
  // dependent types. (It is parameterized with the arrangement type.)
  using Event = Indexed_event<Mgt2, Arr, Allocator>;
  using Subcurve = Arr_construction_subcurve<Mgt2, Event, Allocator>;
  using Helper_tmp = typename Tt::template Construction_helper<Event, Subcurve>;
  using Helper = typename Helper_tmp::template rebind<Mgt2, Arr, Event, Subcurve>::other;
  using Visitor = SweepVisitor<Helper, Arr>;
  using Surface_sweep_2 = Gps_agg_op_surface_sweep_2<Arr, Visitor>;

  using Edges_hash = Unique_hash_map<Halfedge_handle, std::size_t>;
  using Faces_hash = Unique_hash_map<Face_handle, std::size_t>;
  using Bfs_scanner = Gps_bfs_scanner<Arr, Bfs_visitor>;

protected:
  Arr* m_arr;
  Mgt2* m_traits;
  Visitor m_visitor;
  Surface_sweep_2 m_surface_sweep;
  Edges_hash m_edges_hash;      // maps halfedge to its BC (boundary counter)
  Faces_hash m_faces_hash;      // maps face to its IC (inside count)

public:
  /*! constructs. */
  Gps_agg_op(Arr& arr, std::vector<Vertex_handle>& vert_vec, const Gt2& tr) :
    m_arr(&arr),
    m_traits(new Mgt2(tr)),
    m_visitor(&arr, &m_edges_hash, &vert_vec),
    m_surface_sweep(m_traits, &m_visitor)
  {}

  std::pair<std::size_t, std::size_t>
  prepare(std::size_t lower, std::size_t upper, std::size_t jump,
          std::vector<Arr_entry>& arr_vec, std::list<Meta_X_monotone_curve_2>& curves_list) {
    std::size_t n_inf_pgn = 0;  // number of infinite polygons (arrangement
                                // with a contained unbounded face
    std::size_t n_pgn = 0;      // number of polygons (arrangements)
    for (auto i = lower; i <= upper; i += jump, ++n_pgn) {
      // The BFS scan (after the loop) starts in the reference face,
      // so we count the number of polygons that contain the reference face.
      Arr* arr = (arr_vec[i]).first;
      if (arr->reference_face()->contained()) ++n_inf_pgn;

      for (auto itr = arr->edges_begin(); itr != arr->edges_end(); ++itr) {
        // take only relevant edges (which separate between contained and
        // non-contained faces.
        Halfedge_handle he = itr;
        if (he->face()->contained() == he->twin()->face()->contained()) continue;
        if ((Arr_halfedge_direction)he->direction() == ARR_RIGHT_TO_LEFT) he = he->twin();

        Curve_data cv_data(arr, he, 1, 0);
        curves_list.push_back(Meta_X_monotone_curve_2(he->curve(), cv_data));
      }
    }
    return std::make_pair(n_inf_pgn, n_pgn);
  }

  /*! sweeps the plane without interceptions.
   */
  void sweep_arrangements(std::size_t lower, std::size_t upper, std::size_t jump,
                          std::vector<Arr_entry>& arr_vec) {
    std::size_t n_inf_pgn, n_pgn;
    std::list<Meta_X_monotone_curve_2> curves_list;
    std::tie(n_inf_pgn, n_pgn) = prepare(lower, upper, jump, arr_vec, curves_list);
    m_surface_sweep.sweep(curves_list.begin(), curves_list.end(), lower, upper, jump, arr_vec);
    m_faces_hash[m_arr->reference_face()] = n_inf_pgn;
    Bfs_visitor visitor(&m_edges_hash, &m_faces_hash, n_pgn);
    visitor.visit_ubf(m_arr->faces_begin(), n_inf_pgn);
    Bfs_scanner scanner(visitor);
    scanner.scan(*m_arr);
    visitor.after_scan(*m_arr);
  }

  /*! sweeps the plane without interceptions, but stop when an intersection occurs.
   */
  bool sweep_intercept_arrangements(std::size_t lower, std::size_t upper, std::size_t jump,
                                    std::vector<Arr_entry>& arr_vec) {
    std::size_t n_inf_pgn, n_pgn;
    std::list<Meta_X_monotone_curve_2> curves_list;
    std::tie(n_inf_pgn, n_pgn) = prepare(lower, upper, jump, arr_vec, curves_list);
    auto res = m_surface_sweep.sweep_intercept(curves_list.begin(), curves_list.end(), lower, upper, jump, arr_vec);
    if (res) return true;

    m_faces_hash[m_arr->reference_face()] = n_inf_pgn;
    Bfs_visitor visitor(&m_edges_hash, &m_faces_hash, n_pgn);
    visitor.visit_ubf(m_arr->faces_begin(), n_inf_pgn);
    Bfs_scanner scanner(visitor);
    scanner.scan(*m_arr);
    visitor.after_scan(*m_arr);
    return false;
  }

  template <typename InputIterator>
  std::size_t prepare2(InputIterator begin, InputIterator end, std::list<Meta_X_monotone_curve_2>& curves_list) {
    std::size_t n_inf_pgn = 0;  // number of infinite polygons (arrangement
                                // with a contained unbounded face
    for (auto it = begin; it != end; ++it) {
      // The BFS scan (after the loop) starts in the reference face,
      // so we count the number of polygons that contain the reference face.
      Arr* arr = it->first;
      if (arr->reference_face()->contained()) ++n_inf_pgn;

      for (auto ite = arr->edges_begin(); ite != arr->edges_end(); ++ite) {
        // take only relevant edges (which separate between contained and
        // non-contained faces.
        Halfedge_handle he = ite;
        if (he->face()->contained() == he->twin()->face()->contained()) continue;
        if ((Arr_halfedge_direction)he->direction() == ARR_RIGHT_TO_LEFT) he = he->twin();

        Curve_data cv_data(arr, he, 1, 0);
        curves_list.push_back(Meta_X_monotone_curve_2(he->curve(), cv_data));
      }
    }
    return n_inf_pgn;
  }

  /*! sweeps the plane without interceptions, but stop when an intersection occurs.
   */
  template <typename InputIterator>
  bool sweep_intercept_arrangements2(InputIterator begin, InputIterator end) {
    std::list<Meta_X_monotone_curve_2> curves_list;
    auto n_inf_pgn = prepare2(begin, end, curves_list);
    auto res = m_surface_sweep.sweep_intercept2(curves_list.begin(), curves_list.end(), begin, end);
    if (res) return true;

    m_faces_hash[m_arr->reference_face()] = n_inf_pgn;
    std::size_t n_pgn = std::distance(begin, end);      // number of polygons (arrangements)
    Bfs_visitor visitor(&m_edges_hash, &m_faces_hash, n_pgn);
    visitor.visit_ubf(m_arr->faces_begin(), n_inf_pgn);
    Bfs_scanner scanner(visitor);
    scanner.scan(*m_arr);
    visitor.after_scan(*m_arr);
    return false;
  }

  /*! destructs.
   */
  ~Gps_agg_op() { delete m_traits; }
};

} //namespace CGAL

#endif
