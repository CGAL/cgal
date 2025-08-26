// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_BSO_2_GPS_POLYGON_SIMPILFIER_H
#define CGAL_BSO_2_GPS_POLYGON_SIMPILFIER_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/Boolean_set_operations_2/Gps_simplifier_traits.h>
#include <CGAL/Surface_sweep_2.h>
#include <CGAL/Surface_sweep_2/Arr_construction_subcurve.h>
#include <CGAL/Surface_sweep_2/Arr_construction_event.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_scanner.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_join_visitor.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Arr_accessor.h>
#include <CGAL/iterator.h>

namespace CGAL {

namespace Ss2 = Surface_sweep_2;

template <typename Arrangement_>
class Gps_polygon_simplifier {
  using Arrangement_2 = Arrangement_;

  using Geometry_traits_2 = typename Arrangement_2::Geometry_traits_2;
  using Topology_traits = typename Arrangement_2::Topology_traits;

  using Arr = Arrangement_2;
  using Gt2 = Geometry_traits_2;
  using Tt = Topology_traits;

  using Curve_const_iterator = typename Gt2::Curve_const_iterator;
  using Polygon_2 = typename Gt2::Polygon_2;
  using Polygon_with_holes_2 = typename Gt2::Polygon_with_holes_2;
  using Construct_curves_2 = typename Gt2::Construct_curves_2;

  using Mgt2 = Gps_simplifier_traits<Gt2>;
  using Curve_data = typename Mgt2::Curve_data;
  using Meta_X_monotone_curve_2 = typename Mgt2::X_monotone_curve_2;

  using Halfedge_handle = typename Arr::Halfedge_handle;
  using Halfedge_iterator = typename Arr::Halfedge_iterator;
  using Face_handle = typename Arr::Face_handle;
  using Face_iterator = typename Arr::Face_iterator;
  using Edge_iterator = typename Arr::Edge_iterator;
  using Vertex_handle = typename Arr::Vertex_handle;
  using Ccb_halfedge_const_circulator = typename Arr::Ccb_halfedge_const_circulator;
  using Ccb_halfedge_circulator = typename Arr::Ccb_halfedge_circulator;
  using Allocator = typename Arr::Allocator;

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
  using Visitor = Gps_agg_op_base_visitor<Helper, Arr>;
  using Surface_sweep_2 = Ss2::Surface_sweep_2<Visitor>;

  using Edges_hash = Unique_hash_map<Halfedge_handle, std::size_t>;

  using Faces_hash = Unique_hash_map<Face_handle, std::size_t>;
  using Bfs_visitor = Gps_bfs_join_visitor<Arr>;
  using Bfs_scanner = Gps_bfs_scanner<Arr, Bfs_visitor>;

protected:
  Arr* m_arr;
  const Mgt2* m_traits;
  bool m_own_traits;
  Visitor m_visitor;
  Surface_sweep_2 m_surface_sweep;
  Edges_hash m_edges_hash;      // maps halfedge to its BC (boundary counter)
  Faces_hash m_faces_hash;      // maps face to its IC (inside count)

public:
   /*! Construct. */
  Gps_polygon_simplifier(Arr& arr, const Gt2& tr) :
    m_arr(&arr),
    m_traits(new Mgt2(tr)),
    m_own_traits(true),
    m_visitor(&arr, &m_edges_hash),
    m_surface_sweep(m_traits, &m_visitor)
  {}

  /*! Destructor. */
  ~Gps_polygon_simplifier() {
    if (m_own_traits && (m_traits != nullptr)) {
      delete m_traits;
      m_traits = nullptr;
    }
  }

  void simplify(const Polygon_2& pgn) {
    Construct_curves_2 ctr_curves =
      reinterpret_cast<const Gt2*>(m_traits)->construct_curves_2_object();

    std::list<Meta_X_monotone_curve_2> curves_list;

    std::pair<Curve_const_iterator, Curve_const_iterator> itr_pair =
      ctr_curves(pgn);

    std::size_t index = 0;
    for (Curve_const_iterator itr = itr_pair.first; itr != itr_pair.second;
         ++itr, ++index) {
      Curve_data cv_data(1, 0, index);
      curves_list.push_back(Meta_X_monotone_curve_2(*itr, cv_data));
    }
    m_traits->set_polygon_size(curves_list.size());

    m_surface_sweep.sweep(curves_list.begin(), curves_list.end());

    // we use the first face with out outer ccbs. This assumpsion should
    // be fixed when we can make a face with no outer ccb to a face with
    // outer ccb.
    Face_iterator it;
    for (it = m_arr->faces_begin(); it != m_arr->faces_end(); ++it)
      if (it->number_of_outer_ccbs() == 0) break;
    CGAL_assertion(it != m_arr->faces_end());

    m_faces_hash[it] = 0;
    Bfs_visitor visitor(&m_edges_hash, &m_faces_hash, 1);
    visitor.visit_ubf(it, 0);
    Bfs_scanner scanner(visitor);
    scanner.scan(*m_arr);
    visitor.after_scan(*m_arr);
  }

  /*! Obtain the arrangement.
   */
  const Arr& arrangement() const { return (*m_arr); }

  /*! Obtain the arrangement.
   */
  Arr& arrangement() { return (*m_arr); }
};

} // namespace CGAL

#endif
