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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

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
  typedef Arrangement_                                  Arrangement_2;

  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Arrangement_2::Topology_traits       Topology_traits;

  typedef Arrangement_2                                 Arr;
  typedef Geometry_traits_2                             Gt2;
  typedef Topology_traits                               Tt;

  typedef typename Gt2::Curve_const_iterator            Curve_const_iterator;
  typedef typename Gt2::Polygon_2                       Polygon_2;
  typedef typename Gt2::Polygon_with_holes_2            Polygon_with_holes_2;
  typedef typename Gt2::Construct_curves_2              Construct_curves_2;

  typedef Gps_simplifier_traits<Gt2>                    Mgt2;
  typedef typename Mgt2::Curve_data                     Curve_data;
  typedef typename Mgt2::X_monotone_curve_2             Meta_X_monotone_curve_2;

  typedef typename Arr::Halfedge_handle                 Halfedge_handle;
  typedef typename Arr::Halfedge_iterator               Halfedge_iterator;
  typedef typename Arr::Face_handle                     Face_handle;
  typedef typename Arr::Face_iterator                   Face_iterator;
  typedef typename Arr::Edge_iterator                   Edge_iterator;
  typedef typename Arr::Vertex_handle                   Vertex_handle;
  typedef typename Arr::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Arr::Ccb_halfedge_circulator         Ccb_halfedge_circulator;
  typedef typename Arr::Allocator                       Allocator;

  // We obtain a proper helper type from the topology traits of the arrangement.
  // However, the arrangement is parametrized with the Gt2 geometry traits,
  // while we need the Mgt2 geometry traits (which derives from Gt2).
  // Thus, we rebind the helper.
  // We cannot parameterized the arrangement with the Mgt2 geometry
  // traits to start with, because it extends the curve type with arrangement
  // dependent types. (It is parameterized with the arrangement type.)
  typedef Indexed_event<Mgt2, Arr, Allocator>           Event;
  typedef Arr_construction_subcurve<Mgt2, Event, Allocator>
                                                        Subcurve;
  typedef typename Tt::template Construction_helper<Event, Subcurve>
                                                        Helper_tmp;
  typedef typename Helper_tmp::template rebind<Mgt2, Arr, Event, Subcurve>::other
                                                        Helper;
  typedef Gps_agg_op_base_visitor<Helper, Arr>          Visitor;
  typedef Ss2::Surface_sweep_2<Visitor>                 Surface_sweep_2;

  typedef Unique_hash_map<Halfedge_handle, unsigned int>
                                                        Edges_hash;

  typedef Unique_hash_map<Face_handle, unsigned int>    Faces_hash;
  typedef Gps_bfs_join_visitor<Arr>                     Bfs_visitor;
  typedef Gps_bfs_scanner<Arr, Bfs_visitor>             Bfs_scanner;

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
  ~Gps_polygon_simplifier()
  {
    if (m_own_traits && (m_traits != nullptr)) {
      delete m_traits;
      m_traits = nullptr;
    }
  }

  void simplify(const Polygon_2& pgn)
  {
    Construct_curves_2 ctr_curves =
      reinterpret_cast<const Gt2*>(m_traits)->construct_curves_2_object();

    std::list<Meta_X_monotone_curve_2> curves_list;

    std::pair<Curve_const_iterator, Curve_const_iterator> itr_pair =
      ctr_curves(pgn);

    unsigned int index = 0;
    for (Curve_const_iterator itr = itr_pair.first; itr != itr_pair.second;
         ++itr, ++index)
    {
      Curve_data cv_data(1, 0, index);
      curves_list.push_back(Meta_X_monotone_curve_2(*itr, cv_data));
    }
    m_traits->set_polygon_size(static_cast<unsigned int>(curves_list.size()));

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
