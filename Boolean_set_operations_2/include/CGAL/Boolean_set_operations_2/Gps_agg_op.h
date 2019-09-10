// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ophir Setter    <ophir.setter@cs.tau.ac.il>

#ifndef CGAL_BSO_2_GPS_AGG_OP_H
#define CGAL_BSO_2_GPS_AGG_OP_H

#include <CGAL/license/Boolean_set_operations_2.h>

/*! \file Gps_agg_op.h
 *
 * The class Gps_agg_op is responsible for aggregated Boolean set operations
 * depending on a visitor template parameter.  It uses the surface-sweep
 * algorithm from the arrangement packages to overlay all the polygon sets, and
 * then it uses a BFS that determines which of the faces is contained in the
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

template <typename Arrangement_, typename BfsVisitor>
class Gps_agg_op {
  typedef Arrangement_                                  Arrangement_2;
  typedef BfsVisitor                                    Bfs_visitor;

  typedef typename Arrangement_2::Traits_adaptor_2      Geometry_traits_2;
  typedef typename Arrangement_2::Topology_traits       Topology_traits;

  typedef Arrangement_2                                 Arr;
  typedef Geometry_traits_2                             Gt2;
  typedef Topology_traits                               Tt;

  typedef typename Gt2::Curve_const_iterator            Curve_const_iterator;
  typedef Gps_agg_meta_traits<Arr>                      Mgt2;
  typedef typename Mgt2::Curve_data                     Curve_data;
  typedef typename Mgt2::X_monotone_curve_2             Meta_X_monotone_curve_2;

  typedef typename Arr::Halfedge_handle                 Halfedge_handle;
  typedef typename Arr::Halfedge_iterator               Halfedge_iterator;
  typedef typename Arr::Face_handle                     Face_handle;
  typedef typename Arr::Edge_iterator                   Edge_iterator;
  typedef typename Arr::Vertex_handle                   Vertex_handle;
  typedef typename Arr::Allocator                       Allocator;

  typedef std::pair<Arr*, std::vector<Vertex_handle> *> Arr_entry;

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
  typedef Gps_agg_op_visitor<Helper, Arr>               Visitor;
  typedef Gps_agg_op_surface_sweep_2<Arr, Visitor>      Surface_sweep_2;

  typedef Unique_hash_map<Halfedge_handle, unsigned int>
                                                        Edges_hash;

  typedef Unique_hash_map<Face_handle, unsigned int>    Faces_hash;
  typedef Gps_bfs_scanner<Arr, Bfs_visitor>             Bfs_scanner;

protected:
  Arr* m_arr;
  Mgt2* m_traits;
  Visitor m_visitor;
  Surface_sweep_2 m_surface_sweep;
  Edges_hash m_edges_hash;      // maps halfedge to its BC (boundary counter)
  Faces_hash m_faces_hash;      // maps face to its IC (inside count)

public:
  /*! Constructor. */
  Gps_agg_op(Arr& arr, std::vector<Vertex_handle>& vert_vec, const Gt2& tr) :
    m_arr(&arr),
    m_traits(new Mgt2(tr)),
    m_visitor(&arr, &m_edges_hash, &vert_vec),
    m_surface_sweep(m_traits, &m_visitor)
  {}

  void sweep_arrangements(unsigned int lower, unsigned int upper,
                          unsigned int jump, std::vector<Arr_entry>& arr_vec)
  {
    std::list<Meta_X_monotone_curve_2> curves_list;

    unsigned int n_inf_pgn = 0; // number of infinte polygons (arrangement
                                // with a contained unbounded face
    unsigned int n_pgn = 0;     // number of polygons (arrangements)
    unsigned int i;

    for (i = lower; i <= upper; i += jump, ++n_pgn) {
      // The BFS scan (after the loop) starts in the reference face,
      // so we count the number of polygons that contain the reference face.
      Arr* arr = (arr_vec[i]).first;
      if (arr->reference_face()->contained()) ++n_inf_pgn;

      Edge_iterator  itr = arr->edges_begin();
      for(; itr != arr->edges_end(); ++itr) {
        // take only relevant edges (which seperate between contained and
        // non-contained faces.
        Halfedge_iterator he = itr;
        if(he->face()->contained() == he->twin()->face()->contained())
          continue;
        if ((Arr_halfedge_direction)he->direction() == ARR_RIGHT_TO_LEFT)
          he = he->twin();

        Curve_data cv_data(arr, he, 1, 0);
        curves_list.push_back(Meta_X_monotone_curve_2(he->curve(), cv_data));
      }
    }

    m_surface_sweep.sweep(curves_list.begin(), curves_list.end(),
                          lower, upper, jump, arr_vec);

    m_faces_hash[m_arr->reference_face()] = n_inf_pgn;
    Bfs_visitor visitor(&m_edges_hash, &m_faces_hash, n_pgn);
    visitor.visit_ubf(m_arr->faces_begin(), n_inf_pgn);
    Bfs_scanner scanner(visitor);
    scanner.scan(*m_arr);
    visitor.after_scan(*m_arr);
  }

  /*! Destruct.
   */
  ~Gps_agg_op() { delete m_traits; }
};

} //namespace CGAL

#endif
