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
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_GPS_POLYGON_SIMPILFIER_H
#define CGAL_GPS_POLYGON_SIMPILFIER_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <CGAL/Boolean_set_operations_2/Gps_simplifier_traits.h>
#include <CGAL/Sweep_line_2.h>
#include <CGAL/Sweep_line_2/Arr_construction_subcurve.h>
#include <CGAL/Sweep_line_2/Arr_construction_event.h>

#include <CGAL/Boolean_set_operations_2/Gps_agg_op_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_scanner.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_join_visitor.h>

#include <CGAL/Unique_hash_map.h>
#include <CGAL/Arr_accessor.h>
#include <CGAL/iterator.h>

namespace CGAL {

template <class Arrangement_>
class Gps_polygon_simplifier
{
  typedef Arrangement_                                Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2   Traits_2;
  typedef typename Traits_2::Curve_const_iterator     Curve_const_iterator;
  typedef typename Traits_2::Polygon_2                Polygon_2;
  typedef typename Traits_2::Polygon_with_holes_2     Polygon_with_holes_2;
  typedef typename Traits_2::Construct_curves_2       Construct_curves_2;

  typedef Gps_simplifier_traits<Traits_2>             Meta_traits;
  typedef typename Meta_traits::Curve_data            Curve_data;
  typedef typename Meta_traits::X_monotone_curve_2    Meta_X_monotone_curve_2;
  typedef typename Arrangement_2::Halfedge_handle     Halfedge_handle;
  typedef typename Arrangement_2::Halfedge_iterator   Halfedge_iterator;
  typedef typename Arrangement_2::Face_handle         Face_handle;
  typedef typename Arrangement_2::Face_iterator       Face_iterator;
  typedef typename Arrangement_2::Edge_iterator       Edge_iterator;
  typedef typename Arrangement_2::Vertex_handle       Vertex_handle;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Ccb_halfedge_circulator
                                                      Ccb_halfedge_circulator;
  typedef Arr_construction_subcurve<Meta_traits>      Subcurve;
  typedef Arr_construction_event<Meta_traits,
                                 Subcurve,
                                 Arrangement_2>       Event;

  typedef Gps_agg_op_base_visitor<Meta_traits,
                                  Arrangement_2,
                                  Event,
                                  Subcurve>           Visitor;

  typedef CGAL::Sweep_line_2<Meta_traits,
                             Visitor,
                             Subcurve,
                             Event>                   Sweep_line_2;

  typedef Unique_hash_map<Halfedge_handle,
                          unsigned int>               Edges_hash;

  typedef Unique_hash_map<Face_handle,
                          unsigned int>               Faces_hash;
  typedef Gps_bfs_join_visitor<Arrangement_2>         Bfs_visitor;
  typedef Gps_bfs_scanner<Arrangement_2, Bfs_visitor> Bfs_scanner;

protected:
  Arrangement_2* m_arr;
  const Meta_traits* m_traits;
  bool                 m_own_traits;
  Visitor              m_visitor;
  Sweep_line_2         m_sweep_line;
  Edges_hash           m_edges_hash; // maps halfedge to its BC (boundary counter)
  Faces_hash           m_faces_hash;  // maps face to its IC (inside count)

public:
   /*! Constructor. */
  Gps_polygon_simplifier(Arrangement_2& arr, const Traits_2& tr) :
    m_arr(&arr),
    m_traits(new Meta_traits(tr)),
    m_own_traits(true),
    m_visitor(&arr, &m_edges_hash),
    m_sweep_line(m_traits, &m_visitor)
  {}

  /*! Destructor. */
  ~Gps_polygon_simplifier()
  {
    if (m_own_traits && (m_traits != NULL)) {
      delete m_traits;
      m_traits = NULL;
    }
  }

  void simplify(const Polygon_2& pgn)
  {
    Construct_curves_2 ctr_curves =
      reinterpret_cast<const Traits_2*>(m_traits)->construct_curves_2_object();

    std::list<Meta_X_monotone_curve_2> curves_list;

    std::pair<Curve_const_iterator,
              Curve_const_iterator>  itr_pair = ctr_curves(pgn);

    unsigned int index = 0;
    for(Curve_const_iterator itr = itr_pair.first;
        itr != itr_pair.second;
        ++itr, ++index)
    {
      Curve_data cv_data(1, 0, index);
      curves_list.push_back(Meta_X_monotone_curve_2(*itr, cv_data));
    }
    m_traits->set_polygon_size(static_cast<unsigned int>(curves_list.size()));

    m_sweep_line.sweep(curves_list.begin(), curves_list.end());

    // we use the first face with out outer ccbs. This assumpsion should
    // be fixed when we can make a face with no outer ccb to a face with
    // outer ccb.
    Face_iterator it;
    for (it = m_arr->faces_begin(); it != m_arr->faces_end(); ++it)
      if (it->number_of_outer_ccbs() == 0)
        break;
    CGAL_assertion(it != m_arr->faces_end());

    m_faces_hash[it] = 0;
    Bfs_visitor visitor(&m_edges_hash, &m_faces_hash, 1);
    visitor.visit_ubf(it, 0);
    Bfs_scanner scanner(visitor);
    scanner.scan(*m_arr);
    visitor.after_scan(*m_arr);
  }

  const Arrangement_2& arrangement() const
  {
    return (*m_arr);
  }

  Arrangement_2& arrangement()
  {
    return (*m_arr);
  }

};
} //namespace CGAL

#endif
