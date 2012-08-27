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
// $Id$ $Date$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ophir Setter    <ophir.setter@cs.tau.ac.il>

#ifndef CGAL_GPS_AGG_OP_H
#define CGAL_GPS_AGG_OP_H

/*!
  \file   Gps_agg_op.h
  \brief  The class Gps_agg_op is responsible for aggregated Boolean set 
          operations depending on a visitor template parameter.
          It uses the sweep-line algorithm from the arrangement packages
          to overlay all the polygon sets, and then it uses a BFS that 
          determines which of the faces is contained in the result using
          the visitor.
*/


#include <CGAL/Boolean_set_operations_2/Gps_agg_meta_traits.h>
#include <CGAL/Boolean_set_operations_2/Gps_agg_op_sweep.h>
#include <CGAL/Sweep_line_2/Arr_construction_subcurve.h>
#include <CGAL/Sweep_line_2/Arr_construction_event.h>

#include <CGAL/Boolean_set_operations_2/Gps_agg_op_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_scanner.h>
//#include <CGAL/Boolean_set_operations_2/Gps_insertion_meta_traits.h>
#include <CGAL/Unique_hash_map.h> 
#include <CGAL/Arr_accessor.h>
#include <CGAL/iterator.h> 

namespace CGAL {

template <class Arrangement_, class Bfs_visitor_>
class Gps_agg_op
{
  typedef Arrangement_                                Arrangement_2;
  typedef typename Arrangement_2::Traits_adaptor_2    Traits_2;
  typedef typename Traits_2::Curve_const_iterator     Curve_const_iterator;
  typedef Gps_agg_meta_traits<Arrangement_2>          Meta_traits;
  typedef typename Meta_traits::Curve_data            Curve_data;
  typedef typename Meta_traits::X_monotone_curve_2    Meta_X_monotone_curve_2;
  typedef typename Arrangement_2::Halfedge_handle     Halfedge_handle;
  typedef typename Arrangement_2::Halfedge_iterator   Halfedge_iterator;
  typedef typename Arrangement_2::Face_handle         Face_handle;
  typedef typename Arrangement_2::Edge_iterator       Edge_iterator;
  typedef typename Arrangement_2::Vertex_handle       Vertex_handle;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator 
                                                      Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Ccb_halfedge_circulator 
                                                      Ccb_halfedge_circulator;

  typedef std::pair<Arrangement_2 *,
                    std::vector<Vertex_handle> *>     Arr_entry;

  typedef Arr_construction_subcurve<Meta_traits>      Subcurve; 

  typedef Arr_construction_event<Meta_traits,
                                 Subcurve,
//				 Halfedge_handle>       Base_event;
                                 Arrangement_2>       Base_event;

  typedef Indexed_event<Base_event>                   Event;

  typedef Gps_agg_op_visitor<Meta_traits,
                             Arrangement_2,
                             Event,
                             Subcurve>                Visitor;

  typedef Gps_agg_op_sweep_line_2<Arrangement_2,
                                  Meta_traits,
                                  Visitor,
                                  Subcurve,
                                  Event>              Sweep_line_2;

  typedef Unique_hash_map<Halfedge_handle, 
                          unsigned int>               Edges_hash;

  typedef Unique_hash_map<Face_handle, 
                          unsigned int>               Faces_hash;
  typedef Bfs_visitor_                                Bfs_visitor;
  typedef Gps_bfs_scanner<Arrangement_2, Bfs_visitor> Bfs_scanner;
 
protected:
  Arrangement_2*       m_arr;
  Meta_traits*         m_traits;
  Visitor              m_visitor;
  Sweep_line_2         m_sweep_line;
  Edges_hash           m_edges_hash; // maps halfedge to its BC (boundary counter)
  Faces_hash           m_faces_hash;  // maps face to its IC (inside count)
  
public:

  /*! Constructor. */
  Gps_agg_op (Arrangement_2& arr, std::vector<Vertex_handle>& vert_vec,
              const Traits_2 & tr) :
    m_arr (&arr),
    m_traits(new Meta_traits(tr)),
    m_visitor (&arr, &m_edges_hash, &vert_vec),
    m_sweep_line (m_traits, &m_visitor)
  {}

  void sweep_arrangements(unsigned int lower,
                          unsigned int upper,
                          unsigned int jump,
                          std::vector<Arr_entry>& arr_vec)
  {
    std::list<Meta_X_monotone_curve_2> curves_list;

    unsigned int n_inf_pgn = 0; // number of infinte polygons (arrangement 
                                // with a contained unbounded face
    unsigned int n_pgn = 0;     // number of polygons (arrangements)
    unsigned int i;

    for (i = lower; i <= upper; i += jump, ++n_pgn)
    {
      // The BFS scan (after the loop) starts in the reference face,
      // so we count the number of polygons that contain the reference face.
      Arrangement_2* arr = (arr_vec[i]).first;
      if (arr->reference_face()->contained())
        ++n_inf_pgn;

      Edge_iterator  itr = arr->edges_begin();
      for(; itr != arr->edges_end(); ++itr)
      {
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

    m_sweep_line.sweep (curves_list.begin(), curves_list.end(),
                        lower, upper, jump,
                        arr_vec);

    m_faces_hash[m_arr->reference_face()] = n_inf_pgn; 
    Bfs_visitor visitor(&m_edges_hash, &m_faces_hash, n_pgn);
    visitor.visit_ubf(m_arr->faces_begin(), n_inf_pgn);
    Bfs_scanner scanner(visitor);
    scanner.scan(*m_arr);
    visitor.after_scan(*m_arr);
  }

  ~Gps_agg_op()
  {
    delete m_traits;
  }
};

} //namespace CGAL

#endif
