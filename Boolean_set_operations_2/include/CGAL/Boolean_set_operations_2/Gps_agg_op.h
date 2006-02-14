// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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

#ifndef GSP_AGG_OP_H
#define GSP_AGG_OP_H

#include <CGAL/Boolean_set_operations_2/Gps_agg_meta_traits.h>
#include <CGAL/Sweep_line_2.h>
#include <CGAL/Sweep_line_2/Arr_construction_curve.h>
#include <CGAL/Sweep_line_2/Arr_construction_event.h>

#include <CGAL/Boolean_set_operations_2/Gps_agg_op_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_scanner.h>
//#include <CGAL/Boolean_set_operations_2/Gps_insertion_meta_traits.h>
#include <CGAL/Unique_hash_map.h> 
#include <CGAL/Arr_accessor.h>
#include <CGAL/iterator.h> 

CGAL_BEGIN_NAMESPACE

template <class Arrangement_, class Bfs_visitor_>
class Gps_agg_op
{
  typedef Arrangement_                                Arrangement_2;
  typedef typename Arrangement_2::Traits_2            Traits_2;
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
  typedef Arr_construction_curve<Meta_traits>         Subcurve; 
  typedef Arr_construction_event<Meta_traits,
                                 Subcurve,
                                 Halfedge_handle>     Event;

  typedef Gps_agg_op_visitor<Meta_traits,
                             Arrangement_2,
                             Event,
                             Subcurve>                Visitor;

  typedef Sweep_line_2<Meta_traits,
		                   Visitor,
                       Subcurve,
                       Event>                         Sweep_line_2;
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
  Gps_agg_op (Arrangement_2& arr, Traits_2& tr) :
    m_arr (&arr),
    m_traits(new Meta_traits(tr)),
    m_visitor (&arr, &m_edges_hash),
    m_sweep_line (m_traits, &m_visitor)
  {}

 

  void sweep_arrangements(unsigned int lower,
                          unsigned int upper,
                          unsigned int jump,
                          std::vector<Arrangement_2*>&  arr_vec)
  {
    std::list<Meta_X_monotone_curve_2> curves_list;

    typename Meta_traits::Compare_endpoints_xy_2 cmp_endpoints = 
      m_traits->compare_endpoints_xy_2_object();

    unsigned int n_inf_pgn = 0; // number of infinte polygons (arrangement 
                                // with a contained unbounded face
    unsigned int n_pgn = 0; // number of polygons (arrangements)
    for(unsigned int i=lower; i<=upper; i+=jump, ++n_pgn)
    {
      Arrangement_2* arr = arr_vec[i];
      if(arr->unbounded_face()->contained())
        ++n_inf_pgn;

      Edge_iterator itr = arr->edges_begin();
      for(; itr != arr->edges_end(); ++itr)
      {
        // take only relevant edges (which seperate between contained and 
        // non-contained faces.
        Halfedge_iterator he = itr;
        if(he->face()->contained() == he->twin()->face()->contained())
          continue;

        if(he->direction() == LARGER)
          he = he->twin();

        Curve_data cv_data(arr, he, 1, 0);
        curves_list.push_back(Meta_X_monotone_curve_2(he->curve(), cv_data));
      }
    }
    m_sweep_line.sweep(curves_list.begin(), curves_list.end());

    m_faces_hash[m_arr->unbounded_face()] = n_inf_pgn; 
    Bfs_visitor visitor(&m_edges_hash, &m_faces_hash, n_pgn);
    visitor.visit_ubf(m_arr->unbounded_face(), n_inf_pgn);
    Bfs_scanner scanner(visitor);
    scanner.scan(*m_arr);
    visitor.after_scan(*m_arr);
  }

  //Arrangement_2* create_clean_arr()
  //{
  //  typedef Gps_insertion_meta_traits<Arrangement_2>  Insert_meta_traits;
  //  typedef typename Insert_meta_traits::Curve_data            Curve_data;
  //  typedef Arr_construction_curve<Insert_meta_traits>         Subcurve; 
  //  typedef Arr_construction_event<Insert_meta_traits,
  //                                 Subcurve,
  //                                 Halfedge_handle>            Event;
  //  typedef typename Insert_meta_traits::X_monotone_curve_2    Ex_X_monotone_curve_2;
  //  typedef Arr_construction_visitor<Insert_meta_traits,
  //                                   Arrangement_2,
  //                                   Event,
  //                                   Subcurve>                 Insertion_visitor;
  //  typedef Basic_sweep_line_2<Insert_meta_traits,
		//                           Insertion_visitor,
  //                             Subcurve,
  //                             Event>                          Basic_sweep_line_2;

  // 
  //                                   
  //  Arrangement_2* new_arr = new Arrangement_2();
  //  std::vector<Ex_X_monotone_curve_2> xcurves_vec;
  //  for(Edge_iterator itr = m_arr->edges_begin();
  //      itr != m_arr->edges_end();
  //      ++itr)
  //  {
  //    Halfedge_handle he = itr;
  //    if(he->face()->contained() == he->twin()->face()->contained())
  //      continue;  //redundent edge, continue.

  //    if(he->direction() == LARGER)
  //      he = he->twin();

  //    xcurves_vec.push_back(Ex_X_monotone_curve_2(he->curve(), Curve_data(he)));
  //  }

  //  std::cout<<"number of edges in clear arr: " << xcurves_vec.size()<<"\n";
  //  Insertion_visitor visitor(new_arr);
  //  Basic_sweep_line_2  sl (&visitor);
  //  sl.sweep(xcurves_vec.begin(), xcurves_vec.end());
  //}
       
  ~Gps_agg_op()
  {
    delete m_traits;
  }
};

CGAL_END_NAMESPACE

#endif
