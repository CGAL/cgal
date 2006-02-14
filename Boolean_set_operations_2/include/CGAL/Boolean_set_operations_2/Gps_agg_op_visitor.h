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

#ifndef GSP_AGG_OP_VISITOR
#define GSP_AGG_OP_VISITOR

#include <CGAL/Unique_hash_map.h> 
#include <CGAL/Sweep_line_2/Arr_construction_visitor.h>

CGAL_BEGIN_NAMESPACE

template<class Traits, class Arrangement_, class Event,class Subcurve>
class Gps_agg_op_visitor : 
  public Arr_construction_visitor<Traits,Arrangement_,Event,Subcurve>
{
  protected:

  typedef Arrangement_                                     Arrangement;
  
  typedef Arr_construction_visitor<Traits,
                                   Arrangement,
                                   Event,
                                   Subcurve>               Base;

  typedef Gps_agg_op_visitor<Traits,
                             Arrangement,
                             Event,
                             Subcurve>                     Self;

  typedef typename Base::SL_iterator                       SL_iterator;
  typedef typename Base::Halfedge_handle                   Halfedge_handle;
  typedef typename Base::Vertex_handle                     Vertex_handle;
  typedef typename Base::SubCurveIter                      SubCurveIter;
  typedef typename Base::SubCurveRevIter                   SubCurveRevIter;
  typedef typename Traits::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Traits::Point_2                         Point_2;
  
  typedef typename Arrangement::Face_handle                Face_handle;
  typedef typename Arrangement::Face_const_handle          Face_const_handle;
  typedef Unique_hash_map<Halfedge_handle, unsigned int>   Edges_hash;

protected:

  Edges_hash*  m_edges_hash; // maps halfedges to their BC (coundary counter)


public:

  Gps_agg_op_visitor(Arrangement *arr,
                     Edges_hash* hash): Base(arr),
                                        m_edges_hash(hash)
  {}

  virtual Halfedge_handle insert_in_face_interior(const X_monotone_curve_2& cv,
                                          Subcurve* sc)
  {
    Halfedge_handle he = 
      Base::insert_in_face_interior(cv, sc);
    insert_edge_to_hash(he, cv);
    return (he);
  }

  virtual Halfedge_handle insert_at_vertices(const X_monotone_curve_2& cv,
                                             Halfedge_handle hhandle,
                                             Halfedge_handle prev,
                                             Subcurve* sc,
                                             bool &new_face_created)
  {
    Halfedge_handle res_he =
      Base::insert_at_vertices(cv, hhandle, prev, sc, new_face_created);
    insert_edge_to_hash(res_he, cv);
    return (res_he);
  }

  virtual Halfedge_handle insert_from_right_vertex
                          (const X_monotone_curve_2& cv,
                           Halfedge_handle he,
                           Subcurve* sc)
  {
    Halfedge_handle res_he = 
      Base::insert_from_right_vertex(cv, he, sc);
    insert_edge_to_hash(res_he, cv);
    return (res_he);
  }

  virtual Halfedge_handle insert_from_left_vertex
                          (const X_monotone_curve_2& cv,
                           Halfedge_handle he,
                           Subcurve* sc)
  {
    Halfedge_handle res_he = 
      Base::insert_from_left_vertex(cv, he, sc);
    insert_edge_to_hash(res_he, cv);
    return (res_he);
  }



private:

  void insert_edge_to_hash(Halfedge_handle he, const X_monotone_curve_2& cv)
  {
    Comparison_result he_dir = he->direction();
    Comparison_result cv_dir =
      this->m_arr->get_traits()->compare_endpoints_xy_2_object()(cv);

    if(he_dir == cv_dir)
    {
      (*m_edges_hash)[he] = cv.data().bc();
      (*m_edges_hash)[he->twin()] = cv.data().twin_bc();
    }
    else
    {
      (*m_edges_hash)[he] = cv.data().twin_bc();
      (*m_edges_hash)[he->twin()] = cv.data().bc();
    }
  }
 

};

CGAL_END_NAMESPACE

#endif
