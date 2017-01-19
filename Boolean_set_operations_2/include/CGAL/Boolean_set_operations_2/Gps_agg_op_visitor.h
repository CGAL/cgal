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
//                 Ron Wein        <wein@post.tau.ac.il>

#ifndef CGAL_GSP_AGG_OP_VISITOR_H
#define CGAL_GSP_AGG_OP_VISITOR_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <CGAL/Unique_hash_map.h> 
#include <CGAL/Sweep_line_2/Arr_construction_sl_visitor.h>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_topology_traits/Arr_bounded_planar_construction_helper.h>
#include <CGAL/Arr_topology_traits/Arr_unb_planar_construction_helper.h>

namespace CGAL {

template<class Traits, class Arrangement_, class Event,class Subcurve>
class Gps_agg_op_base_visitor :
  public
  Arr_construction_sl_visitor<
    // TODO derive (helper) from topology traits class
    typename boost::mpl::if_< 
    boost::is_same< typename Arr_are_all_sides_oblivious_tag< 
                                     typename Traits::Left_side_category, 
                                     typename Traits::Bottom_side_category, 
                                     typename Traits::Top_side_category, 
                                     typename Traits::Right_side_category 
    >::result, Arr_all_sides_oblivious_tag >,
    Arr_bounded_planar_construction_helper<Traits, 
                                       Arrangement_,
                                       Event,
                                       Subcurve>,
    Arr_unb_planar_construction_helper<Traits,
                                       Arrangement_,
                                       Event,
                                       Subcurve> 
    >::type
  >
{
  protected:
  typedef Arrangement_                                     Arrangement;


  
  typedef typename boost::mpl::if_< 
    boost::is_same< typename Arr_are_all_sides_oblivious_tag< 
                                     typename Traits::Left_side_category, 
                                     typename Traits::Bottom_side_category, 
                                     typename Traits::Top_side_category, 
                                     typename Traits::Right_side_category 
    >::result, Arr_all_sides_oblivious_tag >,
    Arr_bounded_planar_construction_helper<Traits, 
                                       Arrangement,
                                       Event,
                                       Subcurve>,
    Arr_unb_planar_construction_helper<Traits,
                                       Arrangement,
                                       Event,
                                       Subcurve> 
    >::type Construction_helper;
  typedef Arr_construction_sl_visitor<Construction_helper> Base;

  typedef typename Base::Status_line_iterator              SL_iterator;
  typedef typename Base::Halfedge_handle                   Halfedge_handle;
  typedef typename Base::Vertex_handle                     Vertex_handle;
  typedef typename Base::Event_subcurve_iterator           SubCurveIter;
  typedef typename Base::Event_subcurve_reverse_iterator   SubCurveRevIter;
  typedef typename Traits::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Traits::Point_2                         Point_2;
  
  typedef Unique_hash_map<Halfedge_handle, unsigned int>   Edges_hash;

protected:

  Edges_hash*  m_edges_hash; // maps halfedges to their BC (coundary counter)


public:

  Gps_agg_op_base_visitor(Arrangement *arr,
                          Edges_hash* hash): Base(arr),
                                             m_edges_hash(hash)
  {}

  // TODO (IMPORTANT): unbounded helper might be not fully supported
  // TODO add mpl-warning

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
    const Comparison_result he_dir = 
      ((Arr_halfedge_direction)he->direction() == ARR_LEFT_TO_RIGHT) ? SMALLER : LARGER;

    const Comparison_result cv_dir =
      this->m_arr_access.arrangement().geometry_traits()->
            compare_endpoints_xy_2_object()(cv);

    if (he_dir == cv_dir)
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

template <class BaseEvent_>
class Indexed_event : public BaseEvent_
{
private:

  unsigned int    m_index;

public:

  Indexed_event () :
    BaseEvent_(),
    m_index (0)
  {}

  unsigned int index () const
  {
    return (m_index);
  }

  void set_index (unsigned int index)
  {
    m_index = index;
    return;
  }
};

template<class Traits, class Arrangement_, class Event, class Subcurve>
class Gps_agg_op_visitor : 
  public Gps_agg_op_base_visitor<Traits, Arrangement_, Event, Subcurve>
{
protected:

  typedef Arrangement_                                    Arrangement;

  typedef Gps_agg_op_base_visitor<Traits,
                                  Arrangement,
                                  Event,
                                  Subcurve>               Base;

  typedef typename Base::SL_iterator                       SL_iterator;
  typedef typename Base::Halfedge_handle                   Halfedge_handle;
  typedef typename Base::Vertex_handle                     Vertex_handle;
  typedef typename Base::SubCurveIter                      SubCurveIter;
  typedef typename Base::SubCurveRevIter                   SubCurveRevIter;
  typedef typename Traits::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Traits::Point_2                         Point_2;
  
protected:

  unsigned int m_event_count; // The number of events so far.
  std::vector<Vertex_handle>   *m_vertices_vec;  // The vertices, sorted in
                                                 // ascending order.

public:

  Gps_agg_op_visitor (Arrangement *arr,
                      typename Base::Edges_hash* hash,
                      std::vector<Vertex_handle>* vertices_vec) :
    Base (arr, hash),
    m_event_count (0),
    m_vertices_vec (vertices_vec)
  {}

  void before_handle_event (Event* event)
  {
    event->set_index (m_event_count);
    m_event_count++;

    return;
  }

  virtual Halfedge_handle 
  insert_in_face_interior (const X_monotone_curve_2& cv,
                           Subcurve* sc)
  {
    Halfedge_handle res_he = Base::insert_in_face_interior(cv, sc);

    // We now have a halfedge whose source vertex is associated with the
    // last event and whose target vertex is associated with the current event:
    Event *curr_event = reinterpret_cast<Event*>(this->current_event());
    Event *last_event = reinterpret_cast<Event*>((sc)->last_event());

    CGAL_assertion ((Arr_halfedge_direction)res_he->direction() == ARR_LEFT_TO_RIGHT);
    _insert_vertex (curr_event, res_he->target());
    _insert_vertex (last_event, res_he->source());

    return (res_he);
  }

  virtual Halfedge_handle
  insert_from_right_vertex (const X_monotone_curve_2& cv,
                            Halfedge_handle he,
                            Subcurve* sc)
  {
    Halfedge_handle  res_he = Base::insert_from_right_vertex (cv, he, sc);

    // We now have a halfedge whose target vertex is associated with the
    // last event (we have already dealt with its source vertex).
    Event *last_event = reinterpret_cast<Event*>((sc)->last_event());
    CGAL_assertion ((Arr_halfedge_direction)res_he->direction() == ARR_RIGHT_TO_LEFT);
    _insert_vertex (last_event, res_he->target());

    return (res_he);
  }

  virtual Halfedge_handle
  insert_from_left_vertex (const X_monotone_curve_2& cv,
                           Halfedge_handle he,
                           Subcurve* sc)
  {
    Halfedge_handle  res_he = Base::insert_from_left_vertex (cv, he, sc);

    // We now have a halfedge whose target vertex is associated with the
    // current event (we have already dealt with its source vertex).
    Event *curr_event = reinterpret_cast<Event*>(this->current_event());

    CGAL_assertion ((Arr_halfedge_direction)res_he->direction() == ARR_LEFT_TO_RIGHT);
    _insert_vertex (curr_event, res_he->target());

    return (res_he);
  }

private:

  void _insert_vertex (const Event* event,
                       Vertex_handle v)
  {
    const unsigned int    index = event->index();
    
    if (index >= m_vertices_vec->size())
      m_vertices_vec->resize (2 * (index + 1));

    (*m_vertices_vec)[index] = v;
    return;
  }

};

} //namespace CGAL

#endif
