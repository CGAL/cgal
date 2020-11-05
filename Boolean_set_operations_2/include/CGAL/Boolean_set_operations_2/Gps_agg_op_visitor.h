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
//                 Ron Wein        <wein@post.tau.ac.il>

#ifndef CGAL_BSO_2_GSP_AGG_OP_VISITOR_H
#define CGAL_BSO_2_GSP_AGG_OP_VISITOR_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/Unique_hash_map.h>
#include <CGAL/Surface_sweep_2/Arr_construction_ss_visitor.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_topology_traits/Arr_bounded_planar_construction_helper.h>
#include <CGAL/Arr_topology_traits/Arr_unb_planar_construction_helper.h>
#include <CGAL/Default.h>

namespace CGAL {

template <typename Helper_, typename Arrangement_, typename Visitor_ = Default>
class Gps_agg_op_base_visitor :
  public Arr_construction_ss_visitor<
    Helper_,
    typename Default::Get<Visitor_, Gps_agg_op_base_visitor<Helper_,
                                                            Arrangement_,
                                                            Visitor_> >::type>
{
public:
  typedef Helper_                                       Helper;
  typedef Arrangement_                                  Arrangement_2;

  typedef typename Helper::Geometry_traits_2            Geometry_traits_2;
  typedef typename Helper::Event                        Event;
  typedef typename Helper::Subcurve                     Subcurve;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef Arrangement_2                                 Arr;

  typedef Gps_agg_op_base_visitor<Helper, Arr, Visitor_>
                                                        Self;
  typedef typename Default::Get<Visitor_, Self>::type   Visitor;
  typedef Arr_construction_ss_visitor<Helper, Visitor>  Base;

public:
  typedef typename Arr::Halfedge_handle                 Halfedge_handle;
  typedef typename Arr::Vertex_handle                   Vertex_handle;
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Gt2::Point_2                         Point_2;

  typedef Unique_hash_map<Halfedge_handle, unsigned int>
                                                        Edges_hash;

protected:
  Edges_hash* m_edges_hash; // maps halfedges to their BC (coundary counter)

public:
  Gps_agg_op_base_visitor(Arr* arr, Edges_hash* hash) :
    Base(arr),
    m_edges_hash(hash)
  {}

  // TODO (IMPORTANT): unbounded helper might be not fully supported
  // TODO add mpl-warning

  virtual Halfedge_handle insert_in_face_interior(const X_monotone_curve_2& cv,
                                                  Subcurve* sc)
  {
    Halfedge_handle he = Base::insert_in_face_interior(cv, sc);
    insert_edge_to_hash(he, cv);
    return he;
  }

  virtual Halfedge_handle insert_at_vertices(const X_monotone_curve_2& cv,
                                             Halfedge_handle hhandle,
                                             Halfedge_handle prev,
                                             Subcurve* sc,
                                             bool& new_face_created)
  {
    Halfedge_handle res_he =
      Base::insert_at_vertices(cv, hhandle, prev, sc, new_face_created);
    insert_edge_to_hash(res_he, cv);
    return res_he;
  }

  virtual Halfedge_handle insert_from_right_vertex(const X_monotone_curve_2& cv,
                                                   Halfedge_handle he,
                                                   Subcurve* sc)
  {
    Halfedge_handle res_he = Base::insert_from_right_vertex(cv, he, sc);
    insert_edge_to_hash(res_he, cv);
    return res_he;
  }

  virtual Halfedge_handle insert_from_left_vertex(const X_monotone_curve_2& cv,
                                                  Halfedge_handle he,
                                                  Subcurve* sc)
  {
    Halfedge_handle res_he = Base::insert_from_left_vertex(cv, he, sc);
    insert_edge_to_hash(res_he, cv);
    return res_he;
  }

private:
  void insert_edge_to_hash(Halfedge_handle he, const X_monotone_curve_2& cv)
  {
    const Comparison_result he_dir =
      ((Arr_halfedge_direction)he->direction() == ARR_LEFT_TO_RIGHT) ?
      SMALLER : LARGER;

    const Comparison_result cv_dir =
      this->m_arr_access.arrangement().geometry_traits()->
            compare_endpoints_xy_2_object()(cv);

    if (he_dir == cv_dir) {
      (*m_edges_hash)[he] = cv.data().bc();
      (*m_edges_hash)[he->twin()] = cv.data().twin_bc();
    }
    else {
      (*m_edges_hash)[he] = cv.data().twin_bc();
      (*m_edges_hash)[he->twin()] = cv.data().bc();
    }
  }
};

template <typename Helper_, typename Arrangement_, typename Visitor_ = Default>
class Gps_agg_op_visitor :
  public Gps_agg_op_base_visitor<Helper_, Arrangement_,
                                 Gps_agg_op_visitor<Helper_, Arrangement_,
                                                    Visitor_> >
{
public:
  typedef Helper_                                       Helper;
  typedef Arrangement_                                  Arrangement_2;

  typedef typename Helper::Geometry_traits_2            Geometry_traits_2;
  typedef typename Helper::Event                        Event;
  typedef typename Helper::Subcurve                     Subcurve;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef Arrangement_2                                 Arr;

  typedef Gps_agg_op_visitor<Helper, Arr, Visitor_>     Self;
  typedef typename Default::Get<Visitor_, Self>::type   Visitor;
  typedef Gps_agg_op_base_visitor<Helper, Arr, Visitor> Base;

public:
  typedef typename Base::Halfedge_handle                Halfedge_handle;
  typedef typename Base::Vertex_handle                  Vertex_handle;
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Gt2::Point_2                         Point_2;

protected:
  unsigned int m_event_count;                   // The number of events so far.
  std::vector<Vertex_handle>* m_vertices_vec;   // The vertices, sorted in
                                                // ascending order.

public:
  Gps_agg_op_visitor(Arr* arr, typename Base::Edges_hash* hash,
                     std::vector<Vertex_handle>* vertices_vec) :
    Base(arr, hash),
    m_event_count(0),
    m_vertices_vec(vertices_vec)
  {}

  void before_handle_event(Event* event)
  {
    event->set_index(m_event_count);
    m_event_count++;
  }

  virtual Halfedge_handle
  insert_in_face_interior(const X_monotone_curve_2& cv, Subcurve* sc)
  {
    Halfedge_handle res_he = Base::insert_in_face_interior(cv, sc);

    // We now have a halfedge whose source vertex is associated with the
    // last event and whose target vertex is associated with the current event:
    Event* curr_event = this->current_event();
    Event* last_event = (sc)->last_event();

    CGAL_assertion((Arr_halfedge_direction)res_he->direction() ==
                    ARR_LEFT_TO_RIGHT);
    _insert_vertex(curr_event, res_he->target());
    _insert_vertex(last_event, res_he->source());

    return res_he;
  }

  virtual Halfedge_handle insert_from_right_vertex(const X_monotone_curve_2& cv,
                                                   Halfedge_handle he,
                                                   Subcurve* sc)
  {
    Halfedge_handle res_he = Base::insert_from_right_vertex(cv, he, sc);

    // We now have a halfedge whose target vertex is associated with the
    // last event (we have already dealt with its source vertex).
    Event* last_event = (sc)->last_event();
    CGAL_assertion((Arr_halfedge_direction)res_he->direction() ==
                    ARR_RIGHT_TO_LEFT);
    _insert_vertex(last_event, res_he->target());
    return res_he;
  }

  virtual Halfedge_handle insert_from_left_vertex(const X_monotone_curve_2& cv,
                                                  Halfedge_handle he,
                                                  Subcurve* sc)
  {
    Halfedge_handle  res_he = Base::insert_from_left_vertex(cv, he, sc);

    // We now have a halfedge whose target vertex is associated with the
    // current event(we have already dealt with its source vertex).
    Event* curr_event = this->current_event();

    CGAL_assertion((Arr_halfedge_direction)res_he->direction() ==
                    ARR_LEFT_TO_RIGHT);
    _insert_vertex (curr_event, res_he->target());
    return res_he;
  }

private:
  void _insert_vertex(const Event* event, Vertex_handle v)
  {
    const unsigned int index = event->index();
    if (index >= m_vertices_vec->size()) m_vertices_vec->resize(2 * (index + 1));
    (*m_vertices_vec)[index] = v;
  }

};

} // namespace CGAL

#endif
