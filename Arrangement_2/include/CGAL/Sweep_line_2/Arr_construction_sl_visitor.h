// Copyright (c) 1997  Tel-Aviv University (Israel).
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
// $URL: svn+ssh://golubevs@scm.gforge.inria.fr/svn/cgal/trunk/Arrangement_2/include/CGAL/Sweep_line_2/Arr_construction_visitor.h $
// $Id: Arr_construction_visitor.h 37149 2007-03-16 09:16:03Z afabri $
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARR_CONSTRUCTION_SL_VISITOR_H
#define CGAL_ARR_CONSTRUCTION_SL_VISITOR_H


#include <CGAL/Arr_accessor.h>
#include <CGAL/Sweep_line_empty_visitor.h>
#include <CGAL/Unique_hash_map.h> 
#include <CGAL/Sweep_line_2/Integer_hash_function.h>
#include <vector>

CGAL_BEGIN_NAMESPACE

template <class Helper_>
class Arr_construction_sl_visitor :
  public Helper_::Base_visitor

{
protected:
  typedef Helper_                                      Helper;

  typedef typename Helper::Traits_2                    Traits;
  typedef typename Helper::Arrangement_2               Arrangement;
  typedef typename Helper::Base_visitor                Base;
  typedef typename Helper::Event                       Event;
  typedef typename Helper::Subcurve                    Subcurve;

  typedef typename Arrangement::Halfedge_handle        Halfedge_handle;
  typedef typename Arrangement::Vertex_handle          Vertex_handle;
  typedef typename Arrangement::Face_handle            Face_handle;
  typedef typename Arrangement::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
  typedef typename Arrangement::Ccb_halfedge_circulator 
    Ccb_halfedge_circulator;
  typedef typename Arrangement::Inner_ccb_iterator         Hole_iterator;

  typedef Arr_construction_sl_visitor<Helper>          Self;

  typedef typename Traits::X_monotone_curve_2          X_monotone_curve_2;
  typedef typename Traits::Point_2                     Point_2;
 
  typedef typename Base::SubCurveIter                  Event_subcurve_iterator;//SubCurveIter;
  typedef typename Base::SubCurveRevIter               Event_subcurve_reverse_iterator;//SubCurveRevIter;
  typedef typename Base::Status_line_iterator          SL_iterator;

  typedef Unique_hash_map<Halfedge_handle, 
                          std::list<unsigned int> >    Halfedge_indexes_map;
  typedef typename Halfedge_indexes_map::data_type     Indexes_list;
  typedef Unique_hash_map<unsigned int,
                          Vertex_handle,
                          Integer_hash_function>       Iso_vertices_map;

protected:
          
  Arr_accessor<Arrangement>  m_arr_access;
  unsigned int               m_sc_counter;  // Counter for subcurves that may
                                            // represent a hole (the upper
                                            // subcurves that emarge from event
                                            // points with only right curves). 

  std::vector<Halfedge_handle>
                             m_sc_he_table; // A table that maps a subcurve
                                            // index to its halfedhe handle,
                                            // directed from right ot left.
  Iso_vertices_map           m_iso_verts_map; // maps an index to the isolated vertex.

  Halfedge_indexes_map       m_he_indexes_table;

  // additional data members to support construction of unbounded arrangement.
  Halfedge_handle            m_lh;
  Halfedge_handle            m_th;
  Halfedge_handle            m_bh;
  Halfedge_handle            m_rh;
  std::list<unsigned int>    m_subcurves_at_ubf;
  Event*                     m_prev_minus_inf_x_event;
  Event*                     m_prev_plus_inf_y_event;
  
private:

  Arr_construction_sl_visitor (const Self& );
  Self& operator= (const Self& );

public:

  Arr_construction_sl_visitor(Arrangement *arr):
      m_arr_access (*arr),
      m_sc_counter (0),
      m_sc_he_table(1),
      m_prev_minus_inf_x_event(NULL),
      m_prev_plus_inf_y_event(NULL)
  {}

  virtual ~Arr_construction_sl_visitor(){}

  void before_sweep()
  {
    m_lh = m_arr_access.top_left_fictitious_vertex()->incident_halfedges();
    if(m_lh->source() != m_arr_access.bottom_left_fictitious_vertex())
      m_lh = m_lh->next();
    else
      m_lh = m_lh->twin();

    m_bh = m_lh->next();
    m_rh = m_bh->next();
    m_th = m_rh->next();

    CGAL_assertion(m_lh->direction() == LARGER);
    CGAL_assertion(m_lh->face() != m_arr_access.fictitious_face());
    CGAL_assertion(m_lh->source() == m_arr_access.top_left_fictitious_vertex() &&
                   m_lh->target() == m_arr_access.bottom_left_fictitious_vertex());

    CGAL_assertion(m_bh->direction() == SMALLER);
    CGAL_assertion(m_bh->face() != m_arr_access.fictitious_face());
    CGAL_assertion(m_bh->source() == m_arr_access.bottom_left_fictitious_vertex() &&
                   m_bh->target() == m_arr_access.bottom_right_fictitious_vertex());

    CGAL_assertion(m_rh->direction() == SMALLER);
    CGAL_assertion(m_rh->face() != m_arr_access.fictitious_face());
    CGAL_assertion(m_rh->source() == m_arr_access.bottom_right_fictitious_vertex() &&
                   m_rh->target() == m_arr_access.top_right_fictitious_vertex());

    CGAL_assertion(m_th->direction() == LARGER);
    CGAL_assertion(m_th->face() != m_arr_access.fictitious_face());
    CGAL_assertion(m_th->source() == m_arr_access.top_right_fictitious_vertex() &&
                   m_th->target() == m_arr_access.top_left_fictitious_vertex());
  }

  bool after_handle_event(Event* event, SL_iterator iter, bool)
  {
    if(!event->has_left_curves() && !event->has_right_curves())
    {
      //isolated event (no curves)
      Vertex_handle v = insert_isolated_vertex(event->point(), iter);
      m_iso_verts_map[++m_sc_counter] = v;
      insert_index_to_sc_he_table(m_sc_counter, Halfedge_handle());
      if(iter != this->status_line_end())
      {
        Subcurve *sc_above = *iter;
        sc_above->push_back_halfedge_index(m_sc_counter);
      }
      else
        m_subcurves_at_ubf.push_back(m_sc_counter);
      

      return true;
    }

    if(!event->has_left_curves() && event->is_finite())
    {
      CGAL_assertion(event->has_right_curves());
      //its an event that may represent a hole
      (*(event->right_curves_rbegin()))->set_index(++m_sc_counter);
      if(iter != this->status_line_end())
      {
        Subcurve *sc_above = *iter;
        sc_above->push_back_halfedge_index(m_sc_counter);
      }
      else
        m_subcurves_at_ubf.push_back(m_sc_counter);
    }

    for(Event_subcurve_iterator itr = event->left_curves_begin();
        itr != event->left_curves_end();
        ++itr)
    {
      (*itr)->set_last_event(event);
    }

    if(event->number_of_right_curves() == 0)
    {
      set_prev_inf_event_to_null(event);
      return true;
    }

    event->get_is_curve_in_arr().resize(event->number_of_right_curves(),false);
    for(Event_subcurve_iterator itr = event->right_curves_begin();
      itr != event->right_curves_end();
      ++itr)
    {
      (*itr)->set_last_event(event);
    }
    return false; 
  }

  void add_subcurve(const X_monotone_curve_2& cv,Subcurve* sc)
  {
    Event *lastEvent = get_last_event(sc);
    Halfedge_handle res; 
    Halfedge_handle hhandle = this ->current_event()->get_halfedge_handle();

    int jump = lastEvent->get_halfedge_jump_count(sc);
    // if the previous event on the curve is not in the planar map yet
    if ( lastEvent->get_halfedge_handle() == Halfedge_handle(NULL) ) 
    {
      // we have a handle from the previous insert
      if ( hhandle != Halfedge_handle(NULL) )
      { 
        res = this->insert_from_right_vertex(cv, hhandle, sc);
        res = res->twin();
      }
      else
      {
        // if this is the first left curve being inserted
        res = this->insert_in_face_interior(cv, sc);
      }
    } 
    else 
    {
      // the previous event on the curve is already in the planar map. 
      // Let's use it.
      Halfedge_handle prev = lastEvent->get_halfedge_handle();

      // skip to the right halfedge
      for ( int i = 0 ; i < jump ; i++ )
        prev = (prev->next())->twin();
      
      // we have a handle from the previous insert
      if ( hhandle != Halfedge_handle(NULL) ) 
      {
        CGAL_assertion(prev->face() == hhandle->face());
       
        bool dummy;
        

        res = this->insert_at_vertices(cv,hhandle,prev,sc, dummy);
        res = res->twin();
      }
      else
      {
        res = this->insert_from_left_vertex(cv, prev, sc);
      }
    }
    if ( lastEvent->number_of_left_curves() == 0 &&  
      lastEvent->is_curve_largest((Subcurve*)sc) )
    {
      lastEvent->set_halfedge_handle(res->twin());

      // if sc has valid index, insert his index to m_sc_he_table
      if(sc->has_valid_index())
      {
        CGAL_assertion(res->twin()->direction() == LARGER);
        insert_index_to_sc_he_table(sc->index(), res->twin());
      }
      
    }
    this->current_event()->set_halfedge_handle(res);

    if(lastEvent->dec_right_curves_counter() == 0)
    {
      set_prev_inf_event_to_null(lastEvent);
      (this ->deallocate_event(lastEvent));
    }

    //clear the list of indexes of sc
    sc->clear_haldedges_indexes();
  }

 
  virtual Halfedge_handle
    insert_in_face_interior(const X_monotone_curve_2& cv,
			    Subcurve* sc)
  {
    Vertex_handle v1 = 
      m_arr_access.create_vertex(_point(get_last_event(sc)->point()));
    Vertex_handle v2 =
      m_arr_access.create_vertex(_point(this->current_event()->point()));
    Halfedge_handle res =
      m_arr_access.insert_in_face_interior_ex(_curve(cv),
                                              m_th->face(),
                                              v1,
                                              v2,
                                              SMALLER);                                                                  
    if(sc->has_haldedges_indexes())
    {
      CGAL_assertion(res->twin()->direction() == LARGER);
      Indexes_list& list_ref = m_he_indexes_table[res->twin()];
      list_ref.clear();
      list_ref.splice(list_ref.end(), sc->get_haldedges_indexes_list());
    }
    return (res);

  }

  virtual Halfedge_handle insert_at_vertices(const X_monotone_curve_2& cv,
                                             Halfedge_handle hhandle,
                                             Halfedge_handle prev,
                                             Subcurve* sc,
                                             bool &new_face_created)
  {
    const bool      both_unbounded = hhandle->is_fictitious() || prev->is_fictitious();
    Halfedge_handle res;
    bool flip_res = false;
    if(this->current_event()->is_finite_in_x() && 
       this->current_event()->is_plus_boundary_in_y())
    {
      res = m_arr_access.insert_at_vertices_ex(_curve(cv),
                                               prev,
                                               hhandle,
                                               SMALLER,
                                               new_face_created,
                                               both_unbounded);
      flip_res = true;
    }
    else
      res = m_arr_access.insert_at_vertices_ex(_curve(cv),
                                               hhandle,
                                               prev,
                                               LARGER,
                                               new_face_created,
                                               both_unbounded);

     // map the halfedge to the indexes list of all subcurves that are below him
      if(sc->has_haldedges_indexes())
      {
        Halfedge_handle temp = res;
        if(flip_res)
          temp = temp->twin();
        CGAL_assertion(temp->direction() == LARGER);
        Indexes_list& list_ref = m_he_indexes_table[temp];
        list_ref.clear();
        list_ref.splice(list_ref.end(), sc->get_haldedges_indexes_list());
      }
   
    if (new_face_created)
    {
      // In case a new face has been created (pointed by the new halfedge
      // we obtained), we have to examine the holes and isolated vertices
      // in the existing face (pointed be the twin halfedge) and relocate
      // the relevant features in the new face.
      //m_arr_access.relocate_in_new_face (res);
      CGAL_assertion(res->face() != res->twin()->face());
      
      this->relocate_holes_and_iso_verts_in_new_face(res);
    }

    if(flip_res)
      return res->twin();

    return res;
  }

  virtual Halfedge_handle insert_from_right_vertex
                          (const X_monotone_curve_2& cv,
                           Halfedge_handle he,
                           Subcurve* sc)
  {
    Vertex_handle v = 
      m_arr_access.create_vertex(_point(get_last_event(sc)->point()));
    Halfedge_handle res = m_arr_access.insert_from_vertex_ex(_curve(cv), he, v, LARGER);
    if(sc->has_haldedges_indexes())
    {
      CGAL_assertion(res->direction() == LARGER);
      Indexes_list& list_ref = m_he_indexes_table[res];
      list_ref.clear();
      list_ref.splice(list_ref.end(), sc->get_haldedges_indexes_list());
    }
    return (res);
  }

  virtual Halfedge_handle insert_from_left_vertex
                          (const X_monotone_curve_2& cv,
                           Halfedge_handle he,
                           Subcurve* sc)
  {
    Vertex_handle v = 
      m_arr_access.create_vertex(_point(this->current_event()->point()));
    Halfedge_handle res = m_arr_access.insert_from_vertex_ex(_curve(cv), he, v, SMALLER);
    if(sc->has_haldedges_indexes())
    {
      CGAL_assertion(res->twin()->direction() == LARGER);
      Indexes_list& list_ref = m_he_indexes_table[res->twin()];
      list_ref.clear();
      list_ref.splice(list_ref.end(), sc->get_haldedges_indexes_list());
    }
    
    return (res);
  }


  virtual Vertex_handle insert_isolated_vertex(const Point_2& pt,
                                               SL_iterator)
  {
    return (m_arr_access.arrangement().insert_in_face_interior (_point(pt),
								m_th->face()));
  }

  void relocate_holes_and_iso_verts_in_new_face(Halfedge_handle he)
  {
    // We use a constant indexes map so no new entries are added there.
    const Halfedge_indexes_map& const_he_indexes_table = m_he_indexes_table;
    Face_handle                 new_face = he->face();
    Halfedge_handle             curr_he = he;

    do
    {
      // we are intreseted only in halfedges directed from right to left.
      if(curr_he->direction() == LARGER)
      {
        const Indexes_list& indexes_list = const_he_indexes_table[curr_he];
        typename Indexes_list::const_iterator itr;

        for (itr = indexes_list.begin();
             itr != indexes_list.end();
             ++itr)
        {
          CGAL_assertion(*itr != 0 && *itr < m_sc_he_table.size());
          Halfedge_handle he_on_face = m_sc_he_table[*itr];
          //if he_on_face is a null halfedge handle then its index for an
          //isolated vertex.
          if(he_on_face == Halfedge_handle())
          {
            Vertex_handle v = m_iso_verts_map[*itr];
            CGAL_assertion(v != Vertex_handle());
            if(v->face() == new_face)
              continue;
            m_arr_access.move_isolated_vertex(v->face(),
                                              new_face,
                                              v);

          }
          else
          {
            if(he_on_face->twin()->face() == new_face)
              //this hole was already relocated
              continue;

            m_arr_access.move_inner_ccb (he_on_face->twin()->face(),
                                    new_face,
                                    he_on_face->twin()->ccb());
            relocate_holes_and_iso_verts_in_new_face(he_on_face->twin());
          }
        }
      }
      curr_he = curr_he->next();
    }
    while(curr_he != he);

  }

  void before_handle_event(Event* event)
  {
    if(event->is_finite())
      return;

     // if it is an event at infinity, split the corresponding fictitious edge.
    Boundary_type inf_x = event->infinity_at_x();
    Boundary_type inf_y = event->infinity_at_y();

    Vertex_handle v_at_inf = 
      m_arr_access.create_vertex_at_infinity(inf_x, inf_y);
    switch(inf_x)
    {
    case MINUS_INFINITY:
      m_arr_access.split_fictitious_edge(m_lh, v_at_inf);
      event->set_halfedge_handle(m_lh);
      if(m_prev_minus_inf_x_event)
        m_prev_minus_inf_x_event->set_halfedge_handle(m_lh->next());
      m_prev_minus_inf_x_event = event;
      return;

    case PLUS_INFINITY:
      m_arr_access.split_fictitious_edge(m_rh, v_at_inf);
      event->set_halfedge_handle(m_rh);
      m_rh = m_rh->next();
      return;

    case NO_BOUNDARY:
    default:
      break;
    }

    switch(inf_y)
    {
      case MINUS_INFINITY:
        m_arr_access.split_fictitious_edge(m_bh, v_at_inf);
        event->set_halfedge_handle(m_bh);
        m_bh = m_bh->next();
        return;

      case PLUS_INFINITY:
        {
          m_arr_access.split_fictitious_edge(m_th, v_at_inf);
          event->set_halfedge_handle(m_th);
          if(m_prev_plus_inf_y_event != NULL)
            m_prev_plus_inf_y_event->set_halfedge_handle(m_th->next());
          m_prev_plus_inf_y_event = event;
          Indexes_list& list_ref = m_he_indexes_table[m_th->next()];
          list_ref.clear();
          list_ref.splice(list_ref.end(), m_subcurves_at_ubf);
          CGAL_assertion(m_subcurves_at_ubf.empty());
        }
        return;

      case NO_BOUNDARY:
      default:
        // doesn't suppose to reach here at all.
        CGAL_error();
    }
  }
  
  Event* get_last_event(Subcurve* sc)
  {
    return (reinterpret_cast<Event*>((sc)->last_event()));
  }

private:

  /*!
   * Cast a Traits::Point_2 into an Arrangement::Point_2 object.
   * These two types may not be the same when the addition visitor inherits
   * from this base class.
   */
  const typename Arrangement::Point_2& _point (const Point_2& p) const
  {
    return (static_cast<const typename Arrangement::Point_2&> (p));
  }

  /*!
   * Cast a Traits::X_monotone_curve_2 into an Arrangement::X_monotone_curve_2
   * object.
   * These two types may not be the same when the addition visitor inherits
   * from this base class.
   */
  const typename Arrangement::X_monotone_curve_2&
    _curve (const X_monotone_curve_2& cv) const
  {
    return (static_cast<const typename Arrangement::X_monotone_curve_2&> (cv));
  }

  void insert_index_to_sc_he_table(unsigned int i, Halfedge_handle he)
  {
    CGAL_assertion(i!=0);
    if(i >= m_sc_he_table.size())
    {
      m_sc_he_table.resize(2*i);
    }

    m_sc_he_table[i] = he;
  }

  void set_prev_inf_event_to_null(Event* e)
  {
    CGAL_assertion(!this->current_event()->is_minus_boundary_in_x());    
    if(e == m_prev_plus_inf_y_event)
      m_prev_plus_inf_y_event = NULL;
  }

};

CGAL_END_NAMESPACE

#endif
