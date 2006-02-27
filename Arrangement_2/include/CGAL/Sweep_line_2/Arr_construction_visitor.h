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
// $URL$
// $Id$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARR_CONSTRUCTION_VISITOR_H
#define CGAL_ARR_CONSTRUCTION_VISITOR_H


#include <CGAL/Arr_accessor.h>
#include <CGAL/Sweep_line_2_empty_visitor.h>
#include <CGAL/Unique_hash_map.h> 
#include <vector>

CGAL_BEGIN_NAMESPACE

template <class Traits_, class Arrangement, class Event_, class Subcurve_> 
class Arr_construction_visitor : public Empty_visitor<Traits_,
                                                      Subcurve_,
                                                      Event_>
{
protected:

  typedef typename Arrangement::Halfedge_handle        Halfedge_handle;
  typedef typename Arrangement::Vertex_handle          Vertex_handle;
  typedef typename Arrangement::Face_handle            Face_handle;
  typedef typename Arrangement::Ccb_halfedge_circulator 
    Ccb_halfedge_circulator;
  typedef typename Arrangement::Hole_iterator         Hole_iterator;
  typedef Arr_construction_visitor<Traits_,
                                   Arrangement,
                                   Event_,
                                   Subcurve_>          Self;

  typedef Traits_                                      Traits;
  typedef Event_                                       Event;
  typedef Subcurve_                                    Subcurve;

  typedef typename Traits::X_monotone_curve_2          X_monotone_curve_2;
  typedef typename Traits::Point_2                     Point_2;
 
  typedef Empty_visitor<Traits,
                        Subcurve,
                        Event >                        Base;
 
  typedef typename Base::SubCurveIter                  SubCurveIter;
  typedef typename Base::SubCurveRevIter               SubCurveRevIter;
  typedef typename Base::SL_iterator                   SL_iterator;

private:

  Arr_construction_visitor (const Self& );
  Self& operator= (const Self& );

public:

  Arr_construction_visitor(Arrangement *arr):
      m_arr(arr),
      m_arr_access (*arr),
      m_sc_counter (0),
      m_sc_he_table(1)
  {}

  virtual ~Arr_construction_visitor(){}


  bool after_handle_event(Event* event, SL_iterator iter, bool flag)
  {
    if(!event->has_left_curves() && !event->has_right_curves())
    {
      //isolated event (no curves)
      insert_isolated_vertex(event->get_point(), iter);
      return true;
    }

    if(!event->has_left_curves())
    {
      CGAL_assertion(event->has_right_curves());
      //its an event that may represent a hole
      (*(event->right_curves_rbegin()))->set_index(++m_sc_counter);
      if(iter != this->status_line_end())
      {
        Subcurve *sc_above = *iter;
        sc_above->push_back_halfedge_index(m_sc_counter);
      }
    }

    for(SubCurveIter itr = event->left_curves_begin();
        itr != event->left_curves_end();
        ++itr)
    {
      (*itr)->set_last_event(event);
    }

    if(event->get_num_right_curves() == 0)
      return true;

    event->get_is_curve_in_arr().resize(event->get_num_right_curves(),false);
    for(SubCurveIter itr = event->right_curves_begin();
      itr != event->right_curves_end();
      ++itr)
    {
      (*itr)->set_last_event(event);
    }
    return false; 
  }

  void add_subcurve(const X_monotone_curve_2& cv,Subcurve* sc)
  {
    Event *lastEvent = reinterpret_cast<Event*>((sc)->get_last_event());
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
    if ( lastEvent->get_num_left_curves() == 0 &&  
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
      (this ->deallocate_event(lastEvent));
    }

    //clear the list of indexes of sc
    sc->clear_haldedges_indexes();
  }

 
  virtual Halfedge_handle
    insert_in_face_interior(const X_monotone_curve_2& cv,
			    Subcurve* sc)
  {
    Halfedge_handle res =  m_arr->insert_in_face_interior (_curve(cv),
					    m_arr->unbounded_face());
    if(sc->has_haldedges_indexes())
    {
      CGAL_assertion(res->twin()->direction() == LARGER);
      std::list<unsigned int>& list_ref = 
        (m_he_indexes_table[res->twin()] = std::list<unsigned int>());
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
    Halfedge_handle res =
      m_arr_access.insert_at_vertices_ex (_curve(cv),
                                          hhandle, prev,
                                          LARGER,
                                          new_face_created);
     // map the halfedge to the indexes list of all subcurves that are below him
      if(sc->has_haldedges_indexes())
      {
        CGAL_assertion(res->direction() == LARGER);
        std::list<unsigned int>& list_ref = 
          (m_he_indexes_table[res] = std::list<unsigned int>());
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
      CGAL_assertion(res->face() != m_arr->unbounded_face());
      
      this->relocate_holes_in_new_face(res);
      m_arr_access.relocate_isolated_vertices_in_new_face(res);
    }

    return res;
  }

  virtual Halfedge_handle insert_from_right_vertex
                          (const X_monotone_curve_2& cv,
                           Halfedge_handle he,
                           Subcurve* sc)
  {
    Halfedge_handle res = m_arr->insert_from_right_vertex (_curve(cv), he);
    if(sc->has_haldedges_indexes())
    {
      CGAL_assertion(res->direction() == LARGER);
      std::list<unsigned int>& list_ref = 
        (m_he_indexes_table[res] = std::list<unsigned int>());
      list_ref.splice(list_ref.end(), sc->get_haldedges_indexes_list());
    }
    return (res);
  }

  virtual Halfedge_handle insert_from_left_vertex
                          (const X_monotone_curve_2& cv,
                           Halfedge_handle he,
                           Subcurve* sc)
  {
    Halfedge_handle res = m_arr->insert_from_left_vertex (_curve(cv), he);
    if(sc->has_haldedges_indexes())
    {
      CGAL_assertion(res->twin()->direction() == LARGER);
      std::list<unsigned int>& list_ref = 
        (m_he_indexes_table[res->twin()] = std::list<unsigned int>());
      list_ref.splice(list_ref.end(), sc->get_haldedges_indexes_list());
    }
    
    return (res);
  }


  virtual Vertex_handle insert_isolated_vertex(const Point_2& pt,
                                               SL_iterator iter)
  {
    return (m_arr->insert_in_face_interior (_point(pt),
					    m_arr->unbounded_face()));
  }

  void relocate_holes_in_new_face(Halfedge_handle he)
  {
    const Unique_hash_map<Halfedge_handle, std::list<unsigned int> >& 
      const_he_indexes_table = m_he_indexes_table;
    Face_handle new_face = he->face();
    Halfedge_handle curr_he = he;
    do
    {
      // we are intreseted only in halfedges directed from right to left.
      if(curr_he->direction() == LARGER)
      {
        const std::list<unsigned int>&    indexes_list = 
                                         const_he_indexes_table[curr_he];
        std::list<unsigned int>::const_iterator itr;

        for (itr = indexes_list.begin();
             itr != indexes_list.end();
             ++itr)
        {
          CGAL_assertion(*itr != 0 && *itr < m_sc_he_table.size());
          Halfedge_handle he_on_face = m_sc_he_table[*itr];
          if(he_on_face->twin()->face() == new_face)
            //this hole was already relocated
            continue;

          m_arr_access.move_hole (he_on_face->twin()->face(),
                                  new_face,
                                  he_on_face->twin()->ccb());
          relocate_holes_in_new_face(he_on_face->twin());
        }
      }
      curr_he = curr_he->next();
    }
    while(curr_he != he);

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



protected:
          
  Arrangement*               m_arr;
  Arr_accessor<Arrangement>  m_arr_access;
  // counter for Subcurves that may  represent a hole (the upper sc 
  //emarge from  an event with only right curves
  unsigned int               m_sc_counter; 
  std::vector<Halfedge_handle>  m_sc_he_table;  //a table that maps index of a 
                                             // subcurve to his halfedhe handle
                                             // directed from right ot left.

  Unique_hash_map<Halfedge_handle, std::list<unsigned int> >  m_he_indexes_table;
};

CGAL_END_NAMESPACE

#endif
