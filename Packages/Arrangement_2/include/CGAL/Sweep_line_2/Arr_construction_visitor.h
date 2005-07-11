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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef ARR_CONSTRUCTION_VISITOR_H
#define ARR_CONSTRUCTION_VISITOR_H


#include <CGAL/Arrangement_2/Arr_accessor.h>
#include <CGAL/Sweep_line_2_empty_visitor.h>

CGAL_BEGIN_NAMESPACE

template <class _Traits, class Arr, class Event_, class Subcurve_> 
class Arr_construction_visitor : public Empty_visitor<_Traits,
                                                    Subcurve_,
                                                    Event_>
{
protected:

  typedef typename Arr::Halfedge_handle        Halfedge_handle;
  typedef typename Arr::Vertex_handle          Vertex_handle;
  typedef typename Arr::Face_handle            Face_handle;
  typedef Arr_construction_visitor< _Traits,
                                  Arr,
                                  Event_,
                                  Subcurve_>                 Self;

  typedef  _Traits                                          Traits;
  typedef typename Traits::X_monotone_curve_2               X_monotone_curve_2;
  typedef typename Traits::Point_2                          Point_2;

  typedef Event_ Event;
  typedef Subcurve_ Subcurve;
  typedef Empty_visitor<Traits, Subcurve, Event >         Base;
   
  typedef typename Base::SubCurveIter                  SubCurveIter;
  typedef typename Base::SubCurveRevIter               SubCurveRevIter;
  typedef typename Base::SL_iterator                   SL_iterator;


private:

  Arr_construction_visitor (const Self& );
  Self& operator= (const Self& );

public:

  Arr_construction_visitor(Arr *arr):
      m_arr(arr),
      m_arr_access (*arr)
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
      // the previous event on the curve is already in the planar map. 
      // Let's use it.
    {
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
    }
    this->current_event()->set_halfedge_handle(res);

    if(lastEvent->dec_right_curves_counter() == 0)
    {
      (this ->deallocate_event(lastEvent));
    }
  }

 


  virtual Halfedge_handle insert_in_face_interior(const X_monotone_curve_2& cv,
                                          Subcurve* sc)
  {
    return m_arr->insert_in_face_interior(cv, m_arr->unbounded_face());
  }

  virtual Halfedge_handle insert_at_vertices(const X_monotone_curve_2& cv,
                                             Halfedge_handle hhandle,
                                             Halfedge_handle prev,
                                             Subcurve* sc,
                                             bool &new_face_created)
  {
    Halfedge_handle res = 
      m_arr_access.insert_at_vertices_ex (cv, hhandle, prev,
			                                    new_face_created);

    if (new_face_created)
    {
      // In case a new face has been created (pointed by the new halfedge
      // we obtained), we have to examine the holes and isolated vertices
      // in the existing face (pointed be the twin halfedge) and relocate
      // the relevant features in the new face.
      m_arr_access.relocate_in_new_face (res);
    }

    return res;
  }

  virtual Halfedge_handle insert_from_right_vertex
                          (const X_monotone_curve_2& cv,
                           Halfedge_handle he,
                           Subcurve* sc)
  {
    return m_arr->insert_from_right_vertex(cv, he);
  }

  virtual Halfedge_handle insert_from_left_vertex
                          (const X_monotone_curve_2& cv,
                           Halfedge_handle he,
                           Subcurve* sc)
  {
    return m_arr->insert_from_left_vertex(cv, he);
  }


  virtual Vertex_handle insert_isolated_vertex(const Point_2& pt,
                                               SL_iterator iter)
  {
    return (m_arr ->insert_isolated_vertex(pt, m_arr->unbounded_face()));
  }



   protected:
     
     
  Arr               *m_arr;
  Arr_accessor<Arr>  m_arr_access;
};



CGAL_END_NAMESPACE

#endif
