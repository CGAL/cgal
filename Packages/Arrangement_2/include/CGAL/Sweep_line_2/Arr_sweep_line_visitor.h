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

#ifndef ARR_SWEEP_LINE_VISITOR_H
#define ARR_SWEEP_LINE_VISITOR_H


#include <CGAL/Sweep_line_2/Arr_insert_info.h>
#include <CGAL/Arrangement_2/Arr_accessor.h>

CGAL_BEGIN_NAMESPACE

template <class Arr, class Event, class Subcurve>
class Arr_sweep_line_visitor
{
protected:

  typedef typename Arr::Halfedge_handle                            Halfedge_handle;
  typedef Arr_sweep_line_visitor< Arr,
                                  Event,
                                  Subcurve>                        Self;

  typedef typename Arr::Traits_2                                   Traits;
  typedef typename Traits::X_monotone_curve_2                      X_monotone_curve_2;
  typedef typename Traits::Point_2                                 Point_2;

  typedef Sweep_line_2_impl<Traits,
                            Event,
                            Subcurve,
                            Self,
                            CGAL_ALLOCATOR(int)>           Sweep_line;
  typedef typename Sweep_line::StatusLineIter              StatusLineIter;

  typedef Arr_insert_info<Halfedge_handle>                        ArrInsertInfo;
  typedef std::list<Subcurve *>                                    SubcurveContainer;
  typedef typename SubcurveContainer::iterator                     SubCurveIter;

private:

  Arr_sweep_line_visitor (const Self& );
  Self& operator= (const Self& );

public:

  Arr_sweep_line_visitor(Arr *arr):
      m_arr(arr),
      m_arr_access (*arr)
  {}

  void attach(Sweep_line *sl)
  {
    m_sweep_line = sl;
  }

 

  void before_handle_event(Event* event)
  {
    m_currentEvent = event;
  }

  bool after_handle_event(Event* event, StatusLineIter iter, bool flag)
  {
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
    Event *lastEvent = (sc)->get_last_event();
    ArrInsertInfo *insertInfo = lastEvent->get_insert_info();
    Halfedge_handle res; 
    ArrInsertInfo *currentInfo = m_currentEvent -> get_insert_info();
    Halfedge_handle hhandle = currentInfo->get_halfedge_handle();

    int jump = lastEvent->get_halfedge_jump_count(sc);
   
    // if the previous event on the curve is not in the planar map yet
    if ( insertInfo->get_halfedge_handle() == Halfedge_handle(NULL) ) 
    {
      // we have a handle from the previous insert
      if ( hhandle != Halfedge_handle(NULL) )
      {
        res = m_arr->insert_from_vertex(cv, hhandle);
        res = res.twin();
      }
      else
      {
        // if this is the first left curve being inserted
        res = m_arr->insert_in_face_interior(cv, m_arr->unbounded_face());
      }
    } 
    else 
      // the previous event on the curve is already in the planar map. 
      // Let's use it.
    {
      Halfedge_handle prev = insertInfo->get_halfedge_handle();
     
      // skip to the right halfedge
      for ( int i = 0 ; i < jump ; i++ )
        prev = (prev.next()).twin();
      
      // we have a handle from the previous insert
      if ( hhandle != Halfedge_handle(NULL) ) 
      {
        CGAL_assertion(prev.face() == hhandle.face());
       
        //res = m_arr->insert_at_vertices(cv,hhandle,prev);
        bool      new_face_created;

        res = m_arr_access.insert_at_vertices_ex (cv, hhandle, prev,
					                                        new_face_created);

        if (new_face_created)
        {
          // In case a new face has been created (pointed by the new halfedge
          // we obtained), we have to examine the holes in the existing face
          // (pointed be the twin halfedge) and move the relevant holes into
          // the new face.
          typename Arr::Face_handle     new_face = res.face();
          typename Arr::Face_handle     old_face = res.twin().face();
          typename Arr::Holes_iterator  holes_it = old_face.holes_begin();
          typename Arr::Holes_iterator  hole_to_move;

          while (holes_it != old_face.holes_end())
          {
            // Check whether the current hole is inside new face.
            if (m_arr_access.point_is_in((*(*holes_it)).target().point(), res))
            {
              // We increment the holes itrator before moving the hole, because
              // this operation invalidates the iterator.
              hole_to_move  = holes_it;
              ++holes_it;

              // Move the hole.
              m_arr_access.move_hole (old_face, new_face, hole_to_move);
            }
            else
            {
              ++holes_it;
            }
          }
        }

        res = res.twin();
      }
      else
      {
        res = m_arr->insert_from_vertex(cv, prev);
      }
    }
    if ( lastEvent->get_num_left_curves() == 0 &&  
      lastEvent->is_curve_largest((Subcurve*)sc) )
    {
      insertInfo->set_halfedge_handle(res.twin());
    }
    currentInfo->set_halfedge_handle(res);

    if(lastEvent->get_insert_info()->dec_right_curves_counter() == 0)
    {
      m_sweep_line->deallocate_event(lastEvent);
    }
  }


  void init_subcurve(Subcurve* sc)
  {
    sc -> set_last_event((Event*)(sc->get_left_event()));
  }

 
   protected:
     
     
  Arr               *m_arr;
  Arr_accessor<Arr>  m_arr_access;
  Sweep_line        *m_sweep_line;
  Event             *m_currentEvent;
};



CGAL_END_NAMESPACE

#endif




