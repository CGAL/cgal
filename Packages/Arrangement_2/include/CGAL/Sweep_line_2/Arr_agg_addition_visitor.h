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

#ifndef ARR_AGG_ADDITION_VISITOR_H
#define ARR_AGG_ADDITION_VISITOR_H

CGAL_BEGIN_NAMESPACE

template<class Traits, class Arr, class Event,class Subcurve>
class Arr_agg_addition_visitor : 
  public Arr_sweep_line_visitor<Traits,Arr,Event,Subcurve>
{
protected:

  typedef Arr_sweep_line_visitor<Traits,Arr,Event,Subcurve>    Base;
  typedef Arr_agg_addition_visitor<Traits,Arr,Event,Subcurve>  Self;

  typedef typename Base::Sweep_line                        Sweep_line;
  typedef typename Base::StatusLineIter                    StatusLineIter;
  typedef typename Base::Halfedge_handle                   Halfedge_handle;
  typedef Arr_insert_info<Halfedge_handle>                 ArrInsertInfo;
  typedef std::list<Subcurve *>                            SubcurveContainer;
  typedef typename SubcurveContainer::iterator             SubCurveIter;
  typedef typename SubcurveContainer::reverse_iterator     SubCurveRevIter;
  typedef typename Traits::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Traits::Point_2                         Point_2;
  typedef Sweep_line_2_impl<Traits,
                            Event,
                            Subcurve,
                            Self,
                            CGAL_ALLOCATOR(int)>           Sweep_line;

  typedef typename Arr::Face_handle                        Face_handle;
  typedef typename Arr::Face_const_handle                  Face_const_handle;
  
  using Base::m_currentEvent;
  using Base::m_arr;

public:

  Arr_agg_addition_visitor(Arr *arr): Base(arr)
  {}

  

  void before_handle_event(Event* event)
  {
    Base::before_handle_event(event);
    event->get_is_curve_in_arr().resize(event->get_num_right_curves(),false);

    if(!event->has_right_curves())
    {
      // update the event with the highest left halfedge
      for(SubCurveRevIter iter = event->left_curves_rbegin();
          iter != event->left_curves_rend();
          ++iter)
      {
        Halfedge_handle he;
        if((he =(*iter)->get_last_curve().get_halfedge_handle()) !=
         Halfedge_handle(NULL))
        {
          event->get_insert_info()->set_halfedge_handle(he.twin());
          return;
        }
      }
    }

    if(!event->has_left_curves())
    {
      int i = 0;
      // indicates if there's halfedge to the right of the event
      for(SubCurveRevIter iter = event->right_curves_rbegin();
          iter != event->right_curves_rend();
          ++iter, ++i)
      {
        // update the event with the highest right halfedge
        Halfedge_handle he;
        if((he = (*iter)->get_last_curve().get_halfedge_handle()) !=
           Halfedge_handle(NULL))
        {
          event->get_is_curve_in_arr()[i] = true;
          if(event->get_insert_info()->get_halfedge_handle() ==
             Halfedge_handle(NULL))
            event->get_insert_info()->set_halfedge_handle(he);
        }
      }
      return;
    }

    // the event has left and right curves
    int i = 0;
    bool exist_right_halfedge = false; 
    for(SubCurveRevIter iter = event->right_curves_rbegin();
        iter != event->right_curves_rend();
        ++iter, ++i)
    {
      Halfedge_handle he;
      if((he = (*iter)->get_last_curve().get_halfedge_handle()) !=
         Halfedge_handle(NULL))
      {
        exist_right_halfedge = true;
        event->get_is_curve_in_arr()[i] = true;
        /*if(reinterpret_cast<Event*>((*iter)->get_left_event()) ==
           m_currentEvent)*/
        if(!is_split_event(*iter, event))
          // halfedge will not be splitted 
          event->get_insert_info()->set_halfedge_handle(he);
        else
        {
          he = split_edge((*iter)->get_last_curve().get_halfedge_handle(),
                           event->get_point());
          
          // 'he' has the same source as the splitted halfedge
          event->get_insert_info()->set_halfedge_handle(he);
          X_monotone_curve_2& last_curve =
            const_cast<X_monotone_curve_2&>((*iter)->get_last_curve());
          last_curve.set_halfedge_handle(he);
          
          //there cannot be another existing halfedge that need to be splitted
          // because they are disjoint
          return;
        }
      }
    }

    if(exist_right_halfedge)
    {
      return;
    }
    // if we have reached here, there are no halfedges to the right of 
    // the event, but still can be on the left of the event
    for(SubCurveRevIter iter = event->left_curves_rbegin();
        iter != event->left_curves_rend();
        ++iter)
    {
      Halfedge_handle he;
      if((he =(*iter)->get_last_curve().get_halfedge_handle()) !=
        Halfedge_handle(NULL))
      {
        event->get_insert_info()->set_halfedge_handle(he.twin());
        return;
      }
    }
  }


  void add_subcurve(const X_monotone_curve_2& cv, Subcurve* sc)
  {
    if(cv.get_halfedge_handle() == Halfedge_handle(NULL))
      Base::add_subcurve(cv,sc);
    else
    {
      // sc is an overlap Subcurve of existing edge and new curve,
      // which means that the edeg will have to be modified
      if(sc -> get_orig_subcurve1())
        m_arr -> modify_edge(
          m_currentEvent->get_insert_info()->get_halfedge_handle().next().twin(),
          cv);

      Halfedge_handle next_ccw_he = 
        m_currentEvent->get_insert_info()->get_halfedge_handle().next().twin();
                                                                
      m_currentEvent->get_insert_info()->set_halfedge_handle(next_ccw_he);
    }
  }

  bool after_handle_event(Event* event, StatusLineIter iter, bool flag)
  {
    return (Base::after_handle_event(event,iter,flag));
  }


   Halfedge_handle insert_in_face_interior(const X_monotone_curve_2& cv,
                                          Subcurve* sc)
  {
    Halfedge_handle he_above;
    for(StatusLineIter iter = sc->get_hint();
        iter != (static_cast<Sweep_line*>(m_sweep_line))->StatusLine_end();
        ++iter)
    {
      if((*iter)->get_last_curve().get_halfedge_handle() !=
         Halfedge_handle(NULL))
      {
        he_above = (*iter)->get_last_curve().get_halfedge_handle();
        break;
      }
    }
    if(he_above == Halfedge_handle(NULL))
      return m_arr->insert_in_face_interior(cv,  m_arr->unbounded_face());
    return m_arr->insert_in_face_interior(cv, he_above.face());
  }

  Halfedge_handle insert_at_vertices(const X_monotone_curve_2& cv,
                                     Halfedge_handle hhandle,
                                     Halfedge_handle prev)
  {
    return (m_arr->insert_at_vertices(cv,hhandle,prev));
  }

  Halfedge_handle split_edge(Halfedge_handle he , const Point_2& pt)
  {
    // make sure that the halfedge associated with 'sc' is the directed from
    // right to left , since we always 'look' above , and the incident face 
    //is on the left of the  halfedge

    CGAL_assertion(m_traits.compare_xy_2_object()(he.source().point(),
      he.target().point()) == LARGER);

    X_monotone_curve_2 a, b;
    m_traits.split_2_object()(he.curve(), pt, b, a);
    //return m_arr->split_edge(he,a,b);
    return m_arr_access.split_edge_ex(he,pt, a, b);
  }

  // check if the halfedge associated with 'sc' will be splitted at the given
  // event point (its a recursive function since the hirearchy of potential 
  // ovrlap Subcuves
  bool is_split_event(Subcurve* sc, Event* event)
  {
    if(sc->get_last_curve().get_halfedge_handle() == Halfedge_handle(NULL))
      return false;
    if(! sc->get_orig_subcurve1())
    {
      return (reinterpret_cast<Event*>(sc->get_left_event())!=m_currentEvent);
    }
    return
      (is_split_event(reinterpret_cast<Subcurve*>(sc->get_orig_subcurve1()),
                      event) || 
       is_split_event(reinterpret_cast<Subcurve*>(sc->get_orig_subcurve2()),
                      event));
  }

 protected:
   Traits m_traits;



};


CGAL_END_NAMESPACE

#endif
