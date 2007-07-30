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


#ifndef CGAL_ARR_BASIC_ADDITION_VISITOR_H
#define CGAL_ARR_BASIC_ADDITION_VISITOR_H

#include <CGAL/Arr_topology_traits/Arr_bounded_planar_construction_helper.h>

CGAL_BEGIN_NAMESPACE

template<class Traits, class Arrangement_, class Event,class Subcurve>
class Arr_basic_addition_visitor : 
  public  Arr_construction_sl_visitor< 
            Arr_bounded_planar_construction_helper<
                                                Traits, 
                                                Arrangement_, 
                                                Event, 
                                                Subcurve> >
{
protected:

  typedef Arrangement_                                     Arrangement;
  
  typedef Arr_bounded_planar_construction_helper<Traits, 
                                              Arrangement, 
                                              Event, 
                                              Subcurve> Construction_helper;
  typedef Arr_construction_sl_visitor<Construction_helper> Base;

  typedef Arr_basic_addition_visitor<Traits,
                               Arrangement,
                               Event,
                               Subcurve>                   Self;

  typedef typename Base::SL_iterator                       SL_iterator;
  typedef typename Base::Halfedge_handle                   Halfedge_handle;
  typedef typename Base::Vertex_handle                     Vertex_handle;
  typedef typename Base::SubCurveIter                      SubCurveIter;
  typedef typename Base::SubCurveRevIter                   SubCurveRevIter;
  typedef typename Traits::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Traits::Point_2                         Point_2;
  
  typedef typename Arrangement::Face_handle                Face_handle;
  typedef typename Arrangement::Face_const_handle          Face_const_handle;
  
public:

  Arr_basic_addition_visitor(Arrangement *arr): Base(arr)
  {}

  void before_sweep()
  {
    this->m_lh = 
      this->m_arr_access.bottom_left_fictitious_vertex()->incident_halfedges();
    if (this->m_lh->source()->boundary_in_x() != MINUS_INFINITY)
      this->m_lh = this->m_lh->next()->twin();
   
    this->m_bh = this->m_lh->next();
   
    this->m_th = 
      this->m_arr_access.top_left_fictitious_vertex()->incident_halfedges();
    if(this->m_th->source()->boundary_in_x() == MINUS_INFINITY)
      this->m_th = this->m_th->next()->twin();

    this->m_rh = 
      this->m_arr_access.bottom_right_fictitious_vertex()->incident_halfedges();
    if (this->m_rh->source()->boundary_in_x() == PLUS_INFINITY)
      this->m_rh = this->m_rh->twin();
    else
      this->m_rh = this->m_rh->next();

    CGAL_assertion (this->m_lh->direction() == LARGER);
    CGAL_assertion (this->m_lh->face() != 
                    this->m_arr_access.fictitious_face());
    CGAL_assertion (this->m_lh->target() == 
                    this->m_arr_access.bottom_left_fictitious_vertex());

    CGAL_assertion (this->m_bh->direction() == SMALLER);
    CGAL_assertion (this->m_bh->face() != 
                    this->m_arr_access.fictitious_face());
    CGAL_assertion (this->m_bh->source() ==
                    this->m_arr_access.bottom_left_fictitious_vertex());

    CGAL_assertion (this->m_rh->direction() == SMALLER);
    CGAL_assertion (this->m_rh->face() != 
                    this->m_arr_access.fictitious_face());
    CGAL_assertion (this->m_rh->source() == 
                    this->m_arr_access.bottom_right_fictitious_vertex());

    CGAL_assertion (this->m_th->direction() == LARGER);
    CGAL_assertion (this->m_th->face() != 
                    this->m_arr_access.fictitious_face());
    CGAL_assertion (this->m_th->target() == 
                    this->m_arr_access.top_left_fictitious_vertex());
  }

  void before_handle_event(Event* event)
  {
    if(!event->is_finite())
    {
       if(event->get_unbounded_curve().get_halfedge_handle() ==
          Halfedge_handle(NULL))
         Base::before_handle_event(event);
       else
       {
         Boundary_type x_inf = event->infinity_at_x();
         if(x_inf == MINUS_INFINITY)
         {
           this->m_lh = this->m_lh->twin()->next()->twin();
           this->m_prev_minus_inf_x_event = NULL;
         }
         else
           if(x_inf == PLUS_INFINITY)
             this->m_rh = this->m_rh->twin()->prev()->twin();
           else
           {
             Boundary_type y_inf = event->infinity_at_y();
             if(y_inf == MINUS_INFINITY)
               this->m_bh = this->m_bh->twin()->prev()->twin();
             else
             {
               CGAL_assertion(y_inf == PLUS_INFINITY);
               this->m_th = this->m_th->twin()->next()->twin();
               this->m_prev_plus_inf_y_event = NULL;
             }
           }
       }
    }

    event->get_is_curve_in_arr().resize(event->number_of_right_curves(),false);

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
          event->set_halfedge_handle(he->twin());
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
          if(event->get_halfedge_handle() ==
             Halfedge_handle(NULL))
            event->set_halfedge_handle(he);
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
        if(!is_split_event(*iter, event))
          // halfedge will not be splitted 
          event->set_halfedge_handle(he);
        else
        {
          he = split_edge((*iter)->get_last_curve().get_halfedge_handle(),
                          (*iter),
                           event->point());
          
          // 'he' has the same source as the splitted halfedge
          event->set_halfedge_handle(he);
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
        event->set_halfedge_handle(he->twin());
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
      if (sc->get_orig_subcurve1())
      {
        this->m_arr_access.arrangement().modify_edge
          (this->current_event()->get_halfedge_handle()->
           next()->twin(),
           cv.base());
      }

      Halfedge_handle next_ccw_he = 
        this->current_event()->get_halfedge_handle()->next()->twin();
                                                                
      this->current_event()->set_halfedge_handle(next_ccw_he);
    }
  }

  bool after_handle_event(Event* event, SL_iterator iter, bool flag)
  {
    return (Base::after_handle_event(event,iter,flag));
  }


  Halfedge_handle _insert_in_face_interior (const X_monotone_curve_2& cv,
                                            Subcurve* sc)
  {
    Halfedge_handle he_above = ray_shoot_up(sc);
    Vertex_handle v1 = 
      this->m_arr_access.create_vertex(this->get_last_event(sc)->point());
    Vertex_handle v2 =
      this->m_arr_access.create_vertex(this->current_event()->point());
    if(he_above == Halfedge_handle(NULL))
      return this->m_arr_access.insert_in_face_interior_ex(cv.base(),
                                                   this->m_th->face(),
                                                   v1,
                                                   v2,
                                                   SMALLER);

    return this->m_arr_access.insert_in_face_interior_ex(cv.base(),
                                                   he_above->face(),
                                                   v1,
                                                   v2,
                                                   SMALLER);
  }

  Halfedge_handle insert_at_vertices (const X_monotone_curve_2& cv,
                                      Halfedge_handle hhandle,
                                      Halfedge_handle prev,
                                      Subcurve* sc,
                                      bool &new_face_created)
  {
    return (_insert_at_vertices(cv, hhandle, prev, sc, new_face_created));
  }

  virtual Halfedge_handle split_edge(Halfedge_handle /*he*/,
                                     Subcurve* /*sc*/,
                                     const Point_2& /*pt*/)
  {
    return Halfedge_handle();
  }
  

  // check if the halfedge associated with 'sc' will be splitted at the given
  // return false.
  virtual bool is_split_event(Subcurve* /*sc*/, Event* /*event*/)
  {
    return false;
  }


  virtual Vertex_handle insert_isolated_vertex(const Point_2& pt,
                                               SL_iterator iter)
  {
    Vertex_handle res(NULL);
    //the isolated vertex is already at the arrangement
    if(pt.get_vertex_handle() != Vertex_handle(NULL))
      return res;
    if(iter == this->status_line_end())
    {
      res = this->m_arr_access.arrangement().insert_in_face_interior
        (pt.base(), this->m_arr_access.arrangement().unbounded_face());
    }
    else
    {
      Halfedge_handle he = ray_shoot_up(*iter);
      if (he == Halfedge_handle(NULL))
      {
        res = this->m_arr_access.arrangement().insert_in_face_interior
          (pt.base(), this->m_arr_access.arrangement().unbounded_face());
      }
      else
      {
        res = this->m_arr_access.arrangement().insert_in_face_interior (pt.base(), he->face());
      }
    }

    return (res);
  }


  Halfedge_handle ray_shoot_up(Subcurve* sc)
  {
    Halfedge_handle he_above;
    for(SL_iterator iter = this -> status_line_position(sc);
        iter != this -> status_line_end();
        ++iter)
    {
      if((*iter)->get_last_curve().get_halfedge_handle() !=
         Halfedge_handle(NULL))
      {
        he_above = (*iter)->get_last_curve().get_halfedge_handle();
        return (he_above);
      }
    }
  
    return (he_above);
  }




   void update_event(){}

   void update_event(Event*,
                    const Point_2&,
                    const X_monotone_curve_2&,
		     bool /* is_left_end */)
  {}

  void update_event(Event*,
                    Subcurve*,
                    Subcurve*,
                    bool /*created*/ = false)
  {}

  void update_event(Event*,
                    Subcurve*)
  {}


  void update_event(Event* e, const Point_2& pt)
  {
    Vertex_handle vh;
    if(e->point().get_vertex_handle() == vh)
      e->point().set_vertex_handle(pt.get_vertex_handle());
  }

  virtual Halfedge_handle
    insert_in_face_interior(const X_monotone_curve_2& cv,
                            Subcurve* sc)
  {
    Event *lastEvent = this->get_last_event(sc);
    Vertex_handle last_v = lastEvent->point().get_vertex_handle();
    Vertex_handle curr_v =
      this->current_event()->point().get_vertex_handle();
    Vertex_handle null_v;

    if(last_v == null_v && curr_v == null_v)
      return (this->_insert_in_face_interior(cv, sc));
    if(last_v == null_v && curr_v != null_v)
    {
      Halfedge_handle he = this->m_arr_access.arrangement().insert_from_right_vertex (cv.base(),
                                                                  curr_v);
      return (he->twin());
    }
    if(last_v != null_v && curr_v == null_v)
      return (this->m_arr_access.arrangement().insert_from_left_vertex (cv.base(), last_v));

    CGAL_assertion(last_v != null_v && curr_v != null_v);
    return (this->m_arr_access.arrangement().insert_at_vertices (cv.base(), last_v, curr_v));

  }


  virtual Halfedge_handle
    insert_from_right_vertex (const X_monotone_curve_2& cv,
                              Halfedge_handle he,
                              Subcurve* sc)
  {
    Event *lastEvent = this->get_last_event(sc);
    Vertex_handle last_v = lastEvent->point().get_vertex_handle();

    if(last_v != Vertex_handle())
      return (this->m_arr_access.arrangement().insert_at_vertices (cv.base(), he, last_v));
    
    return (_insert_from_right_vertex(cv, he, sc));
  }

  virtual Halfedge_handle
    insert_from_left_vertex (const X_monotone_curve_2& cv,
                             Halfedge_handle he,
                             Subcurve* sc)
  {
    Vertex_handle curr_v =
      this->current_event()->point().get_vertex_handle();

     if(curr_v != Vertex_handle())
       return (this->m_arr_access.arrangement().insert_at_vertices (cv.base(), he, curr_v));

    return (_insert_from_left_vertex(cv, he, sc));
  }

  Halfedge_handle _insert_from_left_vertex(const X_monotone_curve_2& cv,
                                           Halfedge_handle he,
                                           Subcurve* /*sc*/)
  {
    Event* curr_event = this->current_event();
    Vertex_handle v = 
        this->m_arr_access.create_vertex(curr_event->point().base());
      return this->m_arr_access.insert_from_vertex_ex(cv.base(), he, v, SMALLER);
  }

  Halfedge_handle _insert_from_right_vertex(const X_monotone_curve_2& cv,
                                            Halfedge_handle he,
                                            Subcurve* sc)
  {
    Event* last_event = this->get_last_event(sc);
    Vertex_handle v = 
        this->m_arr_access.create_vertex(last_event->point().base());
      return this->m_arr_access.insert_from_vertex_ex(cv.base(), he, v, LARGER);
  }

  Halfedge_handle _insert_at_vertices (const X_monotone_curve_2& cv,
                                       Halfedge_handle hhandle,
                                       Halfedge_handle prev,
                                       Subcurve* /*sc*/,
                                       bool &new_face_created)
  {
 
  bool        prev1_before_prev2 = true;

  if (this->m_arr_access.are_on_same_inner_component(hhandle, prev))
  {
    // If prev1 and prev2 are on different components, the insertion of the
    // new curve does not generate a new face, so the way we send these
    // halfedge pointers to the auxiliary function _insert_at_vertices() does
    // not matter.
    // However, in our case, where the two halfedges are reachable from one
    // another and are located on the same hole, a new face will be created
    // and form a hole inside their current incident face. In this case we
    // have to arrange prev1 and prev2 so that the new face (hole) will be
    // incident to the correct halfedge, directed from prev1's target to
    // prev2's target.
    // To do this, we check whether prev1 lies inside the new face we are
    // about to create (or alternatively, whether prev2 does not lie inside
    // this new face).
    const unsigned int  dist1 = this->m_arr_access.halfedge_distance (hhandle, prev);
    const unsigned int  dist2 = this->m_arr_access.halfedge_distance (prev, hhandle);

    prev1_before_prev2 = (dist1 > dist2) ?
      (this->m_arr_access.is_inside_new_face (hhandle, prev, cv.base())) :
      (! this->m_arr_access.is_inside_new_face (prev, hhandle, cv.base()));
  }

  // Perform the insertion.

  new_face_created = false;
  Halfedge_handle  new_he = (prev1_before_prev2) ?
    this->m_arr_access.insert_at_vertices_ex (cv.base(),
                                        hhandle,
                                        prev,
                                        LARGER,
                                        new_face_created,
                                        false) :
    this->m_arr_access.insert_at_vertices_ex (cv.base(),
                                        prev,
                                        hhandle,
                                        SMALLER,
                                        new_face_created,
                                        false);

  if (new_face_created)
  {
    // In case a new face has been created (pointed by the new halfedge we
    // obtained), we have to examine the holes and isolated vertices in the
    // existing face (pointed by the twin halfedge) and move the relevant
    // holes and isolated vertices into the new face.
    this->m_arr_access.relocate_in_new_face (new_he);
  }

  // Return a handle to the new halfedge directed from prev1's target to
  // prev2's target. Note that this may be the twin halfedge of the one
  // returned by _insert_at_vertices();
  if (! prev1_before_prev2)
    new_he = new_he->twin();

  return (new_he);
  }

protected:

  X_monotone_curve_2   sub_cv1;         // Auxiliary variable (for splitting).
  X_monotone_curve_2   sub_cv2;         // Auxiliary variable (for splitting).

};


CGAL_END_NAMESPACE

#endif
