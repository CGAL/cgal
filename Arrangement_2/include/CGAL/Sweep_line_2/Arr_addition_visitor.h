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

#ifndef CGAL_ARR_ADDITION_VISITOR_H
#define CGAL_ARR_ADDITION_VISITOR_H

#include <CGAL/Sweep_line_2/Arr_basic_addition_visitor.h>

CGAL_BEGIN_NAMESPACE

template<class Traits, class Arrangement_, class Event,class Subcurve>
class Arr_addition_visitor : 
  public Arr_basic_addition_visitor<Traits,Arrangement_,Event,Subcurve>
{
protected:
  typedef Arr_basic_addition_visitor<Traits,
                                     Arrangement_,
                                     Event,
                                     Subcurve>    Base;
  typedef Arrangement_                            Arrangement;
  typedef typename Base::Halfedge_handle          Halfedge_handle;
  typedef typename Base::Point_2                  Point_2;

public:

  Arr_addition_visitor(Arrangement* arr) : Base(arr)
  {}

  /*void before_sweep()
  {}*/

  // check if the halfedge associated with 'sc' will be splitted at the given
  // event point (its a recursive function since the hirearchy of potential 
  // ovrlap Subcuves
  virtual bool is_split_event(Subcurve* sc, Event* event)
  {
    if(sc->get_last_curve().get_halfedge_handle() == Halfedge_handle(NULL))
      return (false);

    if(! sc->get_orig_subcurve1())
    {
      return (reinterpret_cast<Event*>(sc->get_left_event())!= 
              this->current_event());
    }
    return
      (this->is_split_event(reinterpret_cast<Subcurve*>(sc->get_orig_subcurve1()),
                      event) || 
       this->is_split_event(reinterpret_cast<Subcurve*>(sc->get_orig_subcurve2()),
                      event));
  }

  virtual Halfedge_handle split_edge(Halfedge_handle he,
                                     Subcurve* sc,
                                     const Point_2& pt)
  {
    // make sure that the halfedge associated with 'sc' is the directed from
    // right to left , since we always 'look' above , and the incident face 
    //is on the left of the  halfedge
    CGAL_assertion (he->direction() == LARGER);
    
    this->traits()->split_2_object()(he->curve(), pt, this->sub_cv2, this->sub_cv1);
    Halfedge_handle new_he =  
      this->m_arr_access.split_edge_ex (he,
                                        pt.base(),
                                        this->sub_cv1.base(),
                                        this->sub_cv2.base());

    Event* last_event_on_sc = reinterpret_cast<Event*>(sc->last_event());
    if(last_event_on_sc->get_halfedge_handle() == he)
      last_event_on_sc->set_halfedge_handle(new_he->next());
    return (new_he);
  }
     

   
};
CGAL_END_NAMESPACE

#endif
