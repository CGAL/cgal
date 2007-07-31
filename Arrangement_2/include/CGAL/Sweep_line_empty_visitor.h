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
// $URL: svn+ssh://golubevs@scm.gforge.inria.fr/svn/cgal/trunk/Arrangement_2/include/CGAL/Sweep_line_2_empty_visitor.h $
// $Id: Sweep_line_2_empty_visitor.h 37160 2007-03-16 14:33:01Z lrineau $
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>


#ifndef CGAL_SWEEP_LINE_EMPTY_VISITOR_H
#define CGAL_SWEEP_LINE_EMPTY_VISITOR_H

#include <CGAL/Basic_sweep_line_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>

CGAL_BEGIN_NAMESPACE

template <class Traits_,
          class Subcurve_  =  Sweep_line_subcurve<Traits_>,
          typename Event_  =  Sweep_line_event<Traits_, Subcurve_>,
          class Allocator_ =  CGAL_ALLOCATOR(int)>
class Sweep_line_empty_visitor
{
public:

  typedef Event_        Event;
  typedef Traits_       Traits;
  typedef Subcurve_     Subcurve;
  typedef Allocator_    Allocator;

  typedef typename Traits::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Traits::Point_2                         Point_2; 

  typedef Sweep_line_empty_visitor<Traits,
                        Subcurve,
                        Event,
                        Allocator>                         Self;

private:

  // we want to hide the Sweep_line type
  typedef Basic_sweep_line_2<Traits,
                             Self,
                             Subcurve,
                             Event,
                             Allocator>                    Sweep_line;
  


  typedef typename Sweep_line::StatusLineIter              StatusLineIter;

  public:
  class Status_line_iterator : public StatusLineIter
  {
  public:

    typedef Subcurve*                     value_type;
    typedef  value_type&                  reference;
    typedef  value_type*                  pointer;
    //typedef typename StatusLineIter::distance_type      distance_type;
    //typedef typename StatusLineIter::iterator_category  iterator_category;
    
    Status_line_iterator() {}

    Status_line_iterator(StatusLineIter iter) : StatusLineIter(iter)
    {}
    
    //override operator*
    reference operator* ()
    {
      return (reinterpret_cast<reference>(((StatusLineIter*)this)->operator*()));
    }

    //override operator->
    pointer operator-> ()
    {
      return (reinterpret_cast<pointer>(StatusLineIter::operator->()));
    }

  };


  typedef typename Event::SubCurveIter                     SubCurveIter;
  typedef typename Event::SubCurveRevIter                  SubCurveRevIter;

 

protected:

  void*    m_sweep_line;


private:

  Sweep_line* sweep_line()
  {
    return ( reinterpret_cast<Sweep_line*>(m_sweep_line) );
  }



public:

  void attach(void* sl)
  {
    m_sweep_line = sl;
  }

  Sweep_line_empty_visitor(){}

  /*! Destructor */
  virtual ~Sweep_line_empty_visitor() {}

  void before_handle_event(Event*){}

  bool after_handle_event(Event*, Status_line_iterator, bool)
  {
    return true;
  }

  void add_subcurve(X_monotone_curve_2, Subcurve*){}

  
  void init_event(Event* ){}

  void before_sweep(){}
  void after_sweep(){}
  void after_init(){}

  //// ioslated point falls on existing event
  //void update_event(Event* e, const Point_2& pt)
  //{}

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

  void update_event(Event*, const Point_2&)
  {}


  Status_line_iterator status_line_begin()
  {
    return (sweep_line()->status_line_begin());
  }

  Status_line_iterator status_line_end()
  {
    return (sweep_line()->status_line_end());
  }

  Status_line_iterator status_line_position(Subcurve *sc)
  {
    return (sc->get_hint());
  }

  unsigned status_line_size() const
  {
    return (sweep_line()->status_line_size());
  }

  bool is_status_line_empty() const
  {
    return (sweep_line()->is_status_line_empty());
  }

  void deallocate_event(Event* e)
  {
    sweep_line()->deallocate_event(e);
  }

  void stop_sweep()
  {
    sweep_line()->stop_sweep();
  }

  Event* current_event()
  {
    return sweep_line()-> current_event();
  }

  Traits* traits()
  {
    return sweep_line()->traits();
  }


};

CGAL_END_NAMESPACE

#endif
