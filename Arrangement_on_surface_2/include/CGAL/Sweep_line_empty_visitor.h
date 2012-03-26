// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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

#ifndef CGAL_SWEEP_LINE_EMPTY_VISITOR_H
#define CGAL_SWEEP_LINE_EMPTY_VISITOR_H

/*! \file
 * Definition of the Sweep_line_empty_visitor class-template.
 */

#include <CGAL/Basic_sweep_line_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>

namespace CGAL {

/*! \class
 * An empty sweep-line visitor that does nothing. It is used as a base-class
 * for other concrete visitors that produce some output.
 */
template <class Traits_,
          class Subcurve_  =  Sweep_line_subcurve<Traits_>,
          typename Event_  =  Sweep_line_event<Traits_, Subcurve_>,
          class Allocator_ =  CGAL_ALLOCATOR(int)>
class Sweep_line_empty_visitor
{
public:

  typedef Traits_                                    Traits_2;
  typedef Subcurve_                                  Subcurve;
  typedef Event_                                     Event;
  typedef Allocator_                                 Allocator;

  typedef typename Traits_2::X_monotone_curve_2      X_monotone_curve_2;
  typedef typename Traits_2::Point_2                 Point_2; 

  typedef Sweep_line_empty_visitor<Traits_2,
                                   Subcurve,
                                   Event,
                                   Allocator>        Self;

private:

  // we want to hide the Sweep_line type
  typedef Basic_sweep_line_2<Traits_2,
                             Self,
                             Subcurve,
                             Event,
                             Allocator>              Sweep_line;
  
  typedef typename Sweep_line::Status_line_iterator  Base_status_line_iter;

public:
  
  /*! \class
   * A wrapper for the base status-line iterator.
   */
  class Status_line_iterator : public Base_status_line_iter
  {
  public:

    typedef Subcurve*                                        value_type;
    typedef value_type&                                      reference;
    typedef value_type*                                      pointer;
    typedef typename Base_status_line_iter::difference_type  difference_type;
    typedef typename Base_status_line_iter::iterator_category
                                                             iterator_category;
    /*! Constructor. */
    Status_line_iterator()
    {}

    /*! Constructor from a base iterator. */
    Status_line_iterator (Base_status_line_iter iter) :
      Base_status_line_iter (iter)
    {}
    
    // Overriden operator*
    reference operator* ()
    {
      return (reinterpret_cast<reference>
              (((Base_status_line_iter*)this)->operator*()));
    }

    // Overriden operator->
    pointer operator-> ()
    {
      return (reinterpret_cast<pointer>
              (Base_status_line_iter::operator->()));
    }
  };


  typedef typename Event::Subcurve_iterator          Event_subcurve_iterator;
  typedef typename Event::Subcurve_reverse_iterator
    Event_subcurve_reverse_iterator;

protected:

  // Data members:
  void*    m_sweep_line;           // The sweep-line object.

public:

  /*! Constructor. */
  Sweep_line_empty_visitor () :
    m_sweep_line (NULL)
  {}

  /*! Attach the a sweep-line object. */
  void attach (void* sl)
  {
    m_sweep_line = sl;
  }

  /*! Destructor */
  virtual ~Sweep_line_empty_visitor()
  {}

  /*!
   * A notification invoked before the sweep-line starts handling the given
   * event.
   */
  void before_handle_event (Event* /* event */)
  {}

  /*!
   * A notification invoked after the sweep-line finishes handling the given
   * event.
   */
  bool after_handle_event (Event* /* event */,
                           Status_line_iterator /* iter */,
                           bool /* flag */)
  {
    return (true);
  }

  /*! A notification invoked when a new subcurve is created. */
  void add_subcurve (X_monotone_curve_2 /* cv */,
                     Subcurve* /* sc */)
  {}

  /*! A notification issued before the sweep process starts. */
  void before_sweep()
  {}

  /*! A notification issued after the sweep process ends. */
  void after_sweep()
  {}

  /*! Update the event to be the given curve end. */
  void update_event (Event* /* e */,
                     const Point_2& /* end_point */,
                     const X_monotone_curve_2& /* cv */,
                     Arr_curve_end /* cv_end */,
                     bool /* is_new */)
  {}

  /*! Update the event to be the given infinite curve end. */
  void update_event (Event* /* e */,
                     const X_monotone_curve_2& /* cv */,
                     Arr_curve_end /* cv_end */,
                     bool /* is_new */)
  {}

  /*! Update the event to be the intersection point of two subcurves. */
  void update_event (Event* /* e */,
                     Subcurve* /* sc1 */,
                     Subcurve* /* sc2 */,
                     bool /* is_new */)
  {}

  /*! Update the event. */
  void update_event (Event* /* e */,
                     Subcurve* /* sc1 */)
  {}

  /*! Update the event. */
  void update_event (Event* /* e */,
                     const Point_2& /* pt */,
                     bool /* is_new */)
  {}

  /* Found overlap */
  void found_overlap (Subcurve* /* sc1 */,
                      Subcurve* /* sc2 */,
                      Subcurve* /* ov_sc */)
  {}

  /*! Get the first subcurve in the status line. */
  Status_line_iterator status_line_begin()
  {
    return (this->_sweep_line()->status_line_begin());
  }

  /*! Get a past-the-end iterator for the subcurves in the status line. */
  Status_line_iterator status_line_end()
  {
    return (this->_sweep_line()->status_line_end());
  }

  /*! Get the position of the given subcurve in the status line. */
  Status_line_iterator status_line_position (Subcurve *sc)
  {
    return (sc->hint());
  }

  /*! Get the number of subcurves in the status line. */
  unsigned status_line_size() const
  {
    return (this->_sweep_line()->status_line_size());
  }

  /*! Check if the status line is empty. */
  bool is_status_line_empty() const
  {
    return (this->_sweep_line()->is_status_line_empty());
  }

  /*! Deallocate the given event. */
  void deallocate_event (Event* e)
  {
    this->_sweep_line()->deallocate_event(e);
  }

  /*! Stop the sweep-line process. */
  void stop_sweep()
  {
    this->_sweep_line()->stop_sweep();
  }

  /*! Get the sweep-line object. */
  void* sweep_line ()
  {
    return (m_sweep_line);
  }

  /*! Get the current event. */
  Event* current_event()
  {
    return (this->_sweep_line()->current_event());
  }

  /*! Get the geometry-traits class. */
  const Traits_2 * traits ()
  {
    return (this->_sweep_line()->traits());
  }

private:

  /*! Get the sweep-line object. */
  Sweep_line* _sweep_line ()
  {
    return (reinterpret_cast<Sweep_line*>(m_sweep_line));
  }

  const Sweep_line* _sweep_line () const
  {
    return (reinterpret_cast<const Sweep_line*>(m_sweep_line));
  }
};

} //namespace CGAL

#endif
