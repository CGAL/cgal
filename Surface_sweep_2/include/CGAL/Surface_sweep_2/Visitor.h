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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_SURFACE_SWEEP_2_VISITOR_H
#define CGAL_SURFACE_SWEEP_2_VISITOR_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Definition of the Surface_sweep_empty_visitor class-template.
 */

#include <CGAL/No_intersection_surface_sweep_2.h>
#include <CGAL/Surface_sweep_2/Default_event.h>
#include <CGAL/Surface_sweep_2/Default_subcurve.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class Visitor_base
 *
 * An empty surface-sweep visitor that does nothing. It is used as a base-class
 * for other concrete visitors that produce some output.
 */
template <typename GeometryTraits_2, typename Event_, typename Subcurve_,
          typename Allocator_, typename Visitor_>
class Visitor_base {
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Event_                                        Event;
  typedef Subcurve_                                     Subcurve;
  typedef Allocator_                                    Allocator;
  typedef Visitor_                                      Visitor;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef Visitor_base<Gt2, Event, Subcurve, Allocator, Visitor>
                                                        Self;

public:
  typedef typename Subcurve::Status_line_iterator       Status_line_iterator;

  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Gt2::Point_2                         Point_2;

  typedef typename Event::Subcurve_iterator             Event_subcurve_iterator;
  typedef typename Event::Subcurve_reverse_iterator
    Event_subcurve_reverse_iterator;

  typedef No_intersection_surface_sweep_2<Visitor>      Surface_sweep_2;

protected:
  // Data members:
  Surface_sweep_2* m_surface_sweep;           // The sweep-line object.

public:
  /*! Constructor. */
  Visitor_base () : m_surface_sweep(NULL) {}

  /*! Destructor */
  virtual ~Visitor_base() {}

  /*! Attach the a sweep-line object. */
  void attach(Surface_sweep_2* sl) { m_surface_sweep = sl; }

  /*!
   * A notification invoked before the sweep-line starts handling the given
   * event.
   */
  void before_handle_event(Event* /* event */) {}

  /*!
   * A notification invoked after the sweep-line finishes handling the given
   * event.
   */
  bool after_handle_event(Event* /* event */,
                          Status_line_iterator /* iter */,
                          bool /* flag */)
  { return true; }

  /*! A notification invoked when a new subcurve is created. */
  void add_subcurve(X_monotone_curve_2 /* cv */,
                    Subcurve* /* sc */)
  {}

  /*! A notification issued before the sweep process starts. */
  void before_sweep()
  {}

  /*! A notification issued after the sweep process ends. */
  void after_sweep()
  {}

  /*! Update the event to be the given curve end. */
  void update_event(Event* /* e */,
                    const Point_2& /* end_point */,
                    const X_monotone_curve_2& /* cv */,
                    Arr_curve_end /* cv_end */,
                    bool /* is_new */)
  {}

  /*! Update the event to be the given infinite curve end. */
  void update_event(Event* /* e */,
                    const X_monotone_curve_2& /* cv */,
                    Arr_curve_end /* cv_end */,
                    bool /* is_new */)
  {}

  /*! Update the event to be the intersection point of two subcurves. */
  void update_event(Event* /* e */,
                    Subcurve* /* sc1 */,
                    Subcurve* /* sc2 */,
                    bool /* is_new */)
  {}

  /*! Update the event. */
  void update_event(Event* /* e */,
                    Subcurve* /* sc1 */)
  {}

  /*! Update the event. */
  void update_event(Event* /* e */,
                    const Point_2& /* pt */,
                    bool /* is_new */)
  {}

  /* Found overlap */
  void found_overlap(Subcurve* /* sc1 */,
                     Subcurve* /* sc2 */,
                     Subcurve* /* ov_sc */)
  {}

  /*! Obtain the first subcurve in the status line. */
  Status_line_iterator status_line_begin()
  { return surface_sweep()->status_line_begin(); }

  /*! Obtain a past-the-end iterator for the subcurves in the status line. */
  Status_line_iterator status_line_end()
  { return surface_sweep()->status_line_end(); }

  /*! Obtain the number of subcurves in the status line. */
  unsigned status_line_size() const
  { return surface_sweep()->status_line_size(); }

  /*! Check if the status line is empty. */
  bool is_status_line_empty() const
  { return surface_sweep()->is_status_line_empty(); }

  /*! Deallocate the given event. */
  void deallocate_event(Event* e) { surface_sweep()->deallocate_event(e); }

  /*! Stop the sweep-line process. */
  void stop_sweep() { surface_sweep()->stop_sweep(); }

  /*! Obtain the sweep-line object. */
  Surface_sweep_2* surface_sweep() { return m_surface_sweep; }

  /*! Obtain the sweep-line object. */
  const Surface_sweep_2* surface_sweep() const { return m_surface_sweep; }

  /*! Obtain the current event. */
  Event* current_event() { return surface_sweep()->current_event(); }

  /*! Obtain the geometry-traits class. */
  const Gt2* traits() { return surface_sweep()->traits(); }
};

template <typename Visitor_,
          typename GeometryTraits_2,
          typename Event_ = Default_event<GeometryTraits_2>,
          typename Subcurve_ = Default_subcurve<GeometryTraits_2, Event_>,
          typename Allocator_ = CGAL_ALLOCATOR(int)>
class Default_visitor : public Visitor_base<GeometryTraits_2, Event_, Subcurve_,
                                            Allocator_, Visitor_>
{
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Event_                                        Event;
  typedef Subcurve_                                     Subcurve;
  typedef Allocator_                                    Allocator;
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
