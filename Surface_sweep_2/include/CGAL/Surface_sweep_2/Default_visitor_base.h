// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_SURFACE_SWEEP_2_DEFAULT_VISITOR_BASE_H
#define CGAL_SURFACE_SWEEP_2_DEFAULT_VISITOR_BASE_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Definition of the Default_visitor_base class-template.
 */

#include <CGAL/No_intersection_surface_sweep_2.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class Default_visitor_base
 *
 * An empty surface-sweep visitor that does little. It is used as a base-class
 * for other concrete visitors that produce some output.
 *
 * The surface-sweep framework consists of hierarchy of several types, where one
 * derives from the other. The base type is capable of sweeping curves that do
 * not intersect in their interiors. It is implemented as a template that is
 * parameterized by a visitor. The user is allowed to introduce new visitor
 * types that derive from the Default_visitor_base basic type. However, this
 * basic type provides access to the base type of the surface-sweep
 * hierarchy. We use the curiously recurring template pattern (CRTP) idiom to
 * have access to the correct base surface-sweep type.
 * \tparam GeometryTraits_2 the geometry traits trype.
 * \tparam Event_ the sweep event type.
 * \tparam Subcurve_ the sweep subcurve type.
 * \tparam Allocator_ the allocator used to allocate events and subcurves during
 *                    the sweep.
 * \tparam Visitor_ the actual visitor used, which is derived from
 *                  Default_visitor_base.
 */
template <typename GeometryTraits_2, typename Event_, typename Subcurve_,
          typename Allocator_, typename Visitor_>
class Default_visitor_base {
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Event_                                        Event;
  typedef Subcurve_                                     Subcurve;
  typedef Allocator_                                    Allocator;
  typedef Visitor_                                      Visitor;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef Default_visitor_base<Gt2, Event, Subcurve, Allocator, Visitor>
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
  Default_visitor_base () : m_surface_sweep(nullptr) {}

  /*! Destructor */
  virtual ~Default_visitor_base() {}

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

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
