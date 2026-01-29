// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Baruch Zukerman  <baruchzu@post.tau.ac.il>
//             Efi Fogel        <efif@post.tau.ac.il>

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
 * An surface-sweep visitor used as a base-class for other concrete visitors
 * that produce some output.
 *
 * The surface-sweep framework consists of hierarchy of several types. The base
 * type, namely `No_intersection_surface_sweep_2`, is capable of sweeping curves
 * that do not intersect in their interiors.  It is implemented as a template
 * that is parameterized by a visitor. Visitors are also organized in a
 * hierarchy. This is the base type of the visitor hierarchy; it stores a
 * pointer to the actual surface-sweep object. However, the stored object type
 * is hard coded as the base type of the surface-sweep hierarchy, namely
 * `No_intersection_surface_sweep_2`. Therefore, it only provides functionality
 * available through virtual inheritance.
 *
 * We use the curiously recurring template pattern (CRTP) idiom to gain access
 * to the actual visitor type (the `Visitor_` template parameter) and in turn to
 * the correct base surface-sweep type.
 *
 * \tparam GeometryTraits_2 the* geometry traits trype.
 * \tparam Event_ the sweep event type.
 * \tparam Subcurve_ the sweep subcurve type.
 * \tparam Allocator_ the allocator used to allocate events and subcurves during the sweep.
 * \tparam Visitor_ the actual visitor used, which is derived from Default_visitor_base.
 *
 * For future reference, we can extend the use of (CRTP) to gain access to the
 * correct type of the actual surface-sweep object (instead of the base type).
 * Once this is in place, the visitor can call any function of the surface-sweep
 * object (and not only virtual functions). The application should be carried
 * out as follows:
 *
 * 1. Apply CRTP to the type in the surface-sweep hierarchy, for example:
 *
 *   template <typename Visitor_, typename SurfaceSweep_2 = void>
 *   class No_intersection_surface_sweep_2 {
 *   private:
 *     using Visitor = Visitor_;
 *     using Self = No_intersection_surface_sweep_2<Visitor, SurfaceSweep_2>;
 *     using Surface_sweep_2 = typename std::conditional<std::is_void<SurfaceSweep_2>::value, Self, SurfaceSweep>::type;
 *     ...
 *   };
 *
 * 2. Start a new hierarchy of visitors, and add yet another template parameter
 * to the base template class, called `SurfaceSweep` that is the type of the
 * actual surface-sweep object. (It's possible that adding it, makes the
 * `Visitor_` parameter expendable.) For example:
 *
 *   template <typename GeometryTraits_2, typename Event_, typename Subcurve_, typename Allocator_,
 *             typename SurfaceSweep_2>
 *   class Surface_sweep_visitor_base {
 *     using Surface_sweep_2 = SurfaceSweep_2;
 *     using Visitor = typename Surface_sweep_2::Visitor;
 *   };
 *
 * 3. Populate the new hierarchy, and gradually replace the visitors in the old
 * hierarchy with new visitors that derive from `Surface_sweep_visitor_base`. In
 * this new setup, a concrete visitor defines the surface-sweep type in use, and
 * derives from `Surface_sweep_visitor_base`, or better yet, from a new utility
 * class template called `Default_surface_sweep_visitor`:
 *
 *   template <typename SurfaceSweep_2,
 *             typename GeometryTraits_2,
 *             typename Allocator_ = CGAL_ALLOCATOR(int),
 *             typename Event_ = Default_event<GeometryTraits_2, Allocator_>,
 *             typename Subcurve_ = Default_subcurve<GeometryTraits_2, Event_, Allocator_>>
 *   class Default_surface_sweep_visitor :
 *     public Surface_sweep_visitor_base<GeometryTraits_2, Event_, Subcurve_, Allocator_, SurfaceSweep_2> {
 *     using Surface_sweep_2 = SurfaceSweep_2;
 *     using Geometry_traits_2 = GeometryTraits_2;
 *     using Allocator = Allocator_;
 *     using Event = Event_;
 *     using Subcurve = Subcurve_;
 *   };
 *
 * and
 *
 *   template <typename GeometryTraits_2, typename Allocator_ = CGAL_ALLOCATOR(int)>
 *   class Do_intersect_visitor :
 *     public Default_surface_sweep_visitor<
 *       Do_intersect_surface_sweep_2<Do_intersect_visitor<GeometryTraits_2, Allocator_>>,
 *       GeometryTraits_2, Allocator_> {
 *   };
 */
template <typename GeometryTraits_2, typename Event_, typename Subcurve_, typename Allocator_, typename Visitor_>
class Default_visitor_base {
public:
  using Geometry_traits_2 = GeometryTraits_2;
  using Event = Event_;
  using Subcurve = Subcurve_;
  using Allocator = Allocator_;
  using Visitor = Visitor_;

private:
  using Gt2 = Geometry_traits_2;
  using Self = Default_visitor_base<Gt2, Event, Subcurve, Allocator, Visitor>;

public:
  using Status_line_iterator = typename Subcurve::Status_line_iterator;

  using X_monotone_curve_2 = typename Gt2::X_monotone_curve_2;
  using Point_2 = typename Gt2::Point_2;
  using Multiplicity = typename Gt2::Multiplicity;

  using Event_subcurve_iterator = typename Event::Subcurve_iterator;
  using Event_subcurve_reverse_iterator = typename Event::Subcurve_reverse_iterator;

  using Surface_sweep_2 = No_intersection_surface_sweep_2<Visitor>;

protected:
  // Data members
  Surface_sweep_2* m_surface_sweep;           // the sweep-line object.

public:
  /*! constructs. */
  Default_visitor_base () : m_surface_sweep(nullptr) {}

  /*! destructs */
  virtual ~Default_visitor_base() {}

  /*! attaches a sweep-line object.
   */
  void attach(Surface_sweep_2* sl) { m_surface_sweep = sl; }

  /*! A notification invoked before the sweep-line starts handling the given event. */
  void before_handle_event(Event* /* event */) {}

  /*! A notification invoked after the sweep-line finishes handling the given event. */
  bool after_handle_event(Event* /* event */, Status_line_iterator /* iter */, bool /* flag */) { return true; }

  /*! A notification invoked when a new subcurve is created. */
  void add_subcurve(X_monotone_curve_2 /* cv */, Subcurve* /* sc */) {}

  /*! A notification issued before the sweep process starts. */
  void before_sweep() {}

  /*! A notification issued after the sweep process ends. */
  void after_sweep() {}

  /*! updates the event to be the given curve end. */
  void update_event(Event* /* e */, const Point_2& /* end_point */, const X_monotone_curve_2& /* cv */,
                    Arr_curve_end /* cv_end */, bool /* is_new */)
  {}

  /*! updates the event to be the given infinite curve end. */
  void update_event(Event* /* e */, const X_monotone_curve_2& /* cv */, Arr_curve_end /* cv_end */,
                    bool /* is_new */)
  {}

  /*! updates the event to be the intersection point of two subcurves. */
  void update_event(Event* /* e */,
                    Subcurve* /* sc1 */,
                    Subcurve* /* sc2 */,
                    bool /* is_new */,
                    Multiplicity /* multiplicity */)
  {}

  /*! updates the event. */
  void update_event(Event* /* e */, Subcurve* /* sc1 */) {}

  /*! updates the event. */
  void update_event(Event* /* e */, const Point_2& /* pt */, bool /* is_new */) {}

  /*! records that an overlap has been found */
  void found_overlap(Subcurve* /* sc1 */, Subcurve* /* sc2 */, Subcurve* /* ov_sc */) {}

  /*! obtains the first subcurve in the status line. */
  Status_line_iterator status_line_begin() { return m_surface_sweep->status_line_begin(); }

  /*! obtains a past-the-end iterator for the subcurves in the status line. */
  Status_line_iterator status_line_end() { return m_surface_sweep->status_line_end(); }

  /*! obtains the number of subcurves in the status line. */
  unsigned status_line_size() const { return m_surface_sweep->status_line_size(); }

  /*! checkes whether the status line is empty. */
  bool is_status_line_empty() const { return m_surface_sweep->is_status_line_empty(); }

  /*! Deallocate the given event. */
  void deallocate_event(Event* e) { m_surface_sweep->deallocate_event(e); }

  /*! Stop the sweep-line process. */
  void stop_sweep() { m_surface_sweep->stop_sweep(); }

  /*! obtains the current event. */
  Event* current_event() { return m_surface_sweep->current_event(); }

  /*! obtains the geometry-traits class. */
  const Gt2* traits() { return m_surface_sweep->traits(); }
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
