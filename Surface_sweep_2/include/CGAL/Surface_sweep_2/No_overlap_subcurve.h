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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Tali Zvi <talizvi@post.tau.ac.il>,
//             Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Ron Wein <wein@post.tau.ac.il>
//             Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_SURFACE_SWEEP_2_NO_OVERLAP_SUBCURVE_H
#define CGAL_SURFACE_SWEEP_2_NO_OVERLAP_SUBCURVE_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Defintion of the No_overlap_subcurve class, which is an
 * extended curve type, referred to as Subcurve, used by the surface-sweep
 * framework.
 *
 * The surface-sweep framework is implemented as a template that is
 * parameterized by a visitor; the visitor is parameterized, among the other, by
 * the Subcurve and Event types. That is, instance types of Subcurve and Event
 * must be available when the surface-sweep template is instantiated.
 *
 * No_overlap_subcurve is the most basic type. The user is allowed
 * to introduce new types that derive from the basic type. However, some of the
 * fields of the basic type depends on the Subcurve type.  We use the curiously
 * recurring template pattern (CRTP) idiom to force the correct matching of
 * these types.
 */

#include <CGAL/Surface_sweep_2/Curve_comparer.h>
#include <CGAL/Multiset.h>
#include <CGAL/assertions.h>
#include <CGAL/Default.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class No_overlap_subcurve_base
 *
 * This is the base class of the No_overlap_subcurve class
 * template used by the (CRTP) idiom.
 * \tparam GeometryTraits_2 the geometry traits.
 * \tparam Subcurve_ the subcurve actual type.
 *
 * The information contained in this class is:
 * - two event points, which are associated with the left and right end of the
 *   curve.
 * - an iterator that points to the location of the subcurve in the status line.
 */
template <typename GeometryTraits_2, typename Event_, typename Allocator_,
          typename Subcurve_>
class No_overlap_subcurve_base {
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Subcurve_                                     Subcurve;
  typedef Event_                                        Event;
  typedef Allocator_                                    Allocator;

private:
  typedef Geometry_traits_2                             Gt2;

public:
  typedef Curve_comparer<Gt2, Event, Subcurve>          Compare_curves;
  typedef Multiset<Subcurve*, Compare_curves, Allocator>
                                                        Status_line;
  typedef typename Status_line::iterator                Status_line_iterator;

protected:
  Status_line_iterator m_hint;      // The location of the subcurve in the
                                    // status line (the Y-structure).
  Event* m_left_event;              // The event associated with the left end.
  Event* m_right_event;             // The event associated with the right end

public:
  /*! Check whether the given event is the matches the left-end event.
   */
  template <typename SweepEvent>
  bool is_start_point(const SweepEvent* event) const
  { return (m_left_event == event); }

  /*! Check whether the given event is the matches the right-end event.
   */
  template <typename SweepEvent>
  bool is_end_point(const SweepEvent* event) const
  { return (m_right_event == event); }

  /*! Obtain the event that corresponds to the left end of the subcurve.
   */
  Event* left_event() const { return m_left_event; }

  /*! Obtain the event that corresponds to the right end of the subcurve.
   */
  Event* right_event() const { return m_right_event; }

  /*! Set the event that corresponds to the left end of the subcurve. */
  void set_left_event(Event* event) { m_left_event = event; }

  /*! Set the event that corresponds to the right end of the subcurve. */
  void set_right_event(Event* event) { m_right_event = event; }

  /*! Obtain the location of the subcurve in the status line .*/
  Status_line_iterator hint() const { return m_hint; }

  /*! Set the location of the subcurve in the status line .*/
  void set_hint(Status_line_iterator hint) { m_hint = hint; }
};

/*! \class No_overlap_subcurve
 *
 * This is a class template that wraps a traits curve of type
 * X_monotone_curve_2.  It contains data that is used when applying the sweep
 * algorithm on a set of x-monotone curves. This class derives from the
 * No_overlap_subcurve_base class template.
 *
 * The information contained in this class (in addition to the information
 * contaisn in its base) is:
 * - the remaining x-monotone curve that is to the right of the current sweep
 *   line.
 * \tparam GeometryTraits_2 the geometry traits.
 * \tparam Event_ the event type.
 * \tparam Allocator_ a type of an element that is used to acquire/release
 *                    memory for elements of the event queue and the status
 *                    structure, and to construct/destroy the elements in that
 *                    memory. The type must meet the requirements of Allocator.
 * \tparam Subcurve_ the type of the subcurve or Default. If the default is not
 *                   overriden it implies that the type is No_overlap_subcurve.
 */
template <typename GeometryTraits_2, typename Event_,
          typename Allocator_ = CGAL_ALLOCATOR(int),
          typename Subcurve_ = Default>
class No_overlap_subcurve :
  public No_overlap_subcurve_base<
    GeometryTraits_2, Event_, Allocator_,
    typename Default::Get<Subcurve_,
                          No_overlap_subcurve<GeometryTraits_2, Event_,
                                              Allocator_, Subcurve_> >::type>
{
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Event_                                        Event;
  typedef Allocator_                                    Allocator;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef No_overlap_subcurve<Gt2, Event, Allocator, Subcurve_>
                                                        Self;
  typedef typename Default::Get<Subcurve_, Self>::type  Subcurve;
  typedef No_overlap_subcurve_base<Gt2, Event, Allocator, Subcurve>
                                                        Base;

public:
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;

  typedef typename Base::Status_line_iterator           Status_line_iterator;

protected:
  // Data members:
  X_monotone_curve_2 m_last_curve;  // The portion of the curve that lies to
                                    // the right of the last event point
                                    // that occured on the curve.

public:
  /*! Construct default. */
  No_overlap_subcurve() {}

  /*! Construct from a curve.
   * \param curve the input x-monotone curve.
   */
  No_overlap_subcurve(const X_monotone_curve_2& curve) :
    m_last_curve(curve)
  {}

  /*! Initialize the subcurves by setting the curve. */
  void init(const X_monotone_curve_2& curve) { m_last_curve = curve; }

  /*! Destruct. */
  ~No_overlap_subcurve() {}

  /*! Get the last intersecing curve so far (const version). */
  const X_monotone_curve_2& last_curve() const { return m_last_curve; }

  /*! Get the last intersecing curve so far (non-const version). */
  X_monotone_curve_2& last_curve() { return m_last_curve; }

  /*! Set the last intersecing curve so far.
   */
  void set_last_curve(const X_monotone_curve_2& cv) { m_last_curve = cv; }

#ifdef CGAL_SS_VERBOSE
  void Print() const;
#endif
};

#ifdef CGAL_SS_VERBOSE
template <typename Gt2, typename Evt, typename Allocator, typename Scv>
void No_overlap_subcurve<Gt2, Evt, Allocator, Scv>::Print() const
{ std::cout << "Curve " << this << "  (" << last_curve() << ") "; }
#endif

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
