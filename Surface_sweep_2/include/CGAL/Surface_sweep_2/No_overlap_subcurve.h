// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Tali Zvi <talizvi@post.tau.ac.il>,
//             Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Ron Wein <wein@post.tau.ac.il>
//             Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_SURFACE_SWEEP_2_NO_OVERLAP_SUBCURVE_H
#define CGAL_SURFACE_SWEEP_2_NO_OVERLAP_SUBCURVE_H

#include <CGAL/license/Surface_sweep_2.h>

#include <utility>

/*! \file
 *
 * Definition of the No_overlap_subcurve class, which is an
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
template <typename GeometryTraits_2, typename Event_, typename Allocator_, typename Subcurve_>
class No_overlap_subcurve_base {
public:
  using Geometry_traits_2 = GeometryTraits_2;
  using Subcurve = Subcurve_;
  using Event = Event_;
  using Allocator = Allocator_;
  using Handle_overlaps = Tag_false;

private:
  using Gt2 = Geometry_traits_2;

public:
  using Compare_curves = Curve_comparer<Gt2, Event, Subcurve>;
  using Status_line = Multiset<Subcurve*, Compare_curves, Allocator>;
  using Status_line_iterator = typename Status_line::iterator;

protected:
  Status_line_iterator m_hint;      // The location of the subcurve in the
                                    // status line (the Y-structure).
  Event* m_left_event;              // The event associated with the left end.
  Event* m_right_event;             // The event associated with the right end

public:
  /*! checks whether the given event is the matches the left-end event.
   */
  template <typename SweepEvent>
  bool is_start_point(const SweepEvent* event) const
  { return (m_left_event == event); }

  /*! checks whether the given event is the matches the right-end event.
   */
  template <typename SweepEvent>
  bool is_end_point(const SweepEvent* event) const
  { return (m_right_event == event); }

  /*! obtains the event that corresponds to the left end of the subcurve.
   */
  Event* left_event() const { return m_left_event; }

  /*! obtains the event that corresponds to the right end of the subcurve.
   */
  Event* right_event() const { return m_right_event; }

  /*! sets the event that corresponds to the left end of the subcurve. */
  void set_left_event(Event* event) { m_left_event = event; }

  /*! sets the event that corresponds to the right end of the subcurve. */
  void set_right_event(Event* event) { m_right_event = event; }

  /*! obtains the location of the subcurve in the status line .*/
  Status_line_iterator hint() const { return m_hint; }

  /*! sets the location of the subcurve in the status line .*/
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
 * contained in its base) is:
 * - the remaining x-monotone curve that is to the right of the current sweep
 *   line.
 * \tparam GeometryTraits_2 the geometry traits.
 * \tparam Event_ the event type.
 * \tparam Allocator_ a type of an element that is used to acquire/release
 *                    memory for elements of the event queue and the status
 *                    structure, and to construct/destroy the elements in that
 *                    memory. The type must meet the requirements of Allocator.
 * \tparam Subcurve_ the type of the subcurve or Default. If the default is not
 *                   overridden it implies that the type is No_overlap_subcurve.
 */
template <typename GeometryTraits_2, typename Event_,
          typename Allocator_ = CGAL_ALLOCATOR(int),
          typename Subcurve_ = Default>
class No_overlap_subcurve :
  public No_overlap_subcurve_base<GeometryTraits_2, Event_, Allocator_,
                                  typename Default::Get<Subcurve_, No_overlap_subcurve<GeometryTraits_2, Event_,
                                                                                       Allocator_, Subcurve_> >::type> {
public:
  using Geometry_traits_2 = GeometryTraits_2;
  using Event = Event_;
  using Allocator = Allocator_;

private:
  using Gt2 = Geometry_traits_2;
  using Self = No_overlap_subcurve<Gt2, Event, Allocator, Subcurve_>;
  using Subcurve = typename Default::Get<Subcurve_, Self>::type;
  using Base = No_overlap_subcurve_base<Gt2, Event, Allocator, Subcurve>;

public:
  using X_monotone_curve_2 = typename Gt2::X_monotone_curve_2;

  using Status_line_iterator = typename Base::Status_line_iterator;

protected:
  // Data members:
  X_monotone_curve_2 m_last_curve;  // The portion of the curve that lies to
                                    // the right of the last event point
                                    // that occurred on the curve.

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

  /*! Get the last intersecting curve so far (const version). */
  const X_monotone_curve_2& last_curve() const { return m_last_curve; }

  /*! Get the last intersecting curve so far (non-const version). */
  X_monotone_curve_2& last_curve() { return m_last_curve; }

  /*! Set the last intersecting curve so far.
   */
  void set_last_curve(const X_monotone_curve_2& cv) { m_last_curve = cv; }

  /*! Swap the last intersecting curve so far with the given curve. This
   * avoids copying the curve when the given curve is no longer needed, and
   * it hands over the storage of the replaced curve (e.g., the subcurve
   * container of a polycurve) for reuse.
   * \param cv Input: the new last curve; Output: the replaced last curve.
   */
  void swap_last_curve(X_monotone_curve_2& cv) { std::swap(m_last_curve, cv); }

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
