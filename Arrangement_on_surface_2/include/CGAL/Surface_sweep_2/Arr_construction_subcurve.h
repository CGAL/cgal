// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARR_CONSTRUCTION_SUBCURVE_H
#define CGAL_ARR_CONSTRUCTION_SUBCURVE_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_construction_subcurve class-template, which is an
 * extended curve type, referred to as Subcurve, used by the surface-sweep
 * framework.
 *
 * The surface-sweep framework is implemented as a template that is
 * parameterized, among the other, by the Subcurve and Event types. That is,
 * instance types of Subcurve and Event must be available when the
 * surface-sweep template is instantiated.
 *
 * Arr_construction_subcurve derives from an instance of the Default_subcurve
 * class template. The user is allowed to introduce new types that derive from
 * an instance of the Arr_construction_subcurve class template. However, some of
 * the fields of this template depends on the Subcurve type.  We use the
 * curiously recurring template pattern (CRTP) idiom to force the correct
 * matching of these types.
 */

#include <CGAL/Surface_sweep_2/Default_subcurve.h>
#include <CGAL/Default.h>

namespace CGAL {

namespace Ss2 = Surface_sweep_2;

/*! \class Arr_construction_subcurve_base
 *
 * This is the base class of the Arr_construction_subcurve class template used
 * by the (CRTP) idiom.
 * \tparam GeometryTraits_2 the geometry traits.
 * \tparam Event_ the event type.
 * \tparam Allocator_ a type of an element that is used to acquire/release
 *                    memory for elements of the event queue and the status
 *                    structure, and to construct/destroy the elements in that
 *                    memory. The type must meet the requirements of Allocator.
 * \tparam Subcurve_ the subcurve actual type.
 *
 * The information contained in this class last:
 * - ishe  event that was handled on the curve.
 * - The index for a subcurve that may represent a hole
 * - Indices of all halfedge below the curve that may represent a hole.
 */
template <typename GeometryTraits_2, typename Event_, typename Allocator_,
          template <typename, typename, typename, typename>
          class SurfaceSweepBaseCurve,
          typename Subcurve_>
class Arr_construction_subcurve_base :
    public SurfaceSweepBaseCurve<GeometryTraits_2, Event_, Allocator_, Subcurve_>
{
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Subcurve_                                     Subcurve;
  typedef Event_                                        Event;
  typedef Allocator_                                    Allocator;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef SurfaceSweepBaseCurve<Gt2, Event, Allocator, Subcurve>
                                                        Base;

public:
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef Event*                                        Event_ptr;
  typedef std::list<unsigned int>                       Halfedge_indices_list;

  /*! Construct default. */
  Arr_construction_subcurve_base() :
    Base(),
    m_last_event(0),
    m_index(0)
  {}

  /*! Constructor from an x-monotone curve. */
  Arr_construction_subcurve_base(X_monotone_curve_2& curve) :
    Base(curve),
    m_last_event(0),
    m_index(0)
  {}

protected:
  // Data members:
  Event_ptr m_last_event;    // The last event that was handled on the curve.

  /*! index for a subcurve that may represent a hole (emarge from the left
   * most vertex of a hole, and its the upper most curve). other subcurves
   * will have 0 value  (invalid index)
   */
  unsigned int m_index;      // Index for a subcurve that may represent a hole
                             // (emarge from the leftmost vertex of a hole,
                             // and it is the topmost curve). Other subcurves
                             // have a 0 (invalid) index.

  Halfedge_indices_list m_halfedge_indices;
                             // Indices of all halfedge below the curve that
                             // may represent a hole.

public:
  /*! Initialize the curve. */
  void init(const X_monotone_curve_2& curve) { Base::init(curve); }

  /*! Set the event associated with the left end of the subcurve. */
  template <typename SweepEvent>
  void set_left_event(SweepEvent* left)
  {
    Base::set_left_event(left);
    set_last_event(left);
  }

  /*! Set the last event on the subcurve. */
  void set_last_event(Event_ptr e) { m_last_event = e; }

  /*! Obtain the last event. */
  Event_ptr last_event() const { return m_last_event; }

  /*! Obtain the subcurve index. */
  unsigned int index() const { return m_index; }

  /*! Set the subcurve index. */
  void set_index(unsigned int i) { m_index = i; }

  /*! Check whether the index is valid. */
  bool has_valid_index() const { return (m_index != 0); }

  /*! Add an index of a halfedge below the subcurve. */
  void add_halfedge_index(unsigned int i) { m_halfedge_indices.push_back(i); }

  /*! Clear the indices of the halfedges below the subcurve. */
  void clear_halfedge_indices() { m_halfedge_indices.clear(); }

  /*! Check if there are any halfedges below the subcurve. */
  bool has_halfedge_indices() const { return (!m_halfedge_indices.empty()); }

  /*! Obtain the indices of the halfedges below the subcurve. */
  Halfedge_indices_list& halfedge_indices_list() { return m_halfedge_indices; }
};

/*! \class Arr_construction_subcurve
 *
 * This class that holds information about a curve that is added to the
 * arrangement.  In addition to the information that is contained in
 * Surface_sweep_subcurve, when an arrangement is constructed, a pointer to the
 * last handled event on the curve is stored (in the base class). This
 * information is used to retrieve hints when a subcurve of this curve is
 * inserted into the planar map.
 *
 * Inherits from `Surface_sweep_subcurve`
 * \sa `Surface_sweep_subcurve`
 */
template <typename GeometryTraits_2, typename Event_,
          typename Allocator_ = CGAL_ALLOCATOR(int),
          template <typename, typename, typename, typename>
          class SurfaceSweepBaseCurve = Ss2::Default_subcurve,
          typename Subcurve_ = Default>
class Arr_construction_subcurve :
  public Arr_construction_subcurve_base<
    GeometryTraits_2, Event_, Allocator_,
    SurfaceSweepBaseCurve,
    typename Default::Get<Subcurve_,
                          Arr_construction_subcurve<GeometryTraits_2, Event_,
                                                    Allocator_,
                                                    SurfaceSweepBaseCurve,
                                                    Subcurve_> >::type>
{
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Event_                                        Event;
  typedef Allocator_                                    Allocator;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef Arr_construction_subcurve<Gt2, Event, Allocator,
                                    SurfaceSweepBaseCurve, Subcurve_>
                                                        Self;
  typedef typename Default::Get<Subcurve_, Self>::type  Subcurve;
  typedef Arr_construction_subcurve_base<Gt2, Event, Allocator,
                                         SurfaceSweepBaseCurve, Subcurve>
                                                        Base;

public:
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;

  typedef typename Base::Event_ptr                      Event_ptr;
  typedef typename Base::Halfedge_indices_list          Halfedge_indices_list;

  /*! Construct deafult. */
  Arr_construction_subcurve() {}

  /*! Constructor from an x-monotone curve. */
  Arr_construction_subcurve(X_monotone_curve_2& curve) :
    Base(curve)
  {}
};


} // namespace CGAL

#endif
