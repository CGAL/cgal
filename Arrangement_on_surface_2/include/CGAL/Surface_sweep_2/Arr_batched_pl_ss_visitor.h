// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_BATCHED_PL_SS_VISITOR_H
#define CGAL_ARR_BATCHED_PL_SS_VISITOR_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_batched_pl_ss_visitor class-template.
 */

#include <CGAL/Arr_point_location_result.h>
#include <CGAL/Surface_sweep_2/Default_visitor_base.h>
#include <CGAL/Default.h>

namespace CGAL {

namespace Ss2 = Surface_sweep_2;

/*! \class Arr_batched_pl_ss_visitor
 *
 * A surface-sweep visitor for performing batched point-location queries on an
 * arrangement embedded on a surface.
 */
template <typename Helper_, typename OutputIterator, typename Visitor_ = Default>
class Arr_batched_pl_ss_visitor :
  public Ss2::Default_visitor_base<typename Helper_::Geometry_traits_2,
                                   typename Helper_::Event,
                                   typename Helper_::Subcurve,
                                   typename Helper_::Allocator,
                                   typename Default::Get<
                                     Visitor_,
                                     Arr_batched_pl_ss_visitor<
                                       Helper_, OutputIterator,
                                       Visitor_> >::type>
{
public:
  typedef Helper_                                       Helper;
  typedef OutputIterator                                Output_iterator;

  typedef typename Helper::Geometry_traits_2            Geometry_traits_2;
  typedef typename Helper::Event                        Event;
  typedef typename Helper::Subcurve                     Subcurve;
  typedef typename Helper::Allocator                    Allocator;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef Arr_batched_pl_ss_visitor<Helper, Output_iterator, Visitor_>
                                                        Self;
  typedef typename Default::Get<Visitor_, Self>::type   Visitor;
  typedef typename Ss2::Default_visitor_base<Gt2, Event, Subcurve, Allocator,
                                             Visitor>   Base;

public:
  typedef typename Helper::Arrangement_2                Arrangement_2;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;

protected:
  typedef typename Base::Status_line_iterator           Status_line_iterator;

  typedef Arr_point_location_result<Arrangement_2>      Pl_result;
  typedef typename Pl_result::Type                      Pl_result_type;

  // Data members:
  Helper m_helper;              // The helper class.
  Output_iterator& m_out;       // An output iterator for the result.

  template<typename T>
  Pl_result_type pl_make_result(T t) { return Pl_result::make_result(t); }

  inline Pl_result_type pl_default_result()
  { return Pl_result::default_result(); }

public:
  /*! Constructor.
   * \param arr The arrangement.
   * \param oi A pointer to the output iterator that will store the result.
   */
  Arr_batched_pl_ss_visitor(const Arrangement_2* arr, Output_iterator& oi) :
    m_helper(arr),
    m_out(oi)
  {}

  /* A notification issued before the sweep process starts. */
  void before_sweep();

  /*! A notification invoked after the surface-sweep finishes handling the given
   * event.
   * \param event The event.
   * \param above An iterator to the surface-sweep subcurves lying right above
   *              (or on) the event point.
   * \param on_above Whether the event is locates on the subcurve above it.
   */
  bool after_handle_event(Event* event,
                          Status_line_iterator above,
                          bool on_above);
};

//-----------------------------------------------------------------------------
// Member-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// A notification issued before the sweep process starts.
//
  template <typename Hlpr, typename OutIt, typename Vis>
  void Arr_batched_pl_ss_visitor<Hlpr, OutIt, Vis>::before_sweep()
{
  // We just have to notify the helper that the sweep process now starts.
  m_helper.before_sweep();
}

//-----------------------------------------------------------------------------
// A notification invoked after the surface-sweep finishes handling the given
// event.
//
template <typename Hlpr, typename OutIt, typename Vis>
bool Arr_batched_pl_ss_visitor<Hlpr, OutIt, Vis>::
after_handle_event(Event* event, Status_line_iterator above, bool on_above)
{
  // Notify the helper on the event.
  m_helper.after_handle_event(event);

  // We are only interested in events associated with query points:
  if (! event->is_query()) return true;

  // Check on what kind of feature does the current query point lie.
  if (event->is_action()) {
    // The query point coincides with an isolated arrangement vertex:
    Vertex_const_handle  vh = event->point().vertex_handle();
    *m_out++ = std::make_pair(event->point().base(), pl_make_result(vh));
    return true;
  }

  if (event->has_right_curves() || event->has_left_curves()) {
    // The point is located on an arrangement vertex:
    Vertex_const_handle  vh;

    if (event->has_right_curves()) {
      // Take the target vertex of one of the halfedges to the right.
      Subcurve* sc = *(event->right_curves_begin());
      Halfedge_const_handle he = sc->last_curve().halfedge_handle();

      vh = he->target();
    }
    else {
      // Take the source vertex of one of the halfedges to the right.
      Subcurve* sc = *(event->left_curves_begin());
      Halfedge_const_handle he = sc->last_curve().halfedge_handle();

      vh = he->source();
    }

    *m_out++ = std::make_pair(event->point().base(), pl_make_result(vh));
    return true;
  }

  if (above == this->status_line_end()) {
    // There are no valid edges above the query point, so we use the helper
    // class to obtain the current top face.
    *m_out++ = std::make_pair(event->point().base(),
                              pl_make_result(m_helper.top_face()));
    return true;
  }

  if (on_above) {
    // The query point lies on the halfedge associated with the subcurve
    // that the status-line iterator refers to.
    Halfedge_const_handle  he = (*above)->last_curve().halfedge_handle();
    *m_out++ = std::make_pair(event->point().base(), pl_make_result(he));
    return true;
  }

  // If we reached here, the status-line iterator refers to a halfedge above
  // the query point, such that the query point is located in the incident
  // face of this halfedge.
  Halfedge_const_handle he = (*above)->last_curve().halfedge_handle();
  *m_out++ = std::make_pair(event->point().base(), pl_make_result(he->face()));
  return true;
}

} // namespace CGAL

#endif
