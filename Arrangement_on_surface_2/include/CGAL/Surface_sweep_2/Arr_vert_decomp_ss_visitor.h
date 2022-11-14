// Copyright (c) 2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ron Wein <wein@post.tau.ac.il>
//            Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_ARR_VERT_DECOMP_SS_VISITOR_H
#define CGAL_ARR_VERT_DECOMP_SS_VISITOR_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_vert_decomp_ss_visitor class-template.
 */

#include <boost/variant.hpp>

namespace CGAL {

#include <CGAL/Surface_sweep_2/Default_visitor_base.h>
#include <CGAL/Default.h>

/*! \class Arr_vert_decomp_ss_visitor
 *
 * A sweep-line visitor for performing vertical decomposition on an
 * arrangement embedded on a surface.
 */
template <typename Helper_, typename OutputIterator, typename Visitor_ = Default>
class Arr_vert_decomp_ss_visitor :
  public Ss2::Default_visitor_base<typename Helper_::Geometry_traits_2,
                                   typename Helper_::Event,
                                   typename Helper_::Subcurve,
                                   typename Helper_::Allocator,
                                   typename Default::Get<
                                     Visitor_,
                                     Arr_vert_decomp_ss_visitor<
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
  typedef Arr_vert_decomp_ss_visitor<Helper, Output_iterator, Visitor_>
                                                        Self;
  typedef typename Default::Get<Visitor_, Self>::type   Visitor;
  typedef typename Ss2::Default_visitor_base<Gt2, Event, Subcurve, Allocator,
                                             Visitor>   Base;
public:
  typedef typename Helper::Arrangement_2                Arrangement_2;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;

  typedef boost::variant<Vertex_const_handle, Halfedge_const_handle,
                         Face_const_handle>
    Cell_type;
  typedef boost::optional<Cell_type>                    Vert_type;
  typedef std::pair<Vert_type, Vert_type>               Vert_pair;
  typedef std::pair<Vertex_const_handle, Vert_pair>     Vert_entry;

protected:
  typedef typename Base::Status_line_iterator           Status_line_iterator;
  //typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;

  // Data members:
  Helper m_helper;                      // The helper class.

  const typename Arrangement_2::Geometry_traits_2* m_traits;
                                        // The traits class.

  const Vertex_const_handle invalid_vh;
                                        // An invalid vertex handle.

  Vertex_const_handle m_prev_vh;        // The previous vertex.
  Vert_type m_prev_obj_below;           // The object this vertex sees below it.
  Vert_type m_prev_obj_above;           // The object this vertex sees above it.

  Output_iterator* m_out;               // An output iterator for the result.

public:
  /*! Constructor.
   * \param arr The arrangement.
   * \param oi A pointer to the output iterator that will store the result.
   */
  Arr_vert_decomp_ss_visitor(const Arrangement_2* arr, OutputIterator* oi) :
    m_helper(arr),
    m_traits(arr->geometry_traits()),
    invalid_vh(),
    m_out(oi)
  {}

  /* A notification issued before the sweep process starts. */
  void before_sweep();

  /*! A notification invoked after the sweep-line finishes handling the given
   * event.
   * \param event The event.
   * \param above An iterator to the sweep-line subcurves lying right above
   *              (or on) the event point.
   * \param on_above Whether the event is locates on the subcurve above it.
   */
  bool after_handle_event(Event* event,
                          Status_line_iterator above,
                          bool on_above);

  /*! A notification issued when the sweep process is over.
   */
  void after_sweep();
};

//-----------------------------------------------------------------------------
// Member-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// A notification issued before the sweep process starts.
//
template <typename Hlpr, typename OutIt, typename Vis>
void Arr_vert_decomp_ss_visitor<Hlpr, OutIt, Vis>::before_sweep()
{
  // Notify the helper that the sweep process now starts.
  m_helper.before_sweep();

  // Set an invalid previous vertex.
  m_prev_vh = Vertex_const_handle();
}

//-----------------------------------------------------------------------------
// A notification invoked after the sweep-line finishes handling the given
// event.
//
template <class Hlpr, class OutIt, typename Vis>
bool Arr_vert_decomp_ss_visitor<Hlpr, OutIt, Vis>::
after_handle_event(Event* event,
                   Status_line_iterator above, bool /* on_above */)
{
  // Notify the helper on the event.
  m_helper.after_handle_event(event);

  // We are only interested in events associated with valid points:
  if (! event->is_closed()) return true;

  // Get the vertex handle associated with the current event (stored with
  // the point).
  Vertex_const_handle vh = event->point().vertex_handle();
  Vert_type obj_above, obj_below;

  // Check the feature from above.
  if (above == this->status_line_end()) {
    // There is no concrete subcurve above the current event point, so we use
    // the helper class to obtain the object above.
    obj_above = m_helper.top_object();
  }
  else {
    // We have a valid subcurve above the event: get its halfedge handle
    // and associate it with the vertex.
    obj_above = Vert_type((*above)->last_curve().halfedge_handle());
  }

  // Check if the previous vertex we handled has the same x-coordinate
  // as the current one (it lies vertically below the current vertex).
  const bool prev_same_x =
    (m_prev_vh != invalid_vh &&
     m_traits->compare_x_2_object()(vh->point(),
                                    m_prev_vh->point()) == EQUAL);

  // Decrement the status-line iterator to reach the subcurve below the
  // event point. If the number of right subcurves associated with the
  // event is n_right, we decrement the iterator n_right times, then
  // check if it is possible to further decrement it.
  Status_line_iterator below = above;
  const size_t n_right = event->number_of_right_curves();
  size_t k;

  for (k = 0; k < n_right; k++) --below;

  if (below == this->status_line_begin()) {
    if (prev_same_x) {
      // The previous vertex is vertically below the current vertex,
      // so we update their respective entries in the output map.
      // We first check if the two vertices are connected (by a vertical
      // segment) - if so, they "see" empty features.
      bool  vert_connected = false;

      if (! vh->is_isolated()) {
        Halfedge_around_vertex_const_circulator circ, first;

        first = circ = vh->incident_halfedges();
        do {
          if (circ->source() == m_prev_vh) vert_connected = true;
          ++circ;
        } while (! vert_connected && circ != first);
      }

      if (! vert_connected) {
        obj_below = Vert_type(m_prev_vh);
        m_prev_obj_above = Vert_type(vh);
      }
      else {
        obj_below = Vert_type();
        m_prev_obj_above = Vert_type();
      }
    }
    else {
      // There is no concrete subcurve below the current event point, so we
      // use the helper class to obtain the object below.
      obj_below = m_helper.bottom_object();
    }
  }
  else {
    // Decrement the iterator once more in order to reach the subcurve below
    // the current event.
    --below;

    // If the previous event has the same x-coordinate as the current one,
    // we check whether it lies below or above the subcurve below. We do this
    // to distinguish between the following scenarios:
    //
    //
    //           *                      *
    //
    //         --------                 x
    //
    //           x                   ---------
    //
    if (prev_same_x &&
        (m_traits->compare_y_at_x_2_object()
            (m_prev_vh->point(),
             (*below)->last_curve().base()) != SMALLER))
    {
        // The previous vertex is vertically below the current vertex,
        // and above the subcurve that lies below vh. We therefore update
        // the respective entries in the output map accordingly.
        // We first check if the two vertices are connected (by a vertical
        // segment) - if so, they "see" empty features.
        bool vert_connected = false;

        if (! vh->is_isolated()) {
          Halfedge_around_vertex_const_circulator circ, first;

          first = circ = vh->incident_halfedges();
          do {
            if (circ->source() == m_prev_vh)
              vert_connected = true;
            ++circ;
          } while (! vert_connected && circ != first);
        }

        if (! vert_connected) {
          obj_below = Vert_type(m_prev_vh);
          m_prev_obj_above = Vert_type(vh);
        }
        else {
          obj_below = Vert_type();
          m_prev_obj_above = Vert_type();
        }
      }
      else {
        // Get the halfedge handle of the subcurve below the current event and
        // associate it with its vertex.
        obj_below = Vert_type((*below)->last_curve().halfedge_handle());
      }
    }

    // We can now create the entry for the previous vertex, as we are not
    // going to change the identity of the features below or above it.
    if (m_prev_vh != Vertex_const_handle()) {
      *(*m_out) = Vert_entry(m_prev_vh, Vert_pair(m_prev_obj_below,
                                                  m_prev_obj_above));
      ++(*m_out);
    }

    // We are done with the current vertex, but we cannot create the entry
    // yet - so we store it for the next event.
    m_prev_vh = vh;
    m_prev_obj_below = obj_below;
    m_prev_obj_above = obj_above;

    // It is safe to deallocate the event.
    return true;
}

//-----------------------------------------------------------------------------
// A notification issued when the sweep process is over.
//
template <typename Hlpr, typename OutIt, typename Vis>
void Arr_vert_decomp_ss_visitor<Hlpr, OutIt, Vis>::after_sweep()
{
  // Create an entry for the last vertex (the xy-largest one).
  if (m_prev_vh != invalid_vh) {
    *(*m_out) = Vert_entry(m_prev_vh, Vert_pair(m_prev_obj_below,
                                                m_prev_obj_above));
    ++(*m_out);
  }
}

} // namespace CGAL

#endif
