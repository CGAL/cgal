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
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARR_CONSTRUCTION_SUBCURVE_H
#define CGAL_ARR_CONSTRUCTION_SUBCURVE_H

/*! \file
 * Definition of the Arr_construction_subcurve class-template.
 */

#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>

namespace CGAL {

/*! \class Arr_construction_subcurve 
 *  
 * a class that holds information about a curve that is added to 
 * the arrangement.
 * In addition to the information that is contained in Sweep_line_subcurve,
 * when an arrangement is constructed, a pointer to the last handled event  
 * on the curve is stored. This information is used to retrieve
 * hints when a subcurve of this curve is inserted into the planar map.
 *
 * Inherits from `Sweep_line_subcurve`
 * \sa `Sweep_line_subcurve`
 */

template <typename Traits_>
class Arr_construction_subcurve : public Sweep_line_subcurve<Traits_>
{
public:
  typedef Traits_                                    Traits_2;
  typedef typename Traits_2::Point_2                 Point_2;
  typedef typename Traits_2::X_monotone_curve_2      X_monotone_curve_2;

  typedef Sweep_line_subcurve<Traits_2>              Base;
  typedef Arr_construction_subcurve<Traits_2>        Self;

  typedef typename Base::Status_line_iterator        Status_line_iterator;
  typedef void*                                      Event_ptr;
  typedef std::list<unsigned int>                    Halfedge_indices_list;

protected:

  // Data members:
  Event_ptr m_lastEvent;     // The last event that was handled on the curve.

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
  /*! Deafult constructor. */
  Arr_construction_subcurve() :
    Base(),
    m_lastEvent(0),
    m_index(0)
  {}

  /*! Constructor from an x-monotone curve. */
  Arr_construction_subcurve(X_monotone_curve_2& curve) :
    Base( curve),
    m_lastEvent(0),
    m_index(0)
  {}

  /*! Initialize the curve. */
  void init(const X_monotone_curve_2& curve) { Base::init(curve); }

  /*! Set the event associated with the left end of the subcurve. */
  template<class SweepEvent>
  void set_left_event(SweepEvent* left)
  {
    Base::set_left_event(left);
    m_lastEvent = left;
  }

  /*! Set the last event on the subcurve. */
  void set_last_event(Event_ptr e) { m_lastEvent = e; }

  /*! Get the last event. */
  Event_ptr last_event() const { return m_lastEvent; }

  /*! Get the subcurve index. */
  unsigned int index() const { return m_index; }

  /*! Set the subcurve index. */
  void set_index(unsigned int i) { m_index = i; }

  /*! Check if the index is valid. */
  bool has_valid_index() const { return (m_index != 0); }

  /*! Add an index of a halfedge below the subcurve. */
  void add_halfedge_index(unsigned int i) { m_halfedge_indices.push_back(i); }

  /*! Clear the indices of the halfedges below the subcurve. */
  void clear_halfedge_indices() { m_halfedge_indices.clear(); }

  /*! Check if there are any halfedges below the subcurve. */
  bool has_halfedge_indices() const { return (!m_halfedge_indices.empty()); }

  /*! Get the indices of the halfedges below the subcurve. */
  Halfedge_indices_list& halfedge_indices_list() { return m_halfedge_indices; }
};


} //namespace CGAL

#endif 
