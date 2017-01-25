// Copyright (c) 2005,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>

#ifndef CGAL_ARR_WITH_HISTORY_ACCESSOR_H
#define CGAL_ARR_WITH_HISTORY_ACCESSOR_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Definition of the Arr_with_history_accessor<Arrangement> class.
 */

namespace CGAL {

/*! \class
 * A class that provides access to some of the internal methods of the
 * Arrangement_on_surface_with_history_2 class.
 * Used mostly by the global functions that operate on arrangments with
 * history objects.
 */
template <class ArrWithHistory_>
class Arr_with_history_accessor
{
public:

  typedef ArrWithHistory_                        Arrangement_with_history_2;
  typedef Arr_with_history_accessor<Arrangement_with_history_2> Self;

  typedef typename Arrangement_with_history_2::Geometry_traits_2
                                                             Geometry_traits_2;
  typedef typename Arrangement_with_history_2::Topology_traits
                                                             Topology_traits;
  typedef typename Arrangement_with_history_2::Size          Size;
  typedef typename Arrangement_with_history_2::Point_2       Point_2;
  typedef typename Arrangement_with_history_2::Curve_2       Curve_2;
  typedef typename Arrangement_with_history_2::Curve_handle  Curve_handle;
  typedef typename Arrangement_with_history_2::Halfedge_handle
                                                             Halfedge_handle;

private:

  Arrangement_with_history_2  *p_arr;           // The associated arrangement.

public:

  /*! Constructor with an associated arrangement. */
  Arr_with_history_accessor (Arrangement_with_history_2& arr) :
    p_arr (&arr)
  {}

  /// \name Accessing the private insertion and removal functions.
  //@{

  /*!
   * Insert a curve into the arrangement.
   * \param cv The curve to be inserted.
   * \param pl a point-location object.
   * \return A handle to the inserted curve.
   */
  template <class PointLocation>
  Curve_handle insert_curve (const Curve_2& cv,
                             const PointLocation& pl)
  {
    return (p_arr->_insert_curve (cv, pl));
  }

  /*!
   * Insert a curve into the arrangement, using the default "walk"
   * point-location strategy.
   * \param cv The curve to be inserted.
   * \return A handle to the inserted curve.
   */
  Curve_handle insert_curve (const Curve_2& cv)
  {
    return (p_arr->_insert_curve (cv));
  }

  /*!
   * Insert a range of curves into the arrangement.
   * \param begin An iterator pointing to the first curve in the range.
   * \param end A past-the-end iterator for the last curve in the range.
   */
  template <class InputIterator>
  void insert_curves (InputIterator begin, InputIterator end)
  {
    p_arr->_insert_curves (begin, end);
    return;
  }

  /*!
   * Remove a curve from the arrangement (remove all the edges it induces).
   * \param ch A handle to the curve to be removed.
   * \return The number of removed edges.
   */
  Size remove_curve (Curve_handle ch)
  {
    return (p_arr->_remove_curve (ch));
  }
  //@}

  /// \name Assigning input curves.
  //@{

  /*!
   * Create a new curve in the curves' container.
   * \param cv The source curve.
   * \return A handle for the curve.
   */
  Curve_handle new_curve (const Curve_2& cv)
  {
    // Allocate an extended curve (with an initially empty set of edges)
    // and store it in the curves' list.
    typename Arrangement_with_history_2::Curve_halfedges  *p_cv = 
                                           p_arr->m_curves_alloc.allocate (1);
    
    p_arr->m_curves_alloc.construct (p_cv, cv);
    p_arr->m_curves.push_back (*p_cv);
    
    // Return a handle to the inserted curve (the last in the list).
    Curve_handle       ch = p_arr->m_curves.end();
    return (--ch);
  }

  /*!
   * Set the cross-pointers between a curve and an edge.
   * \param ch A handle for the curve.
   * \param he A halfedge handle representing the edge.
   */
  void connect_curve_edge (Curve_handle ch, Halfedge_handle he)
  {
    // Add the edge to the list of cv's induced edges.
    typename Arrangement_with_history_2::Curve_halfedges&  cv = *ch;
    cv._insert (he);

    // Add the curve to the set of he's inducing curves. 
    he->curve().data().insert (&cv);
  
    return;
  }  
};

} //namespace CGAL

#endif
