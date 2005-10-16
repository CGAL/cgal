// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>

#ifndef CGAL_ARR_WITH_HISTORY_ACCESSOR_H
#define CGAL_ARR_WITH_HISTORY_ACCESSOR_H

/*! \file
 * Definition of the Arr_with_history_accessor<Arrangement> class.
 */

CGAL_BEGIN_NAMESPACE

/*! \class
 * A class that provides access to some of the internal methods of the
 * Arrangement_with_history_2 class.
 * Used mostly by the global functions that operate on arrangments with
 * history objects.
 */
template <class Arr_with_hist_>
class Arr_with_history_accessor
{
public:

  typedef Arr_with_hist_                         Arrangement_with_history_2;
  typedef Arr_with_history_accessor<Arrangement_with_history_2> Self;

  typedef typename Arrangement_with_history_2::Traits_2      Traits_2;
  typedef typename Arrangement_with_history_2::Size          Size;
  typedef typename Arrangement_with_history_2::Point_2       Point_2;
  typedef typename Arrangement_with_history_2::Curve_2       Curve_2;
  typedef typename Arrangement_with_history_2::Curve_handle  Curve_handle;

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
};

CGAL_END_NAMESPACE

#endif
