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
  Curve_handle insert (const Curve_2& cv,
                       const PointLocation& pl)
  {
    return (p_arr->_insert (cv, pl));
  }

  /*!
   * Insert a curve into the arrangement, using the default "walk"
   * point-location strategy.
   * \param cv The curve to be inserted.
   * \return A handle to the inserted curve.
   */
  Curve_handle insert (const Curve_2& cv)
  {
    return (p_arr->_insert (cv));
  }

  /*!
   * Insert a range of curves into the arrangement.
   * \param begin An iterator pointing to the first curve in the range.
   * \param end A past-the-end iterator for the last curve in the range.
   */
  template <class InputIterator>
  void insert (InputIterator begin, InputIterator end)
  {
    p_arr->_insert (begin, end);
    return;
  }

  /*!
   * Remove a curve from the arrangement (remove all the edges it induces).
   * \param ch A handle to the curve to be removed.
   * \return The number of removed edges.
   */
  Size remove (Curve_handle ch)
  {
    return (p_arr->_remove (ch));
  }

  /*!
   * Overlay the two given arrangements.
   * \param arr1 The first arrangement.
   * \param arr2 The second arrangement.
   * \param overlay_tr An overlay-traits class.
   */
  /*
  template <class Dcel1, class Dcel2, class OverlayTraits>
  void overlay (const Arrangement_with_history_2<Traits_2, Dcel1>& arr1,
                const Arrangement_with_history_2<Traits_2, Dcel2>& arr2,
                OverlayTraits& overlay_tr)
  {
    p_arr->_overlay (arr1, arr2, overlay_tr);
    return;
  }
  */
  //@}
};

CGAL_END_NAMESPACE

#endif
