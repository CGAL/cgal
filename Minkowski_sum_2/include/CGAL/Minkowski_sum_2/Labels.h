// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein   <wein_r@yahoo.com>

#ifndef CGAL_MINKOWSKI_SUM_LABELS_H
#define CGAL_MINKOWSKI_SUM_LABELS_H

#include <CGAL/license/Minkowski_sum_2.h>


namespace CGAL {

/*! \class
 * A label for endpoints of Minkowski-sum curves.
 */
class Point_label
{
private:

  unsigned int  _component;   // The index of the component this point
                              // belongs to.
  unsigned int  _index;       // Index of the point in the component.

public:

  /*! Default constructor. */
  Point_label () :
    _component (0),
    _index (0)
  {}

  /*!
   * Constructor.
   * \param component The index of the component.
   * \param index Index of the point within the component.
   */
  Point_label (unsigned int component,
               unsigned int index) :
    _component (component),
    _index (index)
  {
    CGAL_precondition (component != 0);
  }

  /*! Check whether the label is valid. */
  bool is_valid () const
  {
    return (_component != 0);
  }

  /*! Check two labels for equality. */
  bool operator== (const Point_label& label) const
  {
    if (_component == 0 || label._component == 0)
      return (false);

    return (_component == label._component &&
	    _index == label._index);
  }

  /*! Get the component. */
  unsigned int component () const
  {
    return (_component);
  }

  /*! Get the index within the component. */
  unsigned int index () const
  {
    return (_index);
  }

  /*! Set the component. */
  void set_component (unsigned int component)
  {
    CGAL_precondition (component != 0);

    _component = component;
    return;
  }

  /*! Set the index. */
  void set_index (unsigned int index)
  {
    _index = index;
    return;
  }
};

/*! \class
 * A label for Minkowski-sum x-monotone curves.
 */
class X_curve_label
{
private:

  unsigned int  _component;   // The index of the component the x-monotone
                              // curve belongs to.
  unsigned int  _index;       // Index of the curve in the component.
  bool          _is_last;     // Is the curve the last one in its component.
  int           _right;       // The number of curves directed from left to
                              // right that are associated with the curve.
  int           _left;        // The number of curves directed from right to
                              // left that are associated with the curve.

public:

  /*! Default constructor. */
  X_curve_label () :
    _component (0),
    _index (0),
    _is_last (false),
    _right (0),
    _left (0)
  {}

  /*!
   * Constructor.
   * \param is_directed_right Is the curve directed from left to right.
   * \param component The index of the component.
   * \param index Index of the curve in the component.
   * \param is_last Is this the last curve of the component.
   */
  X_curve_label (bool is_directed_right,
                 unsigned int component,
                 unsigned int index,
                 bool is_last = false) :
    _component (component),
    _index (index),
    _is_last (is_last),
    _right (is_directed_right ? 1 : 0),
    _left (is_directed_right ? 0 : 1)
  {
    CGAL_precondition (component != 0);
  }

  /*!
   * Construct a merged curve label.
   * \param label1 The first label.
   * \param label2 The second label.
   */
  X_curve_label (const X_curve_label& label1,
                 const X_curve_label& label2) :
    _component (0),
    _index (0),
    _is_last (false),
    _right (label1._right + label2._right),
    _left (label1._left + label2._left)
  {}

  /*! Check whether the label is valid. */
  bool is_valid () const
  {
    return (_component != 0);
  }

  /*! Check two labels for equality. */
  bool operator== (const X_curve_label& label) const
  {
    if (_component == 0)
      return (false);

    return (_component == label._component &&
	    _index == label._index);
  }

  /*! Check whether the given label is the predecessor of this label. */
  bool is_prev (const X_curve_label& label) const
  {
    if (_component == 0)
      return (false);

    return (_component == label._component &&
	    (label._index + 1 == _index ||
	     (label._is_last && _index == 0)));
  }

  /*! Check whether the given label is the succcessor of this label. */
  bool is_next (const X_curve_label& label) const
  {
    if (_component == 0)
      return (false);

    return (_component == label._component &&
	    (_index + 1 == label._index ||
	     (_is_last && label._index == 0)));
  }

  /*!
   * Check if the given label is the successor of this label
   * (or vice-versa).
   */
  bool is_adjacent (const X_curve_label& label) const
  {
    if (_component == 0)
      return (false);

    // Two curves can be adjacent only if they belong to the same component.
    if (_component != label._component)
      return (false);

    // Check if one label is adjacent to the next.
    return ((_index + 1 == label._index) ||
            (label._index + 1 == _index) ||
            (_is_last && label._index == 0) ||
            (label._is_last && _index == 0));
  }

  /*! Get the number of curves directed from left to right. */
  int right_count () const
  {
    return (_right);
  }

  /*! Get the number of curves directed from right to left. */
  int left_count () const
  {
    return (_left);
  }

  /*! Get the component. */
  unsigned int component () const
  {
    return (_component);
  }

  /*! Get the index within the component. */
  unsigned int index () const
  {
    return (_index);
  }

  /*! Check whether this is the last curve in the component. */
  bool is_last () const
  {
    return (_is_last);
  }

  /*! Set the component. */
  void set_component (unsigned int component)
  {
    CGAL_precondition (component != 0);

    _component = component;
    return;
  }

  /*!
   * Set the curve index within the component.
   * \param index The index in the component.
   * \param is_last Is this the last curve of the component.
   */
  void set_index (unsigned int index, bool is_last = false)
  {
    _index = index;
    _is_last = is_last;
    return;
  }

  /*! Get the flag value. */
  bool get_flag () const
  {
    return (_is_last);
  }

  /*! Set the flag value. */
  void set_flag (bool flag)
  {
    _is_last = flag;
    return;
  }
};

} //namespace CGAL

#endif
