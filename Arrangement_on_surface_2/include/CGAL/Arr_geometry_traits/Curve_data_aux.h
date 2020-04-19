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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>

#ifndef CGAL_CURVE_DATA_AUX_H
#define CGAL_CURVE_DATA_AUX_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Definition of auxiliary classes for the Arr_curve_data_traits_2<> template.
 */

namespace CGAL {

/*!
 * \struct A simple merge functor.
 */
template <class TYPE>
struct _Default_merge_func
{
  const TYPE& operator() (const TYPE& obj1, const TYPE& /* obj2 */)
  {
    return (obj1);
  }
};

/*!
 * \struct A simple convertor from one type to another.
 */
template <class TYPE_FROM, class TYPE_TO>
struct _Default_convert_func
{
  TYPE_TO operator() (const TYPE_FROM& obj)
  {
    return (obj);
  }
};

/*! \class
 * A template for extending the base curve type with a data field.
 */
template <class BaseCurveType, class Data>
class _Curve_data_ex : public BaseCurveType
{
private:
  Data    m_data;

public:

  /*! Default constructor. */
  _Curve_data_ex ()
  {}

  /*!
   * Construct an extended curve from a base curve.
   * \param cv The base curve.
   */
  _Curve_data_ex (const BaseCurveType& cv) :
    BaseCurveType (cv)
  {}

  /*!
   * Construct an extended curve from a base curve and a data object.
   * \param cv The base curve.
   * \param data The data object.
   */
  _Curve_data_ex (const BaseCurveType& cv, const Data& data) :
    BaseCurveType (cv),
    m_data (data)
  {}

  /*!
   * Get the data object (const version).
   * \return The data object associated with the curve.
   */
  const Data& data () const
  {
    return (m_data);
  }

  /*!
   * Get the data object (non-const version).
   * \return The data object associated with the curve.
   */
  Data& data ()
  {
    return (m_data);
  }

  /*!
   * Set the curve data.
   * \param data The data object to be associated with the curve.
   */
  void set_data (const Data& data)
  {
    m_data = data;
    return;
  }
};

} //namespace CGAL

#endif
