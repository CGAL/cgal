// Copyright (c) 2020  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_ORTHTREE_TRAITS_D_H
#define CGAL_ORTHTREE_TRAITS_D_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Dimension.h>

namespace CGAL
{

/*!
  \ingroup PkgOrthtreeTraits

  The class `Orthtree_traits_d` can be used as a template parameter of
  the `Orthtree` class.

  \tparam GeomTraits model of `Kernel`.
  \tparam DimensionTag specialization of `CGAL::Dimension_tag`.

  \cgalModels `OrthtreeTraits`
  \sa `CGAL::Orthtree`
  \sa `CGAL::Orthtree_traits_2`
  \sa `CGAL::Orthtree_traits_3`
*/

template <typename GeomTraits, typename DimensionTag>
struct Orthtree_traits_d
{
public:

  /// \name Types
  /// @{

  typedef DimensionTag Dimension; ///< Dimension type.
  typedef typename GeomTraits::FT FT; ///< Number type.
  typedef typename GeomTraits::Point_d Point_d; ///< Point type.
  typedef typename GeomTraits::Sphere_d Sphere_d; ///< Sphere type.
  typedef typename GeomTraits::Cartesian_const_iterator_d Cartesian_const_iterator_d; ///< An iterator over the %Cartesian coordinates.
  typedef std::array<FT, Dimension::value> Array; ///< Array type.

#ifdef DOXYGEN_RUNNING
  typedef unspecified_type Bbox_d; ///< Bounding box type.
#else
  class Bbox_d
  {
    Point_d m_min, m_max;
  public:

    Bbox_d (const Point_d& pmin, const Point_d& pmax)
      : m_min (pmin), m_max (pmax)
    { }

    const Point_d& min BOOST_PREVENT_MACRO_SUBSTITUTION () { return m_min; }
    const Point_d& max BOOST_PREVENT_MACRO_SUBSTITUTION () { return m_max; }
  };
#endif

  /*!
    Adjacency type.

    \note This type is used to identify adjacency directions with
    easily understandable keywords (left, right, up, etc.) and is thus
    mainly useful for `Orthtree_traits_2` and `Orthtree_traits_3`. In
    higher dimensions, such keywords do not exist and this type is
    simply an integer. Conversions from this integer to bitsets still
    works but do not provide any easier API for adjacency selection.
  */
  typedef int Adjacency;


#ifdef DOXYGEN_RUNNING
  /*!
    Functor with an operator to construct a `Point_d` from an `Array` object.
  */
  typedef unspecified_type Construct_point_d_from_array;
#else
  struct Construct_point_d_from_array
  {
    Point_d operator() (const Array& array) const
    {
      return Point_d (array.begin(), array.end());
    }
  };
#endif


#ifdef DOXYGEN_RUNNING
  /*!
    Functor with an operator to construct a `Bbox_d` from two `Array` objects (coordinates of minimum and maximum points).
  */
  typedef unspecified_type Construct_bbox_d;
#else
  struct Construct_bbox_d
  {
    Bbox_d operator() (const Array& min,
                       const Array& max) const
    {
      return Bbox_d (Point_d (min.begin(), min.end()),
                     Point_d (max.begin(), max.end()));
    }
  };
#endif

  /// @}

  /// \name Operations
  /// @{

  /*!
    Function used to construct an object of type `Construct_point_d_from_array`.
  */
  Construct_point_d_from_array construct_point_d_from_array_object() const
  { return Construct_point_d_from_array(); }

  /*!
    Function used to construct an object of type `Construct_bbox_d`.
  */
  Construct_bbox_d construct_bbox_d_object() const
  { return Construct_bbox_d(); }

  /// @}
};

}

#endif // CGAL_ORTHTREE_TRAITS_D_H
