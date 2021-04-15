// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_POLYLINE_2_H
#define CGAL_ARR_POLYLINE_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/Arr_geometry_traits/Polyline_2.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/Arr_geometry_traits/Polycurve_2.h>"
#include <CGAL/internal/deprecation_warning.h>

#include <CGAL/Arr_geometry_traits/Polycurve_2.h>

/*! \file
 * Header file for the polyline classes used by the
 * Arr_polycurve_basic_traits_2, Arr_polycurve_traits_2, and
 * Arr_polyline_traits_2 classes.
 */

#include <CGAL/Arr_geometry_traits/Polycurve_2.h>

namespace CGAL {

namespace polyline {

template <typename SubcurveType_2, typename PointType_2>
class Polyline_2 : public internal::Polycurve_2<SubcurveType_2, PointType_2> {
public:
  typedef SubcurveType_2                                Subcurve_type_2;
  typedef PointType_2                                   Point_type_2;

  typedef internal::Polycurve_2<Subcurve_type_2, Point_type_2> Base;

  typedef typename Base::Subcurve_type_2                Segment_type_2;
  typedef typename Base::size_type                      Segments_size_type;
  typedef typename Base::Subcurve_iterator              Segment_iterator;
  typedef typename Base::Subcurve_const_iterator        Segment_const_iterator;
  typedef typename Base::Subcurve_const_reverse_iterator
    Segment_const_reverse_iterator;

  /*! Construct default. */
  Polyline_2() : Base() {}

  /*! Construct from a subcurve. */
  Polyline_2(const Subcurve_type_2& subcurve) : Base(subcurve) {}

  /*! Construct from a range. */
  template <typename InputIterator>
  Polyline_2(InputIterator begin, InputIterator end) : Base(begin, end) {}

  /*! Obtain an iterator for the polycurve subcurves. */
  Segment_const_iterator begin_segments() const
  { return this->subcurves_begin(); }

  /*! Obtain a past-the-end iterator for the polycurve subcurves. */
  Segment_const_iterator end_segments() const
  { return this->subcurves_end(); }

  /*! Obtain a reverse iterator for the polycurve subcurves. */
  Segment_const_reverse_iterator rbegin_segments() const
  { return this->subcurves_rbegin(); }

  /*! Obtain a reverse past-the-end iterator for the polycurve points. */
  Segment_const_reverse_iterator rend_segments() const
  { return this->subcurves_rend() ; }

  /*! Obtain the number of subcurves that comprise the poyline.
   * \return The number of subcurves.
   */
  Segments_size_type number_of_segments() const
  { return this->number_of_subcurves(); }
};

template <typename SubcurveType_2, typename PointType_2>
class X_monotone_polyline_2 :
    public internal::X_monotone_polycurve_2<SubcurveType_2, PointType_2> {
public:
  typedef SubcurveType_2                                Subcurve_type_2;
  typedef PointType_2                                   Point_type_2;

  typedef internal::X_monotone_polycurve_2<Subcurve_type_2, Point_type_2> Base;

  typedef typename Base::Subcurve_type_2                Segment_type_2;
  typedef typename Base::size_type                      Segments_size_type;
  typedef typename Base::Subcurve_iterator              Segment_iterator;
  typedef typename Base::Subcurve_const_iterator        Segment_const_iterator;
  typedef typename Base::Subcurve_const_reverse_iterator
    Segment_const_reverse_iterator;

  /*! Construct default. */
  X_monotone_polyline_2() : Base() {}

  /*! Construct from a subcurve. */
  X_monotone_polyline_2(Subcurve_type_2 seg) : Base(seg) {}

  /*! Construct from a range.
   * Similar to the constructor of a general polycurve.
   * Like in the case of general polycurve, for the sake of backwards
   * compatibility we have to keep an implementation of construction
   * from a range of points. DO NOT USE THIS CONSTRUCTION.
   */
  template <typename InputIterator>
  X_monotone_polyline_2(InputIterator begin, InputIterator end) :
    Base(begin, end)
  {}

  /*! Obtain an iterator for the polycurve subcurves. */
  Segment_const_iterator begin_segments() const
  { return this->subcurves_begin(); }

  /*! Obtain a past-the-end iterator for the polycurve subcurves. */
  Segment_const_iterator end_segments() const
  { return this->subcurves_end(); }

  /*! Obtain a reverse iterator for the polycurve subcurves. */
  Segment_const_reverse_iterator rbegin_segments() const
  { return this->subcurves_rbegin(); }

  /*! Obtain a reverse past-the-end iterator for the polycurve points. */
  Segment_const_reverse_iterator rend_segments() const
  { return this->subcurves_rend() ; }

  /*! Obtain the number of subcurves that comprise the poyline.
   * \return The number of subcurves.
   */
  Segments_size_type number_of_segments() const
  { return this->number_of_subcurves(); }
};

} // namespace polyline
} //namespace CGAL

#endif
