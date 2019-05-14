// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Philipp Möller

#ifndef CGAL_INTERSECTION_TRAITS_2_H
#define CGAL_INTERSECTION_TRAITS_2_H

#include <CGAL/Intersection_traits.h>

#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <vector>

#if !(CGAL_INTERSECTION_VERSION < 2)

namespace CGAL {

CGAL_INTERSECTION_TRAITS_2(Line_2, Line_2, Point_2, Line_2)

CGAL_INTERSECTION_TRAITS_2(Segment_2, Line_2, Point_2, Segment_2)
CGAL_INTERSECTION_TRAITS_2(Line_2, Segment_2, Point_2, Segment_2)

CGAL_INTERSECTION_TRAITS_2(Segment_2, Segment_2, Point_2, Segment_2)

CGAL_INTERSECTION_TRAITS_2(Ray_2, Line_2, Point_2, Ray_2)
CGAL_INTERSECTION_TRAITS_2(Line_2, Ray_2, Point_2, Ray_2)

CGAL_INTERSECTION_TRAITS_2(Ray_2, Segment_2, Point_2, Segment_2)
CGAL_INTERSECTION_TRAITS_2(Segment_2, Ray_2, Point_2, Segment_2)

CGAL_INTERSECTION_TRAITS_3(Ray_2, Ray_2, Point_2, Segment_2, Ray_2)

CGAL_INTERSECTION_TRAITS_2(Triangle_2, Line_2, Point_2, Segment_2)
CGAL_INTERSECTION_TRAITS_2(Line_2, Triangle_2, Point_2, Segment_2)

CGAL_INTERSECTION_TRAITS_2(Triangle_2, Segment_2, Point_2, Segment_2)
CGAL_INTERSECTION_TRAITS_2(Segment_2, Triangle_2, Point_2, Segment_2)

CGAL_INTERSECTION_TRAITS_2(Triangle_2, Ray_2, Point_2, Segment_2)
CGAL_INTERSECTION_TRAITS_2(Ray_2, Triangle_2, Point_2, Segment_2)

template<typename K>
struct Intersection_traits<K, typename K::Triangle_2, typename K::Triangle_2>  {
  typedef typename 
  boost::variant< typename K::Point_2, typename K::Segment_2,
                  typename K::Triangle_2, typename std::vector< typename K::Point_2 > > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};


CGAL_INTERSECTION_TRAITS_2(Iso_rectangle_2, Line_2, Point_2, Segment_2)
CGAL_INTERSECTION_TRAITS_2(Line_2, Iso_rectangle_2, Point_2, Segment_2)

CGAL_INTERSECTION_TRAITS_2(Iso_rectangle_2, Segment_2, Point_2, Segment_2)
CGAL_INTERSECTION_TRAITS_2(Segment_2, Iso_rectangle_2, Point_2, Segment_2)

CGAL_INTERSECTION_TRAITS_2(Iso_rectangle_2, Ray_2, Point_2, Segment_2)
CGAL_INTERSECTION_TRAITS_2(Ray_2, Iso_rectangle_2, Point_2, Segment_2)

// undocumented

// Variants of one for backwards compatibility
template<typename K>
struct Intersection_traits<K, typename K::Iso_rectangle_2, typename K::Iso_rectangle_2>  {
  typedef typename boost::variant<typename K::Iso_rectangle_2> variant_type;
  typedef boost::optional<variant_type> result_type;
};


// Point_2 is special
template<typename K, typename B>
struct Intersection_traits<K, typename K::Point_2, B> {
  typedef typename boost::variant<typename K::Point_2> variant_type;
  typedef boost::optional<variant_type> result_type;
};

template<typename K, typename A>
struct Intersection_traits<K, A, typename K::Point_2> {
  typedef typename boost::variant<typename K::Point_2> variant_type;
  typedef boost::optional<variant_type> result_type;
};


template<typename K>
struct Intersection_traits<K, typename K::Iso_rectangle_2, typename K::Triangle_2>
{
  typedef typename boost::variant<typename K::Segment_2, typename K::Triangle_2, 
                                  typename K::Point_2, 
                                  typename std::vector< typename K::Point_2 > > variant_type;
  typedef typename boost::optional < variant_type > result_type;
};  

template<typename K>
struct Intersection_traits<K, typename K::Triangle_2, typename K::Iso_rectangle_2>
  : public Intersection_traits<K, typename K::Iso_rectangle_2, typename K::Triangle_2> {};

} // namespace CGAL

#endif

#endif /* CGAL_INTERSECTION_TRAITS_2_H */

