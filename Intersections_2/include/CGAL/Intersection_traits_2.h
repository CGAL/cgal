// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Philipp Möller

#ifndef CGAL_INTERSECTION_TRAITS_H
#define CGAL_INTERSECTION_TRAITS_H

#include <boost/variant.hpp>
#include <boost/optional.hpp>



#define CGAL_INTERSECTION_TRAITS_2(A, B, R1, R2)                          \
  template<typename K>                                                  \
  struct Intersection_traits_2<K, typename K::A, typename K::B>  {        \
    typedef typename boost::variant<typename K::R1, typename K::R2 >    \
                     variant_type;                                      \
    typedef typename boost::optional< variant_type > result_type;       \
  };  

#define CGAL_INTERSECTION_TRAITS_3(A, B, R1, R2, R3)                      \
  template<typename K>                                                  \
  struct Intersection_traits_2<K, typename K::A, typename K::B>  {        \
    typedef typename boost::variant<typename K::R1, typename K::R2,     \
                                    typename K::R3> variant_type;       \
    typedef typename boost::optional< variant_type > result_type;       \
  };

namespace CGAL {

// only declaration
template<typename, typename, typename>
struct Intersection_traits_2;

CGAL_INTERSECTION_TRAITS_2(Line_2, Line_2, Point_2, Line_2);

CGAL_INTERSECTION_TRAITS_2(Segment_2, Line_2, Point_2, Segment_2);
CGAL_INTERSECTION_TRAITS_2(Line_2, Segment_2, Point_2, Segment_2);

CGAL_INTERSECTION_TRAITS_2(Segment_2, Segment_2, Point_2, Segment_2);

CGAL_INTERSECTION_TRAITS_2(Ray_2, Line_2, Point_2, Ray_2);
CGAL_INTERSECTION_TRAITS_2(Line_2, Ray_2, Point_2, Ray_2);

CGAL_INTERSECTION_TRAITS_2(Ray_2, Segment_2, Point_2, Segment_2);
CGAL_INTERSECTION_TRAITS_2(Segment_2, Ray_2, Point_2, Segment_2);

CGAL_INTERSECTION_TRAITS_3(Ray_2, Ray_2, Point_2, Segment_2, Ray_2);

CGAL_INTERSECTION_TRAITS_2(Triangle_2, Line_2, Point_2, Segment_2);
CGAL_INTERSECTION_TRAITS_2(Line_2, Triangle_2, Point_2, Segment_2);

CGAL_INTERSECTION_TRAITS_2(Triangle_2, Segment_2, Point_2, Segment_2);
CGAL_INTERSECTION_TRAITS_2(Segment_2, Triangle_2, Point_2, Segment_2);

CGAL_INTERSECTION_TRAITS_2(Triangle_2, Ray_2, Point_2, Segment_2);
CGAL_INTERSECTION_TRAITS_2(Ray_2, Triangle_2, Point_2, Segment_2);

template<typename K>
struct Intersection_traits_2<K, typename K::Triangle_2, typename K::Triangle_2>  {
  typedef typename 
  boost::variant< typename K::Point_2, typename K::Segment_2,
                  typename K::Triangle_2, typename std::vector< typename K::Point_2 > > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

CGAL_INTERSECTION_TRAITS_2(Iso_rectangle_2, Line_2, Point_2, Segment_2);
CGAL_INTERSECTION_TRAITS_2(Line_2, Iso_rectangle_2, Point_2, Segment_2);

CGAL_INTERSECTION_TRAITS_2(Iso_rectangle_2, Segment_2, Point_2, Segment_2);
CGAL_INTERSECTION_TRAITS_2(Segment_2, Iso_rectangle_2, Point_2, Segment_2);

CGAL_INTERSECTION_TRAITS_2(Iso_rectangle_2, Ray_2, Point_2, Segment_2);
CGAL_INTERSECTION_TRAITS_2(Ray_2, Iso_rectangle_2, Point_2, Segment_2);


// Variants of one for backwards compatibility
template<typename K>
struct Intersection_traits_2<K, typename K::Iso_rectangle_2, typename K::Iso_rectangle_2>  {
  typedef typename boost::variant<typename K::Iso_rectangle_2> variant_type;
  typedef boost::optional<variant_type> result_type;
};

template<typename K, typename B>
struct Intersection_traits_2<K, typename K::Point_2, B> {
  typedef typename boost::variant<typename K::Point_2> variant_type;
  typedef boost::optional<variant_type> result_type;
};

template<typename K, typename A>
struct Intersection_traits_2<K, A, typename K::Point_2> {
  typedef typename boost::variant<typename K::Point_2> variant_type;
  typedef boost::optional<variant_type> result_type;
};

}

#endif /* CGAL_INTERSECTION_TRAITS_H */

