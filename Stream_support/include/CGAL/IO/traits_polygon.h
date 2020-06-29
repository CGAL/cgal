// Copyright (c) 2018  GeometryFactory Sarl (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_IO_TRAITS_POLYGON_H
#define CGAL_IO_TRAITS_POLYGON_H
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Point_2.h>
#include <boost/geometry/io/wkt/write.hpp>
#include <boost/geometry/io/wkt/read.hpp>

namespace boost{
namespace geometry{
namespace traits{
// WKT traits for Polygon
template< typename K > struct tag<CGAL::Polygon_2<K> >
{ typedef ring_tag type; };

template< typename K >
struct tag<CGAL::Polygon_with_holes_2<K> >
{ typedef polygon_tag type; };

template< typename K >
struct ring_const_type<CGAL::Polygon_with_holes_2<K> >
{ typedef const CGAL::Polygon_2<K>& type; };

template< typename K >
struct ring_mutable_type<CGAL::Polygon_with_holes_2<K> >
{ typedef CGAL::Polygon_2<K>& type; };

template< typename K >
struct interior_const_type<CGAL::Polygon_with_holes_2<K> >
{ typedef const typename CGAL::Polygon_with_holes_2<K>::Holes_container& type; };

template< typename K >
struct interior_mutable_type<CGAL::Polygon_with_holes_2<K> >
{ typedef typename CGAL::Polygon_with_holes_2<K>::Holes_container& type; };

template< typename K >
struct exterior_ring<CGAL::Polygon_with_holes_2<K> >
{
  static CGAL::Polygon_2<K>& get(CGAL::Polygon_with_holes_2<K>& p)
  {
    return (p.outer_boundary());
  }
  static CGAL::Polygon_2<K> const& get(CGAL::Polygon_with_holes_2<K> const& p)
  {
    return (p.outer_boundary());
  }
};

template< typename K >
struct interior_rings<CGAL::Polygon_with_holes_2<K> >
{
  static typename CGAL::Polygon_with_holes_2<K>::Holes_container& get(CGAL::Polygon_with_holes_2<K>& p)
  {
    return p.holes();
  }
  static const typename CGAL::Polygon_with_holes_2<K>::Holes_container& get(CGAL::Polygon_with_holes_2<K> const& p)
  {
    return p.holes();
  }
};
}//end traits
}//end geometry

//extra specialization
template< typename K >
struct range_value<CGAL::Polygon_2<K> >
{
  typedef typename CGAL::Polygon_2<K>::Point_2  type;
};

}//end boost

#endif
#endif
