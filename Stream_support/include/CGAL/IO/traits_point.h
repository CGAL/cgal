// Copyright (c) 2018  GeometryFactory Sarl (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_IO_TRAITS_POINT_H
#define CGAL_IO_TRAITS_POINT_H
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <CGAL/number_utils.h>
#include <CGAL/Point_2.h>
#include <boost/geometry/io/wkt/write.hpp>
#include <boost/geometry/io/wkt/read.hpp>

namespace boost{
namespace geometry{
namespace traits{

//WKT traits for Points
template< typename K > struct tag<CGAL::Point_2<K> >
{ typedef point_tag type; };

template< typename K > struct coordinate_type<CGAL::Point_2<K> >
{ typedef double type; };

template< typename K > struct coordinate_system<CGAL::Point_2<K> >
{ typedef cs::cartesian type; };

template< typename K > struct dimension<CGAL::Point_2<K> > : boost::mpl::int_<2> {};

template< typename K >
struct access<CGAL::Point_2<K> , 0>
{
  static double get(CGAL::Point_2<K>  const& p)
  {
    return CGAL::to_double(p.x());
  }
  
  static void set(CGAL::Point_2<K> & p, typename K::FT c)
  {
    p = CGAL::Point_2<K> (c, p.y());
  }
  
};

template< typename K >
struct access<CGAL::Point_2<K> , 1>
{
  static double get(CGAL::Point_2<K>  const& p)
  {
    return CGAL::to_double(p.y());
  }
  
  static void set(CGAL::Point_2<K> & p, typename K::FT c)
  {
    p = CGAL::Point_2<K> (p.x(), c);
  }
  
};

}}}//end namespaces
#endif
#endif
