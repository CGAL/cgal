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

#ifndef CGAL_IO_TRAITS_MULTIPOINT_H
#define CGAL_IO_TRAITS_MULTIPOINT_H
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <CGAL/internal/Geometry_container.h>
#include <boost/geometry/io/wkt/write.hpp>
#include <boost/geometry/io/wkt/read.hpp>

namespace boost{
namespace geometry{
namespace traits{
// WKT traits for MultiPoint
template< typename R >
struct tag<CGAL::internal::Geometry_container<R, multi_point_tag > >
{ typedef multi_point_tag type; };

}//end traits
}//end geometry
}//end boost
#endif // TRAITS_MULTIPOINT_H
#endif
