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
#if BOOST_VERSION >= 105600 || BOOST_GCC < 40500

#ifndef CGAL_IO_TRAITS_LINESTRING_H
#define CGAL_IO_TRAITS_LINESTRING_H
#include <CGAL/internal/Geometry_container.h>
#include <boost/geometry/io/wkt/write.hpp>
#include <boost/geometry/io/wkt/read.hpp>



namespace boost{
namespace geometry{
namespace traits{
template< typename R> struct tag<CGAL::internal::Geometry_container<R, linestring_tag> >
{ typedef linestring_tag type; };

}}} //end namespaces

#endif // TRAITS_LINESTRING_H
#endif
