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

#ifndef CGAL_IO_WKT_TRAITS_LINESTRING_H
#define CGAL_IO_WKT_TRAITS_LINESTRING_H

#include <CGAL/Stream_support/internal/Geometry_container.h>

#include <boost/geometry/io/wkt/write.hpp>
#include <boost/geometry/io/wkt/read.hpp>

namespace boost {
namespace geometry {
namespace traits {

template< typename R>
struct tag<CGAL::internal::Geometry_container<R, linestring_tag> >
{
  typedef linestring_tag type;
};

} // namespace traits
} // namespace geometry
} // namespace boost

#endif // CGAL_IO_WKT_TRAITS_LINESTRING_H
