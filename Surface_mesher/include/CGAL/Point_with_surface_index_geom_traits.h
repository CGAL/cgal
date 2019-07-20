// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_POINT_WITH_SURFACE_INDEX_GEOM_TRAITS_H
#define CGAL_POINT_WITH_SURFACE_INDEX_GEOM_TRAITS_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/Point_with_surface_index.h>

namespace CGAL {

template <class GT>
class Point_with_surface_index_geom_traits : public GT
{
  typedef typename GT::Point_3 Old_point_3;

public:
  typedef Point_with_surface_index<Old_point_3> Point_3;

};  // end Point_with_surface_index_geom_traits

} // end namespace CGAL

#endif // CGAL_POINT_WITH_SURFACE_INDEX_GEOM_TRAITS_H
