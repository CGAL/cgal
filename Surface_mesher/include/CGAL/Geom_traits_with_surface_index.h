// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Laurent RINEAU

#ifndef GEOM_TRAITS_WITH_SURFACE_INDEX_H
#define GEOM_TRAITS_WITH_SURFACE_INDEX_H_H

#include <CGAL/Point_with_surface_index.h>
#include <CGAL/Segment_with_surface_index.h>
#include <CGAL/Triangle_with_surface_index.h>

namespace CGAL {

template <class GT>
class Geom_traits_with_surface_index : public GT
{
  typedef typename GT::Point_3 Old_point_3;
  typedef typename GT::Segment_3 Old_segment_3;
  typedef typename GT::Triangle_3 Old_triangle_3;

public:
  typedef Point_with_surface_index<Old_point_3> Point_3;
  typedef Segment_with_surface_index<Old_triangle_3> Segment_3;
  typedef Triangle_with_surface_index<Old_triangle_3> Triangle_3;

};  // end Geom_traits_with_surface_index

} // end namespace CGAL

#endif // GEOM_TRAITS_WITH_SURFACE_INDEX_H
