// Copyright (c) 1999  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
//                 Monique Teillaud <Monique.Teillaud@sophia.inria.fr>

// geometric traits for a <=3 D triangulation

// OBSOLETE !!!!!!!!!!!!!!!!!

#ifndef CGAL_TRIANGULATION_GEOM_TRAITS_3_H
#define CGAL_TRIANGULATION_GEOM_TRAITS_3_H

#include <CGAL/license/Triangulation_3.h>


#include <CGAL/basic.h>

#include <CGAL/triangulation_assertions.h>

namespace CGAL {

template < class Repres >
class Triangulation_geom_traits_3 : public Repres
{
public:
  Triangulation_geom_traits_3()
  {
    bool The_class_Triangulation_geom_traits_3_is_obsolete;
  }

  typedef Repres Rep;

  typedef typename Rep::Object_3       Object_3;
  typedef typename Rep::Point_3        Point_3;
  typedef typename Rep::Segment_3      Segment_3;
  typedef typename Rep::Triangle_3     Triangle_3;
  typedef typename Rep::Tetrahedron_3  Tetrahedron_3;
  typedef typename Rep::Ray_3          Ray_3;
  typedef typename Rep::Line_3         Line_3;
  typedef typename Rep::Plane_3        Plane_3;

  // The next typedef is there for backward compatibility
  // Some users take their point type from the traits class.
  // Before this type was Point
  typedef Point_3                      Point;

  typedef typename Rep::Compare_x_3                Compare_x_3;
  typedef typename Rep::Compare_y_3                Compare_y_3;
  typedef typename Rep::Compare_z_3                Compare_z_3;
  typedef typename Rep::Equal_3                    Equal_3;
  typedef typename Rep::Orientation_3              Orientation_3;
  typedef typename Rep::Coplanar_orientation_3     Coplanar_orientation_3;
  typedef typename Rep::Side_of_oriented_sphere_3  Side_of_oriented_sphere_3;
  typedef typename Rep::Coplanar_side_of_bounded_circle_3
                                          Coplanar_side_of_bounded_circle_3;

  typedef typename Rep::Construct_segment_3        Construct_segment_3;
  typedef typename Rep::Construct_triangle_3       Construct_triangle_3;
  typedef typename Rep::Construct_tetrahedron_3    Construct_tetrahedron_3;
  typedef typename Rep::Construct_object_3         Construct_object_3;

  // For the hierarchy :
  typedef typename Rep::Compare_distance_3         Compare_distance_3;
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_GEOM_TRAITS_3_H
