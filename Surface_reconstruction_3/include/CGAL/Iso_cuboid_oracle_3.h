// Copyright (c) 2006-2008  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU, Laurent Saboret

#ifndef CGAL_ISO_CUBOID_ORACLE_3_H
#define CGAL_ISO_CUBOID_ORACLE_3_H

#include <CGAL/Surface_mesh_traits_generator_3.h>

CGAL_BEGIN_NAMESPACE


/// Mesher level oracle which tests if a point belongs to a 3D cube.
template <class GT>
class Iso_cuboid_oracle_3
{
public:

  // Public types
  typedef GT Geom_traits;
  typedef typename GT::Point_3 Point;
  typedef typename GT::Iso_cuboid_3 Iso_cuboid_3;

  typedef Iso_cuboid_3 Surface_3;

  typedef Point Intersection_point;

public:

  // Constructors
  Iso_cuboid_oracle_3 ()
  {
  }

  // Predicates 
  bool is_in_volume(const Surface_3& cube, const Point& p) const
  {
    typename GT::Has_on_bounded_side_3 on_bounded_side_of_cube =
                                    GT().has_on_bounded_side_3_object();

    return on_bounded_side_of_cube(cube, p);
  }

};  // end Iso_cuboid_oracle_3


/// Specialization of Surface_mesh_traits_generator_3 for Iso_cuboid_3.
template <typename Kernel>
struct Surface_mesh_traits_generator_3<CGAL::Iso_cuboid_3<Kernel> >
{
  typedef Iso_cuboid_oracle_3<Kernel> Type;
  typedef Type type; // for Boost compatiblity (meta-programming)
};


CGAL_END_NAMESPACE

#endif  // CGAL_ISO_CUBOID_ORACLE_3_H
