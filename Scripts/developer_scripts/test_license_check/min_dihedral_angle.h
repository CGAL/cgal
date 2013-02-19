// Copyright (c) 2007-2009  INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.1-branch/Mesh_3/include/CGAL/Mesh_3/min_dihedral_angle.h $
// $Id: min_dihedral_angle.h 67117 2012-01-13 18:14:48Z lrineau $
// 
//
// Author(s)     : Laurent RINEAU, Stephane Tayeb

#ifndef CGAL_MESH_3_MIN_DIHEDRAL_ANGLE_H
#define CGAL_MESH_3_MIN_DIHEDRAL_ANGLE_H

#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <cmath>

namespace CGAL {

namespace Mesh_3 {
  
  namespace details {
    
    template <typename K>
    typename K::FT
    min_dihedral_angle_aux_compute_quotient(const typename K::Point_3& p0,
                                            const typename K::Point_3& p1,
                                            const typename K::Point_3& p2,
                                            const typename K::Point_3& p3,
                                            K k = K())
    {
      typename K::Construct_triangle_3 make_triangle = 
      k.construct_triangle_3_object();
      typename K::Compute_area_3 area = 
      k.compute_area_3_object();
      typename K::Compute_squared_distance_3 sq_distance = 
      k.compute_squared_distance_3_object();
      
      return CGAL::sqrt(sq_distance(p0, p1))
        / area(make_triangle(p0, p1, p3))
        / area(make_triangle(p0, p1, p2));
    }
    
  } // end namespace details;
