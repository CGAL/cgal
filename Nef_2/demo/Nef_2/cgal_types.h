// Copyright (c) 2002  Max-Planck-Institute Saarbruecken (Germany)
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
//
//
// Author(s)     : Radu Ursu

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Filtered_extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_2.h>



//Cartesian types
typedef CGAL::Cartesian<CGAL::Gmpq>   Rp;
typedef CGAL::Partition_traits_2<Rp>
                                  Traits;
typedef Rp::Point_2               Cartesian_point_2;
typedef Rp::Line_2                Cartesian_line_2;
typedef Traits::Polygon_2         Cartesian_polygon_2;
typedef Cartesian_polygon_2::Vertex_iterator
                                  Vertex_iterator;

typedef CGAL::Cartesian<double>   RP_double;
typedef CGAL::Partition_traits_2<RP_double>
                                  Traits_double;
typedef Traits_double::Polygon_2  Polygon_2_double;
typedef Polygon_2_double::Vertex_iterator
                                  Vertex_iterator_double;

//The Nef_Polyhedron types
typedef CGAL::Gmpz                RT;
typedef CGAL::Filtered_extended_homogeneous<RT>
                                  Extended_kernel;
typedef CGAL::Nef_polyhedron_2<Extended_kernel>
                                  Nef_polyhedron;
typedef Nef_polyhedron::Point     Point_2;
typedef Nef_polyhedron::Line      Line;
