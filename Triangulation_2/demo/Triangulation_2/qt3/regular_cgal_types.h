// Copyright (c) 2003, 2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>


//CGAL headers
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>

typedef CGAL::Simple_cartesian<double> Rp;
typedef Rp::Point_2 Point_2;
typedef Rp::Circle_2 Circle;
typedef double W;
typedef CGAL::Regular_triangulation_euclidean_traits_2<Rp,W>  Gt;
// typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
// typedef CGAL::Regular_triangulation_face_base_2<> Fb;
// typedef CGAL::Triangulation_data_structure_2<Vb,Fb > Tds;
typedef CGAL::Regular_triangulation_2<Gt> Regular_triangulation;
typedef Regular_triangulation::Finite_vertices_iterator
                                          Finite_vertices_iterator;
