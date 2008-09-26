// Copyright (c) 2003, 2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>


//CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Timer.h>

typedef double Coord_type;
typedef CGAL::Exact_predicates_inexact_constructions_kernel	    Rep;

typedef CGAL::Point_2<Rep>                  Point_2;
typedef CGAL::Segment_2<Rep>                Segment;
typedef CGAL::Line_2<Rep>                   Line;
typedef CGAL::Triangle_2<Rep>               Triangle;
typedef CGAL::Circle_2<Rep>                 Circle;
typedef CGAL::Iso_rectangle_2<Rep>          Iso_rectangle_2;

typedef CGAL::Triangulation_2<Rep>          Triangulation;
typedef CGAL::Delaunay_triangulation_2<Rep> Delaunay;


typedef Delaunay::Vertex_iterator           Vertex_iterator;
typedef Delaunay::Face_handle               Face_handle;
typedef Delaunay::Vertex_handle             Vertex_handle;
typedef Delaunay::Edge                      Edge;
typedef Delaunay::Line_face_circulator      Line_face_circulator;
