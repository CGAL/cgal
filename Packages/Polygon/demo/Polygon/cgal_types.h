// Copyright (c) 2002  ETH Zurich (Switzerland).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Radu Ursu

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/point_generators_2.h>


typedef double                                   Coord_type;
typedef CGAL::Cartesian<Coord_type>              K;
typedef K::Point_2                               Point_2;
typedef std::vector<Point_2>                     Container;
typedef CGAL::Polygon_2<K, Container>            Cgal_Polygon;
typedef CGAL::Random_points_in_square_2<Point_2> Point_generator;
