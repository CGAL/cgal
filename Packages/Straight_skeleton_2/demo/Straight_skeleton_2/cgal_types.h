// Copyright (c) 2002  Max Planck Institut fuer Informatik (Germany).
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

#include <boost/shared_ptr.hpp>
#include <vector>

typedef double                     NT;
typedef CGAL::Cartesian<NT>        K;
typedef K::Point_2                 Point;
typedef CGAL::Polygon_2<K>         Polygon;
typedef boost::shared_ptr<Polygon> PolygonPtr;
typedef CGAL::Segment_2<K>         Segment;
typedef std::vector<PolygonPtr>    PolygonalRegion ;

