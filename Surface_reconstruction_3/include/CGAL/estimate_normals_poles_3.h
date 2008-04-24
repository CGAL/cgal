// Copyright (c) 2007-08  INRIA Sophia-Antipolis (France).
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
// Author(s) : Pierre Alliez and Laurent Saboret

#ifndef CGAL_ESTIMATE_NORMALS_POLES_3_H
#define CGAL_ESTIMATE_NORMALS_POLES_3_H

#include <CGAL/basic.h>

#include <iterator>
#include <list>

CGAL_BEGIN_NAMESPACE

// Estimate normal directions using poles of a 3D Voronoi diagram
// 8 vertices of a loose bounding box are added in order to bound
// all Voronoi cells. The last parameter specifies if these points
// should be removed.
// return past-the-end iterator of output
template < typename InputIterator, 
           typename OutputIterator,
           typename Kernel, 
           typename DT> // Delaunay triangulation
OutputIterator
estimate_normals_poles_3(InputIterator first,      // input points
                         InputIterator beyond,   
                         OutputIterator normals,   // output normals
                         DT& dt,                   // Delaunay triangulation
                         const Kernel::FT,         // inflating factor of bounding box
                         const bool remove_bounding_points)
{
  // basic kernel object types
  typedef typename Kernel::FT       FT;
  typedef typename Kernel::Point_3  Point;
  typedef typename Kernel::Vector_3 Vector;

  // precondition: at least one element in the container.
  CGAL_surface_reconstruction_precondition(first != beyond);
}

CGAL_END_NAMESPACE

#endif // CGAL_ESTIMATE_NORMALS_POLES_3_H

