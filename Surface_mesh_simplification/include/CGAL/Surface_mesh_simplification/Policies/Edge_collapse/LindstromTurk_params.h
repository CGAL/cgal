// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_PARAMS_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_PARAMS_H

#include <CGAL/license/Surface_mesh_simplification.h>


namespace CGAL {

namespace Surface_mesh_simplification 
{

struct LindstromTurk_params
{
  LindstromTurk_params()
    :
    VolumeWeight  (0.5)
   ,BoundaryWeight(0.5)
   ,ShapeWeight   (0)
  {}
  
  LindstromTurk_params( double aVolumeWeight, double aBoundaryWeight, double aShapeWeight )
    :
    VolumeWeight  (aVolumeWeight)
   ,BoundaryWeight(aBoundaryWeight)
   ,ShapeWeight   (aShapeWeight)
  {}
    
  double VolumeWeight ;
  double BoundaryWeight ;
  double ShapeWeight ;
};
  
} // namespace Surface_mesh_simplification

} //namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_PARAMS_H
// EOF //
 
