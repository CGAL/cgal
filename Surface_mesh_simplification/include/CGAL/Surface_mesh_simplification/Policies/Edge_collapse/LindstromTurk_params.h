// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Surface_mesh_simplification/include/CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_collapse_data.h $
// $Id: LindstromTurk_collapse_data.h 33643 2006-08-24 13:30:51Z fcacciola $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_PARAMS_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_PARAMS_H

#include <CGAL/Surface_mesh_simplification/TSMS_common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Minimal_collapse_data.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification { namespace Edge_collapse 
{

struct LinstromTurk_params
{
  LinstromTurk_params()
    :
    VolumeWeight  (0.5)
   ,BoundaryWeight(0.5)
   ,ShapeWeight   (0)
  {}
  
  LinstromTurk_params( double aVolumeWeight, double aBoundaryWeight, double aShapeWeight )
    :
    VolumeWeight  (aVolumeWeight)
   ,BoundaryWeight(aBoundaryWeight)
   ,ShapeWeight   (aShapeWeight)
  {}
    
  double VolumeWeight ;
  double BoundaryWeight ;
  double ShapeWeight ;
};
  
} } } // namespace Triangulated_surface_mesh::Simplification::Edge_collapse

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_PARAMS_H
// EOF //
 
