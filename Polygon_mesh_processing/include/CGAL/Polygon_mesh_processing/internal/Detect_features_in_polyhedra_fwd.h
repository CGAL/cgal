// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYHEDRA_FWD_H
#define CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYHEDRA_FWD_H

#include <CGAL/license/Polygon_mesh_processing.h>


namespace CGAL {
namespace Polygon_mesh_processing {
  
  template <typename Polyhedron, typename PatchId_pmap>
  class Detect_features_in_polyhedra;

  template <typename Polyhedron, typename FT, typename PatchId_pmap>
  void detect_features(Polyhedron& p,
                       FT angle_in_deg,
                       PatchId_pmap pid_map);
  
} // end namespace Polygon_mesh_processing
  
} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYHEDRA_FWD_H
