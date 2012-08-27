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

#ifndef CGAL_MESH_3_DETECT_FEATURES_IN_POLYHEDRA_FWD_H
#define CGAL_MESH_3_DETECT_FEATURES_IN_POLYHEDRA_FWD_H

namespace CGAL {
namespace Mesh_3 {
  
  template <typename Polyhedron>
  class Detect_features_in_polyhedra;

  template <typename Polyhedron>
  void detect_features(Polyhedron& p);
  
} // end namespace Mesh_3
  
} // end namespace CGAL

#endif // CGAL_MESH_3_DETECT_FEATURES_IN_POLYHEDRA_FWD_H
