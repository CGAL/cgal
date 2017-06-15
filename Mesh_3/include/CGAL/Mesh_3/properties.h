// Copyright (c) 2017 GeometryFactory (France).
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_POLYGON_MESH_PROCESSING_PROPERTIES_H
#define CGAL_POLYGON_MESH_PROCESSING_PROPERTIES_H

#include <CGAL/license/Polygon_mesh_processing.h>

namespace CGAL
{

enum vertex_num_feature_edges_t { vertex_num_feature_edges };
enum halfedge_is_feature_t      { halfedge_is_feature };

enum vertex_time_stamp_t        { vertex_time_stamp};
enum halfedge_time_stamp_t      { halfedge_time_stamp};
enum face_time_stamp_t          { face_time_stamp};

template <typename ID>
struct vertex_incident_patches_t {
  typedef ID type;
};

template <typename ID>
struct face_patch_id_t {
  typedef ID type;
};
} //end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_PROPERTIES_H
