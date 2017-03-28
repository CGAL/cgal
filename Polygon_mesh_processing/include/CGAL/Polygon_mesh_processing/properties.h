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
enum face_patch_id_t            { face_patch_id };
enum vertex_num_feature_edges_t { vertex_num_feature_edges };
enum halfedge_is_feature_t      { halfedge_is_feature };
enum vertex_selection_t         { vertex_selection};
enum edge_selection_t           { edge_selection};
enum face_selection_t           { face_selection};
} //end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_PROPERTIES_H
