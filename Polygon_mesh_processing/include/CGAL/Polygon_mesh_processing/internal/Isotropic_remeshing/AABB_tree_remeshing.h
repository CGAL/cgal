// Copyright (c) 2015 GeometryFactory (France).
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
// Author(s)     : Jane Tournois

#ifndef CGAL_POLYGON_MESH_PROCESSING_AABB_TREE_REMESHING_H
#define CGAL_POLYGON_MESH_PROCESSING_AABB_TREE_REMESHING_H

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/AABB_filtered_projection_traits.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>

namespace PMP = CGAL::Polygon_mesh_processing;

template<typename PolygonMesh
       , typename VertexPointMap>
class AABB_facet_primitive_remeshing
  : public CGAL::AABB_face_graph_triangle_primitive<PolygonMesh, VertexPointMap>
{
  typedef CGAL::AABB_face_graph_triangle_primitive<PolygonMesh, VertexPointMap> Base;
  typedef typename boost::graph_traits<PolygonMesh>::face_iterator face_iterator;

public:
  AABB_facet_primitive_remeshing(face_iterator it,
                                 const PolygonMesh &pmesh,
                                 VertexPointMap vpm)
    : Base(it, pmesh, vpm)
  {}

  void set_patch_id(const std::size_t& patch_id) {
    patch_id_ = patch_id;
  }
  std::size_t patch_id() const {
    return patch_id_;
  }

private:
  std::size_t patch_id_;
};

template<typename PolygonMesh
       , typename VertexPointMap
       , typename GeomTraits
>
class AABB_tree_remeshing
  : public CGAL::AABB_tree<
      CGAL::AABB_traits<
        GeomTraits,
        AABB_facet_primitive_remeshing<PolygonMesh, VertexPointMap> > >
{
  typedef PolygonMesh    PM;
  typedef VertexPointMap VPMap;

  typedef typename boost::graph_traits<PM>::face_iterator face_iterator;
  typedef typename boost::graph_traits<PM>::face_descriptor face_descriptor;

public:
  typedef AABB_facet_primitive_remeshing<PM, VPMap>       AABB_primitive;
  typedef CGAL::AABB_traits<GeomTraits, AABB_primitive>   AABB_traits;
  typedef CGAL::AABB_tree<AABB_traits>                    AABB_tree;

public:
  template<typename FaceIndexMap>
  void build(face_iterator fbegin,
             face_iterator fend,
             const PM& pmesh,
             VertexPointMap vpmap,
             FaceIndexMap fccmap)//connected components
  {
    for (face_iterator fit = fbegin; fit != fend; ++fit)
    {
      face_descriptor f = *fit;
      AABB_primitive prim(f, pmesh, vpmap);
      prim.set_patch_id(get(fccmap, f));
      this->insert(prim);
    }
  }
};



#endif //CGAL_POLYGON_MESH_PROCESSING_AABB_TREE_REMESHING_H