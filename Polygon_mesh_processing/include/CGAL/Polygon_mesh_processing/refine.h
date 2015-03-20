// Copyright (c) 2013 GeometryFactory (France).
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
// Author(s)     : Ilker O. Yaz

#ifndef CGAL_POLYGON_MESH_PROCESSING_REFINE_H
#define CGAL_POLYGON_MESH_PROCESSING_REFINE_H

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>

#include <CGAL/Polygon_mesh_processing/internal/refine_impl.h>

namespace CGAL {

namespace Polygon_mesh_processing {

  /*!
  \ingroup PkgPolygonMeshProcessing
  @brief refines a region of a polygon mesh

  @tparam PolygonMesh model of `MutableFaceGraph`
  @tparam FaceRange range of face descriptors, model of `SinglePassRange`
  @tparam FaceOutputIterator model of `OutputIterator`
    holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces
  @tparam VertexOutputIterator model of `OutputIterator`
    holding `boost::graph_traits<PolygonMesh>::%vertex_descriptor` for patch vertices

  @param pmesh polygon mesh with patches to be refined
  @param faces the range of faces defining the patches to refine
  @param faces_out output iterator into which descriptors of new faces are put
  @param vertices_out output iterator into which descriptors of new vertices are put

  The function accepts named parameters
   - `CGAL::parameters::density_control_factor`  which defaults to `CGAL::sqrt(2.)`

  @param density_control_factor factor for density where larger values cause denser refinements

  @return pair of `faces_out` and `vertices_out`

  \todo SUBMISSION: missing VertexPointMap
  \todo SUBMISSION: better document density_control_factor 
  @todo current algorithm iterates 10 times at most, since (I guess) there is no termination proof.
  */
  template<class PolygonMesh,
           class FaceRange,
           class FaceOutputIterator,
           class VertexOutputIterator,
           class P, class T, class R>
  std::pair<FaceOutputIterator, VertexOutputIterator>
    refine(PolygonMesh& pmesh,
           FaceRange faces,
           FaceOutputIterator faces_out,
           VertexOutputIterator vertices_out,
           const pmp_bgl_named_params<P, T, R>& p)
  {
    using boost::choose_param;
    using boost::get_param;

    internal::Refine_Polyhedron_3<PolygonMesh> refine_functor(pmesh);
    refine_functor.refine(faces,
      faces_out,
      vertices_out,
      choose_param(get_param(p, density_control_factor), CGAL::sqrt(2.)));
    return std::make_pair(faces_out, vertices_out);
  }

}//end namespace Polygon_mesh_processing

}//end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_REFINE_H
