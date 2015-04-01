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
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/Polygon_mesh_processing/internal/refine_impl.h>

namespace CGAL {

namespace Polygon_mesh_processing {

  /*!
  \ingroup PkgPolygonMeshProcessing
  @brief refines a region of a polygon mesh

  @tparam PolygonMesh model of `MutableFaceGraph`
          that has an internal property map for `CGAL::vertex_point_t`
  @tparam FaceRange range of face descriptors, model of `SinglePassRange`
  @tparam FaceOutputIterator model of `OutputIterator`
    holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces
  @tparam VertexOutputIterator model of `OutputIterator`
    holding `boost::graph_traits<PolygonMesh>::%vertex_descriptor` for patch vertices
  @tparam NamedParameters a sequence of \ref namedparameters

  @param pmesh polygon mesh with patches to be refined
  @param faces the range of faces defining the patches to refine
  @param faces_out output iterator into which descriptors of new faces are put
  @param vertices_out output iterator into which descriptors of new vertices are put
  @param np optional sequence of \ref namedparameters among the ones listed below

  \b Named \b parameters
  <ul>
  <li>\b vertex_point_map the property map with the points associated to the vertices of `pmesh`.
  Instance of a class model of `ReadWritePropertyMap`
  <li>\b density_control_factor factor to control density of the ouput mesh, where larger values cause denser refinements. The density of vertices of `faces_out` is this factor times higher than the vertices of `faces.`
  </ul>

  @return pair of `faces_out` and `vertices_out`

  \todo code: in refine_impl, is density_control_factor used? it seems that scale_attribute
               is filled with zeros
  @todo current algorithm iterates 10 times at most, since (I guess) there is no termination proof.
  */
  template<typename PolygonMesh,
           typename FaceRange,
           typename FaceOutputIterator,
           typename VertexOutputIterator,
           typename NamedParameters>
  std::pair<FaceOutputIterator, VertexOutputIterator>
    refine(PolygonMesh& pmesh,
           FaceRange faces,
           FaceOutputIterator faces_out,
           VertexOutputIterator vertices_out,
           const NamedParameters& np)
  {
    using boost::choose_pmap;
    using boost::choose_param;
    using boost::get_param;

    typedef typename GetVertexPointMap<PolygonMesh,NamedParameters>::type VPmap;
    VPmap vpm = choose_pmap(get_param(np, boost::vertex_point),
                            pmesh,
                            boost::vertex_point);

    internal::Refine_Polyhedron_3<PolygonMesh, VPmap> refine_functor(pmesh, vpm);
    refine_functor.refine(faces,
      faces_out,
      vertices_out,
      choose_param(get_param(np, density_control_factor), CGAL::sqrt(2.)));
    return std::make_pair(faces_out, vertices_out);
  }

}//end namespace Polygon_mesh_processing

}//end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_REFINE_H
