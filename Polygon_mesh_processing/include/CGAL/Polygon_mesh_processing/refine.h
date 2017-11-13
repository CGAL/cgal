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

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>


#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/Polygon_mesh_processing/internal/refine_impl.h>

namespace CGAL {

namespace Polygon_mesh_processing {

  /*!
  \ingroup PMP_meshing_grp
  @brief refines a region of a triangle mesh

  @tparam TriangleMesh model of `MutableFaceGraph`
  @tparam FaceRange range of face descriptors, model of `Range`.
          Its iterator type is `InputIterator`.
  @tparam FaceOutputIterator model of `OutputIterator`
    holding `boost::graph_traits<TriangleMesh>::%face_descriptor` for patch faces
  @tparam VertexOutputIterator model of `OutputIterator`
    holding `boost::graph_traits<TriangleMesh>::%vertex_descriptor` for patch vertices
  @tparam NamedParameters a sequence of \ref namedparameters

  @param tmesh triangle mesh with patches to be refined
  @param faces the range of faces defining the patches to refine
  @param faces_out output iterator into which descriptors of new faces are recorded
  @param vertices_out output iterator into which descriptors of new vertices are recorded
  @param np optional sequence of \ref namedparameters among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `tmesh`
      Instance of a class model of `ReadWritePropertyMap`.
      If this parameter is omitted, an internal property map for
     `CGAL::vertex_point_t` should be available in `TriangleMesh`
     \cgalParamEnd
    \cgalParamBegin{density_control_factor} factor to control density of the output mesh,
      where larger values lead to denser refinements.
      The density of vertices of `faces_out` is this factor times higher than the vertices of `faces.` \cgalParamEnd
  \cgalNamedParamsEnd


  @return pair of `faces_out` and `vertices_out`

  \pre `is_triangle_mesh(tmesh)`

  @todo current algorithm iterates 10 times at most, since (I guess) there is no termination proof.
  */
  template<typename TriangleMesh,
           typename FaceRange,
           typename FaceOutputIterator,
           typename VertexOutputIterator,
           typename NamedParameters>
  std::pair<FaceOutputIterator, VertexOutputIterator>
    refine(TriangleMesh& tmesh,
           const FaceRange& faces,
           FaceOutputIterator faces_out,
           VertexOutputIterator vertices_out,
           const NamedParameters& np)
  {
    using boost::choose_param;
    using boost::get_param;

    CGAL_precondition(is_triangle_mesh(tmesh) );

    typedef typename GetVertexPointMap<TriangleMesh,NamedParameters>::type VPmap;
    VPmap vpm = choose_param(get_param(np, internal_np::vertex_point),
                             get_property_map(vertex_point, tmesh));

    internal::Refine_Polyhedron_3<TriangleMesh, VPmap> refine_functor(tmesh, vpm);
    refine_functor.refine(faces,
      faces_out,
      vertices_out,
      choose_param(get_param(np, internal_np::density_control_factor), CGAL::sqrt(2.)));
    return std::make_pair(faces_out, vertices_out);
  }

///\cond SKIP_IN_MANUAL
  template<typename TriangleMesh,
    typename FaceRange,
    typename FaceOutputIterator,
    typename VertexOutputIterator>

  std::pair<FaceOutputIterator, VertexOutputIterator>
    refine(TriangleMesh& tmesh,
           const FaceRange& faces,
           FaceOutputIterator faces_out,
           VertexOutputIterator vertices_out)
  {
    return refine(tmesh, faces, faces_out, vertices_out,
      CGAL::Polygon_mesh_processing::parameters::all_default());
  }
///\endcond
}//end namespace Polygon_mesh_processing

}//end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_REFINE_H
