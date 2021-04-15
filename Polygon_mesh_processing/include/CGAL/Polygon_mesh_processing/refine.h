// Copyright (c) 2013 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
  @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  @param tmesh triangle mesh with patches to be refined
  @param faces the range of faces defining the patches to refine
  @param faces_out output iterator into which descriptors of new faces are recorded
  @param vertices_out output iterator into which descriptors of new vertices are recorded
  @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `TriangleMesh`.}
    \cgalParamNEnd

    \cgalParamNBegin{density_control_factor}
      \cgalParamDescription{a factor to control density of the output mesh, where larger values lead to denser refinements}
      \cgalParamType{double}
      \cgalParamDefault{\f$ \sqrt{2}\f$}
      \cgalParamExtra{The density of vertices of `faces_out` is this factor times higher than the vertices of `faces`.}
    \cgalParamNEnd
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
    using parameters::choose_parameter;
    using parameters::get_parameter;

    CGAL_precondition(is_triangle_mesh(tmesh) );

    typedef typename GetVertexPointMap<TriangleMesh,NamedParameters>::type VPmap;
    VPmap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(vertex_point, tmesh));

    internal::Refine_Polyhedron_3<TriangleMesh, VPmap> refine_functor(tmesh, vpm);
    refine_functor.refine(faces,
      faces_out,
      vertices_out,
      choose_parameter(get_parameter(np, internal_np::density_control_factor), CGAL::sqrt(2.)));
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
