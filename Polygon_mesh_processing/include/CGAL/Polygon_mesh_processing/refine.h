#ifndef CGAL_POLYGON_MESH_PROCESSING_REFINE_H
#define CGAL_POLYGON_MESH_PROCESSING_REFINE_H

#include <CGAL/Polygon_mesh_processing/internal/refine_impl.h>

namespace CGAL {

namespace Polygon_mesh_processing {

  /*!
  \ingroup PkgPolygonMeshProcessing
  @brief Function refining a region on a polygon mesh

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
  @param density_control_factor factor for density where larger values cause denser refinements

  @return pair of @a faces_out and @a vertices_out

  \todo SUBMISSION: missing VertexPointMap
  @todo current algorithm iterates 10 times at most, since (I guess) there is no termination proof.
  */
  template<class PolygonMesh,
           class FaceRange,
           class FaceOutputIterator,
           class VertexOutputIterator>
  std::pair<FaceOutputIterator, VertexOutputIterator>
    refine(PolygonMesh& pmesh,
           FaceRange faces,
           FaceOutputIterator faces_out,
           VertexOutputIterator vertices_out,
           double density_control_factor = std::sqrt(2.0))
  {
    internal::Refine_Polyhedron_3<PolygonMesh> refine_functor(pmesh);
    refine_functor.refine(faces,
                          faces_out,
                          vertices_out,
                          density_control_factor);
    return std::make_pair(faces_out, vertices_out);
  }

}//end namespace Polygon_mesh_processing

}//end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_REFINE_H
