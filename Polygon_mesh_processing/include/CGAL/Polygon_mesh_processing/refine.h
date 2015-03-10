#ifndef CGAL_POLYGON_MESH_PROCESSING_REFINE_H
#define CGAL_POLYGON_MESH_PROCESSING_REFINE_H

#include <CGAL/internal/Meshing_functions/Refine_Polyhedron_3.h>

namespace CGAL {

namespace Polygon_mesh_processing {

  /*!
  \ingroup PkgPolygonMeshProcessing
  @brief Function refining a region on polygon mesh

  @tparam PolygonMesh must be a model of `MutableFaceGraph`
  @tparam FacetRange range of input facets, model of `SinglePassRange`
  @tparam FacetOutputIterator iterator holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch facets
  @tparam VertexOutputIterator iterator holding `boost::graph_traits<PolygonMesh>::%vertex_descriptor` for patch vertices

  @param pmesh mesh to be refined
  @param facets the range of facets to be refined
  @param facet_out iterator over newly created facets
  @param vertex_out iterator over newly created vertices
  @param density_control_factor factor for density where larger values cause denser refinements

  @return pair of @a facet_out and @a vertex_out

  @todo current algorithm iterates 10 times at most, since (I guess) there is no termination proof.
  */
  template<class PolygonMesh,
           class FacetRange,
           class FacetOutputIterator,
           class VertexOutputIterator>
  std::pair<FacetOutputIterator, VertexOutputIterator>
    refine(PolygonMesh& pmesh,
           FacetRange facets,
           FacetOutputIterator facet_out,
           VertexOutputIterator vertex_out,
           double density_control_factor = std::sqrt(2.0))
  {
    internal::Refine_Polyhedron_3<PolygonMesh> refine_functor(pmesh);
    refine_functor.refine(facets,
                          facet_out,
                          vertex_out,
                          density_control_factor);
    return std::make_pair(facet_out, vertex_out);
  }

}//end namespace Polygon_mesh_processing

}//end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_REFINE_H
