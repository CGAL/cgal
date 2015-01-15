
#include <CGAL/internal/Meshing_functions/Refine_Polyhedron_3.h>

namespace CGAL {

namespace Polygon_mesh_processing {

  /*!
  \ingroup PkgPolygonMeshProcessing
  @brief Function refining a region on polygon mesh

  @tparam Polyhedron must be a model of `MutableFaceGraph`
  @tparam InputIterator iterator over input facets
  @tparam FacetOutputIterator iterator holding `boost::graph_traits<PolygonMesh>::face_descriptor` for patch facets
  @tparam VertexOutputIterator iterator holding `boost::graph_traits<PolygonMesh>::vertex_descriptor` for patch vertices

  @param pmesh mesh to be refined
  @param facet_begin first iterator of the range of facets
  @param facet_end past-the-end iterator of the range of facets
  @param facet_out iterator over newly created facets
  @param vertex_out iterator over newly created vertices
  @param density_control_factor factor for density where larger values cause denser refinements

  @return pair of @a facet_out and @a vertex_out

  @todo current algorithm iterates 10 times at most, since (I guess) there is no termination proof.
  */
  template<class PolygonMesh,
           class InputIterator,
           class FacetOutputIterator,
           class VertexOutputIterator>
  std::pair<FacetOutputIterator, VertexOutputIterator>
    refine(PolygonMesh& pmesh,
           InputIterator facet_begin,
           InputIterator facet_end,
           FacetOutputIterator facet_out,
           VertexOutputIterator vertex_out,
           double density_control_factor = std::sqrt(2.0))
  {
    internal::Refine_Polyhedron_3<PolygonMesh> refine_functor(pmesh);
    refine_functor.refine(facet_begin,
                          facet_end,
                          facet_out,
                          vertex_out,
                          density_control_factor);
    return std::make_pair(facet_out, vertex_out);
  }

}//end namespace Polygon_mesh_processing

}//end namespace CGAL
