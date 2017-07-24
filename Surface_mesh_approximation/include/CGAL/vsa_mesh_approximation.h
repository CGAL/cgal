#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H

#include <CGAL/internal/Surface_mesh_approximation/VSA.h>
#include <CGAL/property_map.h>

namespace CGAL
{
/*!
 * \ingroup PkgTSMA
 * @brief variational shape approximation a triangulated mesh.
 * This function approximate the input triangulated mesh by fitting it with proxies.
 *
 * @tparam TriangleMesh model of `FaceGraph`.
 * @tparam FacetProxyMap a property map containing the approximated facet patch id,
           and `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type,
           std::size_t as value type
 * @tparam VertexPointMap a property map containing the input mesh vertex point map,
           and `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type,
           `TriangleMesh::Point_3` as value type
 * @tparam FacetAreaMap a property map containing the input mesh area map,
           and `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type,
           `GeomTraits::FT` as value type
 * @tparam AnchorIndexContainer a container of approximated indexed triangle soup
 * @tparam AnchorPositionContainer a container of extracted anchor position
 * @tparam AnchorVertexContainer a container of extracted anchor vertex
 * @tparam BoundaryContainer a container of proxy patch boundary
 * @tparam ApproximationTrait an approximation trait

 * @param init select seed initialization
 * @param tm a triangle mesh
 * @param number_of_segments target number of approximation patches
 * @param number_of_iterations number of fitting iterations
 * @param f_proxy_pmap facet proxy patch id property map
 * @param v_point_pmap mesh vertex point property map
 * @param f_area_pmap facet area property map
 * @param tris approximation indexed triangle soup
 * @param pos anchor position container
 * @param vtx anchor vertex container
 * @param bdrs proxy patch boundary container
 * @param app_trait shape approximation trait
 */
template<typename TriangleMesh,
  typename FacetProxyMap,
  typename VertexPointMap,
  typename FacetAreaMap,
  typename AnchorIndexContainer,
  typename AnchorPositionContainer,
  typename AnchorVertexContainer,
  typename BoundaryContainer,
  typename ApproximationTrait>
  void vsa_mesh_approximation(
    const int init,
    const TriangleMesh &tm,
    const std::size_t number_of_segments,
    const std::size_t number_of_iterations,
    FacetProxyMap f_proxy_pmap,
    const VertexPointMap &v_point_pmap,
    const FacetAreaMap &f_area_pmap,
    AnchorIndexContainer &tris,
    AnchorPositionContainer &pos,
    AnchorVertexContainer &vtx,
    BoundaryContainer &bdrs,
    const ApproximationTrait &app_trait) {
  // CGAL_precondition(is_pure_triangle(tm));

  typedef CGAL::internal::VSA_approximation<
    TriangleMesh,
    FacetProxyMap,
    ApproximationTrait> VSA_approximation;

  VSA_approximation algorithm(tm, app_trait);

  switch (init) {
    case VSA_approximation::RandomInit:
      algorithm.partition(number_of_segments, number_of_iterations, f_proxy_pmap);
      break;
    case VSA_approximation::IncrementalInit:
      algorithm.partition_incre(number_of_segments, number_of_iterations, f_proxy_pmap);
      break;
    case VSA_approximation::HierarchicalInit:
      algorithm.partition_hierarchical(number_of_segments, number_of_iterations, f_proxy_pmap);
      break;
  }

  typedef CGAL::internal::VSA_mesh_extraction<
    TriangleMesh,
    ApproximationTrait,
    VertexPointMap,
    FacetProxyMap> VSA_mesh_extraction;

  VSA_mesh_extraction extractor(tm, app_trait, v_point_pmap, f_proxy_pmap);

  extractor.extract_mesh(tris);
  BOOST_FOREACH(const typename VSA_mesh_extraction::Anchor &a, extractor.collect_anchors()) {
    vtx.push_back(a.vtx);
    pos.push_back(a.pos);
  }

  bdrs = extractor.collect_borders();
}
}

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
