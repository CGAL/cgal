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
 * @tparam SegmentPropertyMap a property map containing the approximated facet patch id,
           and `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type,
           std::size_t as value type
 * @tparam PointPropertyMap a property map containing the input mesh vertex point map,
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
 * @tparam GeomTraits geometric kernel

 * @param init select seed initialization
 * @param tm a triangle mesh
 * @param number_of_segments target number of approximation patches
 * @param number_of_iterations number of fitting iterations
 * @param segment_ids facet proxy patch id property map
 * @param ppmap mesh vertex point property map
 * @param area_pmap facet area property map
 * @param tris approximation indexed triangle soup
 * @param pos anchor position container
 * @param vtx anchor vertex container
 * @param bdrs proxy patch boundary container
 * @param app_trait shape approximation trait
 * @param traits kernel traits
 */
template<typename TriangleMesh,
  typename SegmentPropertyMap,
  typename PointPropertyMap,
  typename FacetAreaMap,
  typename AnchorIndexContainer,
  typename AnchorPositionContainer,
  typename AnchorVertexContainer,
  typename BoundaryContainer,
  typename ApproximationTrait,
  typename GeomTraits>
  void vsa_mesh_approximation(
    const int init,
    const TriangleMesh &tm,
    const std::size_t number_of_segments,
    const std::size_t number_of_iterations,
    SegmentPropertyMap segment_ids,
    const PointPropertyMap &ppmap,
    const FacetAreaMap &area_pmap,
    AnchorIndexContainer &tris,
    AnchorPositionContainer &pos,
    AnchorVertexContainer &vtx,
    BoundaryContainer &bdrs,
    const ApproximationTrait &app_trait,
    GeomTraits traits) {
  // CGAL_precondition(is_pure_triangle(tm));

  typedef CGAL::internal::VSA<
    TriangleMesh,
    SegmentPropertyMap,
    ApproximationTrait> VSA;

  VSA algorithm(tm, app_trait);

  switch (init) {
    case VSA::RandomInit:
      algorithm.partition(number_of_segments, number_of_iterations, segment_ids);
      break;
    case VSA::IncrementalInit:
      algorithm.partition_incre(number_of_segments, number_of_iterations, segment_ids);
      break;
    case VSA::HierarchicalInit:
      algorithm.partition_hierarchical(number_of_segments, number_of_iterations, segment_ids);
      break;
  }

  typedef CGAL::internal::VSA_mesh_extraction<
    TriangleMesh,
    ApproximationTrait,
    PointPropertyMap,
    SegmentPropertyMap,
    FacetAreaMap> VSA_mesh_extraction;

  VSA_mesh_extraction extractor(tm, app_trait, ppmap, segment_ids, area_pmap);

  extractor.extract_mesh(tris);
  BOOST_FOREACH(const typename VSA_mesh_extraction::Anchor &a, extractor.collect_anchors()) {
    vtx.push_back(a.vtx);
    pos.push_back(a.pos);
  }

  bdrs = extractor.collect_borders();
}
}

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
