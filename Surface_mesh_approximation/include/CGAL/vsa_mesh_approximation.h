#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H

#include <CGAL/internal/Surface_mesh_approximation/VSA.h>
#include <CGAL/property_map.h>

namespace CGAL
{
/// @cond SKIP_IN_MANUAL
/*!
 * \ingroup PkgTSMA
 * @brief variational shape approximation a triangulated mesh.
 * This function approximate the input triangulated mesh by fitting it with proxies.
 * Mainly used for debugging
 *
 * @tparam TriangleMesh model of `FaceGraph`.
 * @tparam FacetProxyMap a property map containing the approximated facet patch id,
           and `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type,
           std::size_t as value type
 * @tparam VertexPointMap a property map containing the input mesh vertex point map,
           and `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type,
           `TriangleMesh::Point_3` as value type
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
 * @param tris approximation indexed triangle soup
 * @param pos anchor position container
 * @param vtx anchor vertex container
 * @param bdrs proxy patch boundary container
 * @param app_trait shape approximation trait
 */
template<typename TriangleMesh,
  typename FacetProxyMap,
  typename VertexPointMap,
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
    AnchorIndexContainer &tris,
    AnchorPositionContainer &pos,
    AnchorVertexContainer &vtx,
    BoundaryContainer &bdrs,
    const ApproximationTrait &app_trait) {
  // CGAL_precondition(is_pure_triangle(tm));

  vsa_mesh_approximation(tm,
    f_proxy_pmap,
    app_trait,
    init,
    number_of_segments,
    number_of_iterations);

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
/// @endcond

/*!
 * \ingroup PkgTSMA
 * @brief Using the VSA algorithm to approximate a triangle mesh.
 * This function approximate the input mesh by fitting it with proxies.
 * Facet proxy index output only.
 *
 * @tparam TriangleMesh model of `FaceGraph`.
 * @tparam FacetProxyMap a property map containing the approximated facet patch id,
           and `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type,
           std::size_t as value type
 * @tparam ApproximationTrait an approximation trait

 * @param tm a triangle mesh
 * @param[out] f_proxy_pmap facet proxy patch id property map
 * @param app_trait shape approximation trait
 * @param init select seed initialization
 * @param number_of_segments target number of approximation patches
 * @param number_of_iterations number of fitting iterations
 */
template<typename TriangleMesh,
  typename FacetProxyMap,
  typename ApproximationTrait>
void vsa_mesh_approximation(
    const TriangleMesh &tm,
    FacetProxyMap f_proxy_pmap,
    const ApproximationTrait &app_trait,
    const int init,
    const std::size_t number_of_segments,
    const std::size_t number_of_iterations) {
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
}

/*!
 * \ingroup PkgTSMA
 * @brief Using the VSA algorithm to extract the approximated surface mesh.
 * This function approximate the input triangulated mesh by fitting it with proxies.
 * Approximated triangle mesh output only.
 *
 * @tparam TriangleMesh model of `FaceGraph`.
 * @tparam AnchorIndexContainer a container of approximated indexed triangle soup
 * @tparam AnchorPositionContainer a container of extracted anchor position
 * @tparam ApproximationTrait an approximation trait

 * @param tm a triangle mesh
 * @param[out] tris approximation indexed triangle soup
 * @param[out] pos anchor position container
 * @param app_trait shape approximation trait
 * @param init select seed initialization
 * @param number_of_segments target number of approximation patches
 * @param number_of_iterations number of fitting iterations
 */
template<typename TriangleMesh,
  typename AnchorIndexContainer,
  typename AnchorPositionContainer,
  typename ApproximationTrait>
void vsa_mesh_approximation(
    const TriangleMesh &tm,
    AnchorIndexContainer &tris,
    AnchorPositionContainer &pos,
    const ApproximationTrait &app_trait,
    const int init,
    const std::size_t number_of_segments,
    const std::size_t number_of_iterations) {
  // CGAL_precondition(is_pure_triangle(tm));

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::associative_property_map<std::map<face_descriptor, std::size_t> > FacetProxyMap;
  std::map<face_descriptor, std::size_t> facet_proxy_map;
  BOOST_FOREACH(face_descriptor f, faces(tm))
    facet_proxy_map.insert(std::pair<face_descriptor, std::size_t>(f, 0));
  FacetProxyMap f_proxy_pmap(facet_proxy_map);

  vsa_mesh_approximation(tm,
    f_proxy_pmap,
    app_trait,
    init,
    number_of_segments,
    number_of_iterations);

  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type VertexPointMap;
  typedef CGAL::internal::VSA_mesh_extraction<
    TriangleMesh,
    ApproximationTrait,
    VertexPointMap,
    FacetProxyMap> VSA_mesh_extraction;

  VSA_mesh_extraction extractor(tm,
    app_trait,
    get(boost::vertex_point, const_cast<TriangleMesh &>(tm)),
    f_proxy_pmap);

  extractor.extract_mesh(tris);
  BOOST_FOREACH(const typename VSA_mesh_extraction::Anchor &a, extractor.collect_anchors())
    pos.push_back(a.pos);
}

/*!
 * \ingroup PkgTSMA
 * @brief Triangle mesh variational shape approximation.
 * This function approximate the input triangulated mesh by fitting it with proxies.
 * Facet proxy index and approximated triangle mesh output.
 *
 * @tparam TriangleMesh model of `FaceGraph`.
 * @tparam FacetProxyMap a property map containing the approximated facet patch id,
           and `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type,
           std::size_t as value type
 * @tparam AnchorIndexContainer a container of approximated indexed triangle soup
 * @tparam AnchorPositionContainer a container of extracted anchor position
 * @tparam ApproximationTrait an approximation trait

 * @param tm a triangle mesh
 * @param[out] f_proxy_pmap facet proxy patch id property map
 * @param[out] tris approximation indexed triangle soup
 * @param[out] pos anchor position container
 * @param app_trait shape approximation trait
 * @param init select seed initialization
 * @param number_of_segments target number of approximation patches
 * @param number_of_iterations number of fitting iterations
 */
template<typename TriangleMesh,
  typename FacetProxyMap,
  typename AnchorIndexContainer,
  typename AnchorPositionContainer,
  typename ApproximationTrait>
void vsa_mesh_approximation(
    const TriangleMesh &tm,
    FacetProxyMap f_proxy_pmap,
    AnchorIndexContainer &tris,
    AnchorPositionContainer &pos,
    const ApproximationTrait &app_trait,
    const int init,
    const std::size_t number_of_segments,
    const std::size_t number_of_iterations) {
  // CGAL_precondition(is_pure_triangle(tm));

  vsa_mesh_approximation(tm,
    f_proxy_pmap,
    app_trait,
    init,
    number_of_segments,
    number_of_iterations);

  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type VertexPointMap;
  typedef CGAL::internal::VSA_mesh_extraction<
    TriangleMesh,
    ApproximationTrait,
    VertexPointMap,
    FacetProxyMap> VSA_mesh_extraction;

  VSA_mesh_extraction extractor(tm,
    app_trait,
    get(boost::vertex_point, const_cast<TriangleMesh &>(tm)),
    f_proxy_pmap);

  extractor.extract_mesh(tris);
  BOOST_FOREACH(const typename VSA_mesh_extraction::Anchor &a, extractor.collect_anchors())
    pos.push_back(a.pos);
}
} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
