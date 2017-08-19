#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H

#include <CGAL/internal/Surface_mesh_approximation/VSA.h>
#include <CGAL/VSA_approximation.h>
#include <CGAL/internal/Surface_mesh_approximation/named_function_params.h>
#include <CGAL/vsa_mesh_approximation_traits.h>
#include <CGAL/property_map.h>
#include <boost/graph/graph_traits.hpp>

namespace CGAL
{
/*!
 * \ingroup PkgTSMA
 * @brief variational shape approximation a triangulated mesh.
 * This function approximate the input triangulated mesh by fitting it with proxies.
 *
 * @tparam TriangleMesh model of `FaceGraph`.
 *         The descriptor types `boost::graph_traits<TriangleMesh>::%face_descriptor`
 *         and `boost::graph_traits<TriangleMesh>::%halfedge_descriptor` must be
 *         models of `Hashable`.
 *         If `TriangleMesh` has an internal property map for `CGAL::face_index_t`,
 *         and no `face_index_map` is given
 *         as a named parameter, then the internal one should be initialized
 * @tparam FaceRange range of `boost::graph_traits<TriangleMesh>::%face_descriptor`,
          model of `Range`. Its iterator type is `ForwardIterator`.
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param pmesh a polygon mesh with triangulated surface patches to be remeshed
 * @param faces the range of triangular faces defining one or several surface patches to be remeshed
 * @param target_edge_length the edge length that is targetted in the remeshed patch
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
 *    Exact constructions kernels are not supported by this function.
 *  \cgalParamEnd
 *  \cgalParamBegin{vertex_point_map} the property map with the points associated
 *    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
 *  \cgalParamEnd
 *  \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh`
 *  \cgalParamEnd
 *  \cgalParamBegin{number_of_iterations} the number of iterations for the
 *    sequence of atomic operations performed (listed in the above description)
 *  \cgalParamEnd
 *  \cgalParamBegin{face_patch_map} a property map with the patch id's associated to the
     faces of `faces`. Instance of a class model of `ReadWritePropertyMap`. It gets
     updated during the remeshing process while new faces are created.
 *  \cgalParamEnd
 * \cgalNamedParamsEnd
 */
template <typename TriangleMesh,
  typename NamedParameters>
bool vsa_mesh_approximation(const TriangleMesh &tm_in,
  TriangleMesh &tm_out,
  const NamedParameters &np)
{
  using boost::get_param;
  using boost::choose_param;

  // typedef typename GetGeomTraits<PM, NamedParameters>::type GeomTraits;
  typedef typename TriangleMesh::Traits GeomTraits;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Vector_3 Vector_3;

  // typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type VPMap;
  // VPMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
  //   get_property_map(vertex_point, tm_in));

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef boost::associative_property_map<std::map<face_descriptor, Vector_3> > FacetNormalMap;
  typedef boost::associative_property_map<std::map<face_descriptor, FT> > FacetAreaMap;
  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type VertexPointMap;

  typedef CGAL::PlaneProxy<TriangleMesh> PlaneProxy;
  typedef CGAL::L21Metric<TriangleMesh, FacetNormalMap, FacetAreaMap> L21Metric;
  typedef CGAL::L21ProxyFitting<TriangleMesh, FacetNormalMap, FacetAreaMap> L21ProxyFitting;
  typedef CGAL::VSA_approximation<TriangleMesh, PlaneProxy, L21Metric, L21ProxyFitting> VSAL21;

  VertexPointMap point_pmap = get(boost::vertex_point, const_cast<TriangleMesh &>(tm_in));
  // construct facet normal & area map
  std::map<face_descriptor, Vector_3> facet_normals;
  std::map<face_descriptor, FT> facet_areas;
  BOOST_FOREACH(face_descriptor f, faces(tm_in)) {
    const halfedge_descriptor he = halfedge(f, tm_in);
    const Point_3 p0 = point_pmap[source(he, tm_in)];
    const Point_3 p1 = point_pmap[target(he, tm_in)];
    const Point_3 p2 = point_pmap[target(next(he, tm_in), tm_in)];
    Vector_3 normal = CGAL::unit_normal(p0, p1, p2);
    facet_normals.insert(std::pair<face_descriptor, Vector_3>(f, normal));
    FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    facet_areas.insert(std::pair<face_descriptor, FT>(f, area));
  }
  FacetNormalMap normal_pmap(facet_normals);
  FacetAreaMap area_pmap(facet_areas);

  L21Metric l21_metric(normal_pmap, area_pmap);
  L21ProxyFitting l21_fitting(normal_pmap, area_pmap);

  VSAL21 vsa_l21(l21_metric, l21_fitting);
  vsa_l21.set_mesh(tm_in);

  std::size_t num_proxies = choose_param(get_param(np, internal_np::number_of_proxies),
    num_faces(tm_in) / 100);
  std::size_t num_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 10);
  std::cout << "#px = " << num_proxies << ", #itr = " << num_iterations << std::endl;
  // bool pca_plane = choose_param(get_param(np, internal_np::pca_plane), false);

  vsa_l21.init_proxies(num_proxies, VSAL21::RandomInit);
  for (std::size_t i = 0; i < num_iterations; ++i)
    vsa_l21.run_one_step();

  // vsa_l21.get_proxy_map();

  return vsa_l21.meshing(tm_out);
  // vsa_l21.get_anchor_points();
}

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
 * @tparam ErrorMetric error metric type
 * @tparam ProxyFitting proxy fitting type

 * @param tm a triangle mesh
 * @param[out] f_proxy_pmap facet proxy patch id property map
 * @param fit_error error metric functor
 * @param proxy_fitting proxy fitting functor
 * @param init select seed initialization
 * @param number_of_segments target number of approximation patches
 * @param number_of_iterations number of fitting iterations
 */
template<typename TriangleMesh,
  typename FacetProxyMap,
  typename ErrorMetric,
  typename ProxyFitting>
void vsa_approximate(
    const TriangleMesh &tm,
    FacetProxyMap &f_proxy_pmap,
    const ErrorMetric &fit_error,
    const ProxyFitting &proxy_fitting,
    const int init,
    const std::size_t number_of_segments,
    const std::size_t number_of_iterations) {
  // CGAL_precondition(is_pure_triangle(tm));

  typedef CGAL::internal::VSA_approximation<
    TriangleMesh,
    FacetProxyMap,
    ErrorMetric,
    ProxyFitting> VSA_approximation;

  VSA_approximation algorithm(tm, fit_error, proxy_fitting);

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
 * @tparam PlaneFitting a plane fitting functor
 * @tparam ErrorMetric error metric type
 * @tparam ProxyFitting proxy fitting type

 * @param tm a triangle mesh
 * @param[out] tris approximation indexed triangle soup
 * @param[out] pos anchor position container
 * @param plane_fitting plane fitting functor
 * @param fit_error error metric functor
 * @param proxy_fitting proxy fitting functor
 * @param init select seed initialization
 * @param number_of_segments target number of approximation patches
 * @param number_of_iterations number of fitting iterations
 */
template<typename TriangleMesh,
  typename AnchorIndexContainer,
  typename AnchorPositionContainer,
  typename PlaneFitting,
  typename ErrorMetric,
  typename ProxyFitting>
void vsa_extract(
    const TriangleMesh &tm,
    AnchorIndexContainer &tris,
    AnchorPositionContainer &pos,
    const PlaneFitting &plane_fitting,
    const ErrorMetric &fit_error,
    const ProxyFitting &proxy_fitting,
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

  vsa_approximate(tm,
    f_proxy_pmap,
    fit_error,
    proxy_fitting,
    init,
    number_of_segments,
    number_of_iterations);

  typedef CGAL::internal::VSA_mesh_extraction<
    TriangleMesh,
    FacetProxyMap,
    PlaneFitting> VSA_mesh_extraction;

  VSA_mesh_extraction extractor(tm, f_proxy_pmap, plane_fitting);

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
 * @tparam PlaneFitting a plane fitting functor
 * @tparam ErrorMetric error metric type
 * @tparam ProxyFitting proxy fitting type

 * @param tm a triangle mesh
 * @param[out] f_proxy_pmap facet proxy patch id property map
 * @param[out] tris approximation indexed triangle soup
 * @param[out] pos anchor position container
 * @param plane_fitting plane fitting functor
 * @param fit_error error metric functor
 * @param proxy_fitting proxy fitting functor
 * @param init select seed initialization
 * @param number_of_segments target number of approximation patches
 * @param number_of_iterations number of fitting iterations
 */
template<typename TriangleMesh,
  typename FacetProxyMap,
  typename AnchorIndexContainer,
  typename AnchorPositionContainer,
  typename PlaneFitting,
  typename ErrorMetric,
  typename ProxyFitting>
void vsa_approximate_and_extract(
    const TriangleMesh &tm,
    FacetProxyMap f_proxy_pmap,
    AnchorIndexContainer &tris,
    AnchorPositionContainer &pos,
    const PlaneFitting &plane_fitting,
    const ErrorMetric &fit_error,
    const ProxyFitting &proxy_fitting,
    const int init,
    const std::size_t number_of_segments,
    const std::size_t number_of_iterations) {
  // CGAL_precondition(is_pure_triangle(tm));

  vsa_approximate(tm,
    f_proxy_pmap,
    fit_error,
    proxy_fitting,
    init,
    number_of_segments,
    number_of_iterations);

  typedef CGAL::internal::VSA_mesh_extraction<
    TriangleMesh,
    FacetProxyMap,
    PlaneFitting> VSA_mesh_extraction;

  VSA_mesh_extraction extractor(tm, f_proxy_pmap, plane_fitting);

  extractor.extract_mesh(tris);
  BOOST_FOREACH(const typename VSA_mesh_extraction::Anchor &a, extractor.collect_anchors())
    pos.push_back(a.pos);
}

/*!
 * \ingroup PkgTSMA
 * @brief Using the VSA algorithm to approximate a triangle mesh.
 * This function approximate the input mesh by fitting it with proxies.
 * Facet proxy index output only.
 * Use default L21 metric.
 *
 * @tparam TriangleMesh model of `FaceGraph`.
 * @tparam FacetProxyMap a property map containing the approximated facet patch id,
           and `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type,
           std::size_t as value type

 * @param tm a triangle mesh
 * @param[out] f_proxy_pmap facet proxy patch id property map
 * @param init select seed initialization
 * @param number_of_segments target number of approximation patches
 * @param number_of_iterations number of fitting iterations
 */
template<typename TriangleMesh,
  typename FacetProxyMap>
void vsa_approximate(
    const TriangleMesh &tm,
    FacetProxyMap f_proxy_pmap,
    const int init,
    const std::size_t number_of_segments,
    const std::size_t number_of_iterations) {
  // CGAL_precondition(is_pure_triangle(tm));

  typedef typename TriangleMesh::Traits GeomTraits;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef boost::associative_property_map<std::map<face_descriptor, Vector_3> > FacetNormalMap;
  typedef boost::associative_property_map<std::map<face_descriptor, FT> > FacetAreaMap;
  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type VertexPointMap;

  typedef CGAL::L21Metric<TriangleMesh, FacetNormalMap, FacetAreaMap> L21Metric;
  typedef CGAL::L21ProxyFitting<TriangleMesh, FacetNormalMap, FacetAreaMap> L21ProxyFitting;

  VertexPointMap point_pmap = get(boost::vertex_point, const_cast<TriangleMesh &>(tm));
  // construct facet normal & area map
  std::map<face_descriptor, Vector_3> facet_normals;
  std::map<face_descriptor, FT> facet_areas;
  BOOST_FOREACH(face_descriptor f, faces(tm)) {
    const halfedge_descriptor he = halfedge(f, tm);
    const Point_3 p0 = point_pmap[source(he, tm)];
    const Point_3 p1 = point_pmap[target(he, tm)];
    const Point_3 p2 = point_pmap[target(next(he, tm), tm)];
    Vector_3 normal = CGAL::unit_normal(p0, p1, p2);
    facet_normals.insert(std::pair<face_descriptor, Vector_3>(f, normal));
    FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    facet_areas.insert(std::pair<face_descriptor, FT>(f, area));
  }
  FacetNormalMap normal_pmap(facet_normals);
  FacetAreaMap area_pmap(facet_areas);

  vsa_approximate(tm,
    f_proxy_pmap,
    L21Metric(normal_pmap, area_pmap),
    L21ProxyFitting(normal_pmap, area_pmap),
    init,
    number_of_segments,
    number_of_iterations);
}

/*!
 * \ingroup PkgTSMA
 * @brief Using the VSA algorithm to extract the approximated surface mesh.
 * This function approximate the input triangulated mesh by fitting it with proxies.
 * Approximated triangle mesh output only.
 * Use default L21 metric.
 *
 * @tparam TriangleMesh model of `FaceGraph`.
 * @tparam AnchorIndexContainer a container of approximated indexed triangle soup
 * @tparam AnchorPositionContainer a container of extracted anchor position

 * @param tm a triangle mesh
 * @param[out] tris approximation indexed triangle soup
 * @param[out] pos anchor position container
 * @param init select seed initialization
 * @param number_of_segments target number of approximation patches
 * @param number_of_iterations number of fitting iterations
 */
template<typename TriangleMesh,
  typename AnchorIndexContainer,
  typename AnchorPositionContainer>
void vsa_extract(
    const TriangleMesh &tm,
    AnchorIndexContainer &tris,
    AnchorPositionContainer &pos,
    const int init,
    const std::size_t number_of_segments,
    const std::size_t number_of_iterations) {
  // CGAL_precondition(is_pure_triangle(tm));

  typedef typename TriangleMesh::Traits GeomTraits;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef boost::associative_property_map<std::map<face_descriptor, Vector_3> > FacetNormalMap;
  typedef boost::associative_property_map<std::map<face_descriptor, FT> > FacetAreaMap;
  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type VertexPointMap;

  typedef CGAL::L21Metric<TriangleMesh, FacetNormalMap, FacetAreaMap> L21Metric;
  typedef CGAL::L21ProxyFitting<TriangleMesh, FacetNormalMap, FacetAreaMap> L21ProxyFitting;
  typedef CGAL::PlaneFitting<TriangleMesh> PlaneFitting;

  VertexPointMap point_pmap = get(boost::vertex_point, const_cast<TriangleMesh &>(tm));
  // construct facet normal & area map
  std::map<face_descriptor, Vector_3> facet_normals;
  std::map<face_descriptor, FT> facet_areas;
  BOOST_FOREACH(face_descriptor f, faces(tm)) {
    const halfedge_descriptor he = halfedge(f, tm);
    const Point_3 p0 = point_pmap[source(he, tm)];
    const Point_3 p1 = point_pmap[target(he, tm)];
    const Point_3 p2 = point_pmap[target(next(he, tm), tm)];
    Vector_3 normal = CGAL::unit_normal(p0, p1, p2);
    facet_normals.insert(std::pair<face_descriptor, Vector_3>(f, normal));
    FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    facet_areas.insert(std::pair<face_descriptor, FT>(f, area));
  }
  FacetNormalMap normal_pmap(facet_normals);
  FacetAreaMap area_pmap(facet_areas);

  typedef typename boost::associative_property_map<std::map<face_descriptor, std::size_t> > FacetProxyMap;
  std::map<face_descriptor, std::size_t> facet_proxy_map;
  BOOST_FOREACH(face_descriptor f, faces(tm))
    facet_proxy_map.insert(std::pair<face_descriptor, std::size_t>(f, 0));
  FacetProxyMap f_proxy_pmap(facet_proxy_map);

  vsa_extract(tm,
    tris,
    pos,
    PlaneFitting(tm),
    L21Metric(normal_pmap, area_pmap),
    L21ProxyFitting(normal_pmap, area_pmap),
    init,
    number_of_segments,
    number_of_iterations);
}

/*!
 * \ingroup PkgTSMA
 * @brief Triangle mesh variational shape approximation.
 * This function approximate the input triangulated mesh by fitting it with proxies.
 * Facet proxy index and approximated triangle mesh output.
 * Use default L21 metric.
 *
 * @tparam TriangleMesh model of `FaceGraph`.
 * @tparam FacetProxyMap a property map containing the approximated facet patch id,
           and `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type,
           std::size_t as value type
 * @tparam AnchorIndexContainer a container of approximated indexed triangle soup
 * @tparam AnchorPositionContainer a container of extracted anchor position

 * @param tm a triangle mesh
 * @param[out] f_proxy_pmap facet proxy patch id property map
 * @param[out] tris approximation indexed triangle soup
 * @param[out] pos anchor position container
 * @param init select seed initialization
 * @param number_of_segments target number of approximation patches
 * @param number_of_iterations number of fitting iterations
 */
template<typename TriangleMesh,
  typename FacetProxyMap,
  typename AnchorIndexContainer,
  typename AnchorPositionContainer>
void vsa_approximate_and_extract(
    const TriangleMesh &tm,
    FacetProxyMap f_proxy_pmap,
    AnchorIndexContainer &tris,
    AnchorPositionContainer &pos,
    const int init,
    const std::size_t number_of_segments,
    const std::size_t number_of_iterations) {
  // CGAL_precondition(is_pure_triangle(tm));

  typedef typename TriangleMesh::Traits GeomTraits;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef boost::associative_property_map<std::map<face_descriptor, Vector_3> > FacetNormalMap;
  typedef boost::associative_property_map<std::map<face_descriptor, FT> > FacetAreaMap;
  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type VertexPointMap;

  typedef CGAL::L21Metric<TriangleMesh, FacetNormalMap, FacetAreaMap> L21Metric;
  typedef CGAL::L21ProxyFitting<TriangleMesh, FacetNormalMap, FacetAreaMap> L21ProxyFitting;
  typedef CGAL::PlaneFitting<TriangleMesh> PlaneFitting;

  VertexPointMap point_pmap = get(boost::vertex_point, const_cast<TriangleMesh &>(tm));
  // construct facet normal & area map
  std::map<face_descriptor, Vector_3> facet_normals;
  std::map<face_descriptor, FT> facet_areas;
  BOOST_FOREACH(face_descriptor f, faces(tm)) {
    const halfedge_descriptor he = halfedge(f, tm);
    const Point_3 p0 = point_pmap[source(he, tm)];
    const Point_3 p1 = point_pmap[target(he, tm)];
    const Point_3 p2 = point_pmap[target(next(he, tm), tm)];
    Vector_3 normal = CGAL::unit_normal(p0, p1, p2);
    facet_normals.insert(std::pair<face_descriptor, Vector_3>(f, normal));
    FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    facet_areas.insert(std::pair<face_descriptor, FT>(f, area));
  }
  FacetNormalMap normal_pmap(facet_normals);
  FacetAreaMap area_pmap(facet_areas);

  vsa_approximate_and_extract(tm,
    f_proxy_pmap,
    tris,
    pos,
    PlaneFitting(tm),
    L21Metric(normal_pmap, area_pmap),
    L21ProxyFitting(normal_pmap, area_pmap),
    init,
    number_of_segments,
    number_of_iterations);
}
} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
