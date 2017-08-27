#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H

#include <CGAL/VSA_metrics.h>
#include <CGAL/VSA_approximation.h>
#include <CGAL/internal/Surface_mesh_approximation/named_function_params.h>
#include <CGAL/internal/Surface_mesh_approximation/named_params_helper.h>

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
 * @tparam PolyhedronSurface should be `CGAL::Polyhedron_3`
 * @tparam FaceRange range of `boost::graph_traits<TriangleMesh>::%face_descriptor`,
          model of `Range`. Its iterator type is `ForwardIterator`.
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param tm_in a polygon mesh with triangulated surface patches to be remeshed
 * @param[out] tm_out approximated polyhedron mesh
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
 *    Exact constructions kernels are not supported by this function.
 *  \cgalParamEnd
 *  \cgalParamBegin{vertex_point_map} the property map with the points associated
 *    to the vertices of `tm_in`. Instance of a class model of `ReadWritePropertyMap`.
 *  \cgalParamEnd
 *  \cgalParamBegin{init_method} the selection of seed initialization method.
 *  \cgalParamEnd
 *  \cgalParamBegin{number_of_proxies} the number of proxies to approximate the geometry.
 *  \cgalParamEnd
 *  \cgalParamBegin{number_of_iterations} the relaxation iterations.
 *  \cgalParamEnd
 *  \cgalParamBegin{chord_subdivide} the threshold of chord subdivision.
 *  \cgalParamEnd
 *  \cgalParamBegin{face_proxy_map} a property map containing the assigned proxy index of each face of `tm_in`
 *  \cgalParamEnd
 *  \cgalParamBegin{anchor_vertex} the anchor verteices output iterator
 *  \cgalParamEnd
 *  \cgalParamBegin{anchor_point} the anchor points output iterator
 *  \cgalParamEnd
 *  \cgalParamBegin{indexed_triangles} the indexed triangles output iterator
 *  \cgalParamEnd
 * \cgalNamedParamsEnd
 */
template <typename TriangleMesh,
  typename PolyhedronSurface,
  typename NamedParameters>
bool vsa_mesh_approximation(const TriangleMesh &tm_in,
  PolyhedronSurface &tm_out,
  const NamedParameters &np)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  using boost::get_param;
  using boost::choose_param;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GeomTraits;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type VPMap;
  VPMap point_pmap = choose_param(get_param(np, internal_np::vertex_point),
    get(boost::vertex_point, const_cast<TriangleMesh &>(tm_in)));
    // get_property_map(vertex_point, tm_in));

  typedef CGAL::VSA_approximation<TriangleMesh, VPMap> VSAL21;
  typedef typename VSAL21::ErrorMetric L21Metric;
  typedef typename VSAL21::ProxyFitting L21ProxyFitting;

  VSAL21 vsa_l21(tm_in, point_pmap);
  L21Metric l21_metric(tm_in);
  L21ProxyFitting l21_fitting(tm_in);
  vsa_l21.set_metric(l21_metric, l21_fitting);

  std::size_t num_proxies = choose_param(get_param(np, internal_np::number_of_proxies),
    num_faces(tm_in) / 100);
  std::size_t num_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 10);
  std::cout << "#px = " << num_proxies << ", #itr = " << num_iterations << std::endl;

  int init = choose_param(get_param(np, internal_np::init_method), 0);
  if (init < 0 || init > 2)
    return false;

  vsa_l21.init_proxies(num_proxies, static_cast<typename VSAL21::Initialization>(init));
  for (std::size_t i = 0; i < num_iterations; ++i)
    vsa_l21.run_one_step();

  bool pca_plane = choose_param(get_param(np, internal_np::pca_plane) , false);
  FT split_criterion = choose_param(get_param(np, internal_np::chord_subdivide), FT(1));
  bool is_manifold = vsa_l21.meshing(tm_out, split_criterion, pca_plane);

  // vsa_l21.get_proxy_map();

  typedef typename boost::lookup_named_param_def <
      internal_np::anchor_vertex_t,
      NamedParameters,
      std::back_insert_iterator<std::vector<vertex_descriptor> >
    > ::type AnchorVertexOutItr;
  AnchorVertexOutItr avtx_out_itr = choose_param(get_param(np, internal_np::anchor_vertex)
    , std::back_inserter(*(new std::vector<vertex_descriptor>())));
  vsa_l21.get_anchor_vertices(avtx_out_itr);

  typedef typename boost::lookup_named_param_def <
      internal_np::anchor_point_t,
      NamedParameters,
      std::back_insert_iterator<std::vector<Point_3> >
    > ::type AnchorPointOutItr;
  AnchorPointOutItr apts_out_itr = choose_param(get_param(np, internal_np::anchor_point)
    , std::back_inserter(*(new std::vector<Point_3>())));
  vsa_l21.get_anchor_points(apts_out_itr);

  typedef typename boost::lookup_named_param_def <
      internal_np::indexed_triangles_t,
      NamedParameters,
      std::back_insert_iterator<std::vector<int> >
    > ::type IndexedTrisOutItr;
  IndexedTrisOutItr itris_out_itr = choose_param(get_param(np, internal_np::indexed_triangles)
    , std::back_inserter(*(new std::vector<int>())));
  vsa_l21.get_indexed_triangles(itris_out_itr);

  
  // typedef typename boost::lookup_named_param_def <
  //     internal_np::proxies_t,
  //     NamedParameters,
  //     std::back_insert_iterator<std::vector<PlaneProxy> >
  //   > ::type ProxiesOutItr;
  // ProxiesOutItr pxies_out_itr = choose_param(get_param(np, internal_np::proxies)
  //   , std::back_inserter(*(new std::vector<PlaneProxy>())));
  // vsa_l21.get_proxies(pxies_out_itr);

  return is_manifold;
}

} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
