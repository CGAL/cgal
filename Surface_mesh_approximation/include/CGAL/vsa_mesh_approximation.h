#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H

#include <CGAL/license/Surface_mesh_approximation.h>


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
 *  \cgalParamBegin{init_by_number} the number of proxies to approximate the geometry.
 *  \cgalParamEnd
 *  \cgalParamBegin{init_by_error} the error drop of the approximation.
 *  \cgalParamEnd
 *  \cgalParamBegin{iterations} the relaxation iterations after seeding.
 *  \cgalParamEnd
 *  \cgalParamBegin{inner_iterations} the relaxation iterations when seeding.
 *  \cgalParamEnd
 *  \cgalParamBegin{chord_subdivide} the threshold of chord subdivision.
 *  \cgalParamEnd
 *  \cgalParamBegin{face_proxy_map} a property map containing the assigned proxy index of each face of `tm_in`
 *  \cgalParamEnd
 *  \cgalParamBegin{proxies} the plane proxies
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
  using boost::get_param;
  using boost::choose_param;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GeomTraits;
  typedef typename GeomTraits::FT FT;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type VPMap;
  VPMap point_pmap = choose_param(get_param(np, internal_np::vertex_point),
    get_property_map(vertex_point, const_cast<TriangleMesh &>(tm_in)));

  typedef CGAL::VSA_approximation<TriangleMesh, VPMap> VSAL21;
  typedef typename VSAL21::ErrorMetric L21Metric;
  typedef typename VSAL21::ProxyFitting L21ProxyFitting;

  VSAL21 vsa_l21(tm_in, point_pmap);
  L21Metric l21_metric(tm_in);
  L21ProxyFitting l21_fitting(tm_in);
  vsa_l21.set_metric(l21_metric, l21_fitting);

  // default random initialization
  int init = choose_param(get_param(np, internal_np::init_method), 0);
  std::size_t num_proxies = choose_param(get_param(np, internal_np::init_by_number), 0);
  std::size_t inner_iterations = choose_param(get_param(np, internal_np::inner_iterations), 10);
  if (num_proxies == 0 || num_proxies > num_faces(tm_in)) {
    FT drop = choose_param(get_param(np, internal_np::init_by_error), FT(0.01));
    vsa_l21.init_by_error(
      static_cast<typename VSAL21::Method>(init), drop, inner_iterations);
  }
  else {
    vsa_l21.init_by_number(
      static_cast<typename VSAL21::Method>(init), num_proxies, inner_iterations);
  }

  std::size_t iterations = choose_param(get_param(np, internal_np::iterations), 10);
  for (std::size_t i = 0; i < iterations; ++i)
    vsa_l21.run_one_step();

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  std::cout << "#px = " << num_proxies
    << ", #itr = " << iterations
    << ", #inner_itr = " << inner_iterations << std::endl;
#endif

  typedef typename boost::lookup_named_param_def<
    internal_np::facet_proxy_map_t,
    NamedParameters,
    internal_np::vsa_no_output_t>::type FPMap;
  FPMap fproxymap = choose_param(
    get_param(np, internal_np::facet_proxy_map), internal_np::vsa_no_output);
  get_proxy_map(vsa_l21, fproxymap);

  FT split_criterion = choose_param(get_param(np, internal_np::chord_subdivide), FT(1));
  bool is_manifold = vsa_l21.meshing(tm_out, split_criterion);

  typedef typename boost::lookup_named_param_def<
    internal_np::anchor_vertex_t,
    NamedParameters,
    internal_np::vsa_no_output_t>::type AnchorVertexOutItr;
  AnchorVertexOutItr avtx_out_itr = choose_param(
    get_param(np, internal_np::anchor_vertex) , internal_np::vsa_no_output);
  get_anchor_vertices(vsa_l21, avtx_out_itr);

  typedef typename boost::lookup_named_param_def<
    internal_np::anchor_point_t,
    NamedParameters,
    internal_np::vsa_no_output_t >::type AnchorPointOutItr;
  AnchorPointOutItr apts_out_itr = choose_param(
    get_param(np, internal_np::anchor_point), internal_np::vsa_no_output);
  get_anchor_points(vsa_l21, apts_out_itr);

  typedef typename boost::lookup_named_param_def <
    internal_np::indexed_triangles_t,
    NamedParameters,
    internal_np::vsa_no_output_t>::type IndexedTrisOutItr;
  IndexedTrisOutItr itris_out_itr = choose_param(
    get_param(np, internal_np::indexed_triangles), internal_np::vsa_no_output);
  get_indexed_triangles(vsa_l21, itris_out_itr);

  typedef typename boost::lookup_named_param_def <
    internal_np::proxies_t,
    NamedParameters,
    internal_np::vsa_no_output_t>::type ProxiesOutItr;
  ProxiesOutItr pxies_out_itr = choose_param(
    get_param(np, internal_np::proxies), internal_np::vsa_no_output);
  get_proxies(vsa_l21, pxies_out_itr);

  return is_manifold;
}

} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
