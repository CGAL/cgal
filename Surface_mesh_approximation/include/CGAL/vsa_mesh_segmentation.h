#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_SEGMENTATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_SEGMENTATION_H

#include <CGAL/license/Surface_mesh_approximation.h>


#include <CGAL/vsa_metrics.h>
#include <CGAL/vsa_approximation.h>
#include <CGAL/internal/Surface_mesh_approximation/named_function_params.h>
#include <CGAL/internal/Surface_mesh_approximation/named_params_helper.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/optional.hpp>

namespace CGAL {

/*!
 * \ingroup PkgTSMA
 * @brief Function segments the input mesh into planar parts.
 * This function segments the input triangle surface mesh by the Variational Shape Approximation algorithm.
 * It fills a property map which associates a segment-id (in [0, number_of_segments -1])
 * or a proxy-id (in [0, number_of_proxies - 1]) to each facet.
 * A segment is a set of connected facets which are placed under the same proxy patch (see \cgalFigureRef{iterations}).
 *
 * @tparam TriangleMesh model of `FaceGraph`.
 *         The descriptor types `boost::graph_traits<TriangleMesh>::%face_descriptor`
 *         and `boost::graph_traits<TriangleMesh>::%halfedge_descriptor` must be
 *         models of `Hashable`.
 *         If `TriangleMesh` has an internal property map for `CGAL::face_index_t`,
 *         and no `face_index_map` is given
 *         as a named parameter, then the internal one should be initialized
 * @tparam SegmentPropertyMap a ReadWritePropertyMap with
 * `boost::graph_traits<TriangleMesh>::%face_descriptor` as key and `std::size_t` as value type
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param tm_in a triangulated surface mesh to be segmented
 * @param[out] segment_ids the segment or proxy id of each facet
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
 *    Exact constructions kernels are not supported by this function.
 *  \cgalParamEnd
 *  \cgalParamBegin{vertex_point_map} property map with the points associated
 *    to the vertices of `tm_in`. Instance of a class model of `ReadWritePropertyMap`.
 *  \cgalParamEnd
 *  \cgalParamBegin{seeding_method} selection of seeding method.
 *  \cgalParamEnd
 *  \cgalParamBegin{max_nb_proxies} maximum number of proxies to approximate the geometry.
 *  \cgalParamEnd
 *  \cgalParamBegin{min_error_drop} minimum error drop of the approximation.
 *  \cgalParamEnd
 *  \cgalParamBegin{nb_of_iterations} number of partitioning and fitting iterations after seeding.
 *  \cgalParamEnd
 *  \cgalParamBegin{nb_of_relaxations} number of relaxations interleaved within seeding.
 *  \cgalParamEnd
 *  \cgalParamBegin{proxies} output iterator over proxies
 *  \cgalParamEnd
 * \cgalNamedParamsEnd
 */
template <typename TriangleMesh,
  typename SegmentPropertyMap,
  typename NamedParameters>
void vsa_mesh_segmentation(const TriangleMesh &tm_in,
  const SegmentPropertyMap segment_ids, const NamedParameters &np)
{
  using boost::get_param;
  using boost::choose_param;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;
  typedef typename Geom_traits::FT FT;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type Vertex_point_map;
  Vertex_point_map point_pmap = choose_param(get_param(np, internal_np::vertex_point),
    get_property_map(vertex_point, const_cast<TriangleMesh &>(tm_in)));

  typedef CGAL::VSA::Mesh_approximation<TriangleMesh, Vertex_point_map> L21_approx;
  typedef typename L21_approx::Error_metric L21_metric;
  typedef typename L21_approx::Proxy_fitting L21_proxy_fitting;

  L21_approx approx(tm_in, point_pmap);
  L21_metric l21_metric(tm_in);
  L21_proxy_fitting l21_fitting(tm_in);
  approx.set_metric(l21_metric, l21_fitting);

  // default hierarchical seeding
  CGAL::VSA::Seeding method = choose_param(
    get_param(np, internal_np::seeding_method), CGAL::VSA::Hierarchical);
  boost::optional<std::size_t> max_nb_proxies = choose_param(
    get_param(np, internal_np::max_nb_proxies), boost::optional<std::size_t>());
  boost::optional<FT> min_error_drop = choose_param(
    get_param(np, internal_np::min_error_drop), boost::optional<FT>());
  std::size_t nb_of_relaxations = choose_param(get_param(np, internal_np::nb_of_relaxations), 5);
  approx.seeding(method, max_nb_proxies, min_error_drop, nb_of_relaxations);

  // reasonable default number of iterations
  std::size_t nb_of_iterations_default = max_nb_proxies ? num_faces(tm_in) / *max_nb_proxies : 30;
  nb_of_iterations_default = (std::min)((std::max)(
    nb_of_iterations_default, static_cast<std::size_t>(20)), static_cast<std::size_t>(60));
  const std::size_t nb_of_iterations = choose_param(
    get_param(np, internal_np::nb_of_iterations), nb_of_iterations_default);

  approx.run(nb_of_iterations);

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  std::cout << "#px = " << approx.get_proxies_sizes()
    << ", #itr = " << nb_of_iterations
    << ", #relx = " << nb_of_relaxations << std::endl;
#endif

  approx.proxy_map(segment_ids);

  typedef typename boost::lookup_named_param_def <
    internal_np::proxies_t,
    NamedParameters,
    internal_np::vsa_no_output_t>::type ProxiesOutItr;
  ProxiesOutItr pxies_out_itr = choose_param(
    get_param(np, internal_np::proxies), internal_np::vsa_no_output);
  get_proxies(approx, pxies_out_itr);
}

} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_SEGMENTATION_H
