#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H

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
namespace VSA {

/*!
 * \ingroup PkgTSMA
 * @brief variational shape approximation of a triangle surface mesh.
 * This function approximates the input mesh by fitting it with proxies.
 *
 * @tparam TriangleMesh model of `FaceGraph`.
 *         The descriptor types `boost::graph_traits<TriangleMesh>::%face_descriptor`
 *         and `boost::graph_traits<TriangleMesh>::%halfedge_descriptor` must be
 *         models of `Hashable`.
 *         If `TriangleMesh` has an internal property map for `CGAL::face_index_t`,
 *         and no `face_index_map` is given as a named parameter, then the internal one should be initialized.
 * @tparam AnchorPointOutItr a class model of `OutputIterator` with `GeomTraits::Point_3` value type
 * @tparam IndexedTrisOutItr a class model of `OutputIterator` with `std::vector<std::size_t>` value type
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param tm_in a triangle surface mesh to be approximated
 * @param[out] apts_out_itr output iterator over anchor points 
 * @param[out] tris_out_itr output iterator over indexed triangles (triplets of integers referring to anchor points)
 * @param np optional sequence of \ref namedparameters among the ones listed below
 * @return true if the indexed triangles represent a valid 2-manifold, oriented surface mesh, and false otherwise. 
  *
 * \cgalNamedParamsBegin
 *  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
 *    Exact constructions kernels are not supported by this function.
 *  \cgalParamEnd
 *  \cgalParamBegin{vertex_point_map} the property map with the points associated
 *    to the vertices of `tm_in`. Instance of a class model of `ReadWritePropertyMap`.
 *  \cgalParamEnd
 *  \cgalParamBegin{init_method} selection of initialization method.
 *  \cgalParamEnd
 *  \cgalParamBegin{max_nb_proxies} maximum number of proxies to approximate the geometry.
 *  \cgalParamEnd
 *  \cgalParamBegin{min_error_drop} minimum error drop of the approximation, expressed in ratio between two iterations of proxy addition.
 *  \cgalParamEnd
 *  \cgalParamBegin{nb_of_iterations} number of partitioning and fitting iterations after initialization.
 *  \cgalParamEnd
 *  \cgalParamBegin{nb_of_relaxations} number of relaxations interleaved within initialization seeding.
 *  \cgalParamEnd
 *  \cgalParamBegin{mesh_chord_error} maximum chord approximation error use for mesh construction.
 *  \cgalParamEnd
 *  \cgalParamBegin{face_proxy_map} property map containing the assigned proxy index of each face of input mesh `tm_in`.
 *  \cgalParamEnd
 *  \cgalParamBegin{proxies} output iterator over proxies.
 *  \cgalParamEnd
 *  \cgalParamBegin{anchor_vertices} output iterator over anchor vertices, defined on the input mesh `tm_in`. 
 *  \cgalParamEnd
 *  \cgalParamBegin{output_mesh} polyhedral surface mesh derived from the indexed facet set. The polyhedron is not empty only 
 *  when the indexed face set represents a 2-manifold, oriented surface triangle mesh.
 *  \cgalParamEnd
 * \cgalNamedParamsEnd
 */
template <typename TriangleMesh,
  typename AnchorPointOutItr,
  typename IndexedTrisOutItr,
  typename NamedParameters>
bool mesh_approximation(const TriangleMesh &tm_in,
  AnchorPointOutItr apts_out_itr,
  IndexedTrisOutItr tris_out_itr,
  const NamedParameters &np)
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

  // default hierarchical initialization
  CGAL::VSA::Seeding method = choose_param(
    get_param(np, internal_np::init_method), CGAL::VSA::Hierarchical);
  boost::optional<std::size_t> max_nb_proxies = choose_param(
    get_param(np, internal_np::max_nb_proxies), boost::optional<std::size_t>());
  boost::optional<FT> min_error_drop = choose_param(
    get_param(np, internal_np::min_error_drop), boost::optional<FT>());
  std::size_t nb_of_relaxations = choose_param(get_param(np, internal_np::nb_of_relaxations), 5);
  approx.init(method, max_nb_proxies, min_error_drop, nb_of_relaxations);

  const std::size_t nb_of_iterations = choose_param(get_param(np, internal_np::nb_of_iterations), 10);
  approx.run(nb_of_iterations);

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  std::cout << "#px = " << approx.get_proxies_size()
    << ", #itr = " << nb_of_iterations
    << ", #relx = " << nb_of_relaxations << std::endl;
#endif

  typedef typename boost::lookup_named_param_def<
    internal_np::facet_proxy_map_t,
    NamedParameters,
    internal_np::vsa_no_output_t>::type FPMap;
  FPMap fproxymap = choose_param(
    get_param(np, internal_np::facet_proxy_map), internal_np::vsa_no_output);
  get_proxy_map(approx, fproxymap);

  typedef CGAL::Polyhedron_3<Geom_traits> PolyhedronSurface;
  PolyhedronSurface tmp_poly;
  PolyhedronSurface * const tm_out = choose_param(get_param(np, internal_np::output_mesh), &tmp_poly);
  const FT chord_error = choose_param(get_param(np, internal_np::mesh_chord_error), FT(1.0));
  const bool is_manifold = approx.extract_mesh(*tm_out, chord_error);

  typedef typename boost::lookup_named_param_def<
    internal_np::anchor_vertices_t,
    NamedParameters,
    internal_np::vsa_no_output_t>::type AnchorVertexOutItr;
  AnchorVertexOutItr avtx_out_itr = choose_param(
    get_param(np, internal_np::anchor_vertices) , internal_np::vsa_no_output);
  get_anchor_vertices(approx, avtx_out_itr);

  // get anchor points
  get_anchor_points(approx, apts_out_itr);

  // get indexed triangles
  get_indexed_triangles(approx, tris_out_itr);

  typedef typename boost::lookup_named_param_def <
    internal_np::proxies_t,
    NamedParameters,
    internal_np::vsa_no_output_t>::type ProxiesOutItr;
  ProxiesOutItr pxies_out_itr = choose_param(
    get_param(np, internal_np::proxies), internal_np::vsa_no_output);
  get_proxies(approx, pxies_out_itr);

  return is_manifold;
}

} // end namespace VSA
} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
