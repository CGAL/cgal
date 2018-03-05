#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H

#include <CGAL/license/Surface_mesh_approximation.h>


#include <CGAL/VSA_approximation.h>
#include <CGAL/internal/Surface_mesh_approximation/named_function_params.h>
#include <CGAL/internal/Surface_mesh_approximation/named_params_helper.h>

#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/optional.hpp>

namespace CGAL {

/*!
 * \ingroup PkgTSMA
 * @brief Function approximates the input mesh by fitting it with proxies.
 * This function use the Variational Shape Approximation algorithm to approximate the shape of a triangle surface mesh.
 * With indexed triangles as output.
 *
 * @tparam TriangleMesh model of `FaceListGraph`.
 *         If `TriangleMesh` has an internal property map for `CGAL::face_index_t`,
 *         and no `face_index_map` is given as a named parameter, then the internal one should be initialized.
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param tm a triangle surface mesh to be approximated
 * @param np optional sequence of \ref namedparameters among the ones listed below
 * @return `true` if the indexed triangles represent a valid 2-manifold, oriented surface mesh, and `false` otherwise. 
 *
 * \cgalNamedParamsBegin
 *  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
 *    Exact constructions kernels are not supported by this function.
 *  \cgalParamEnd
 *  \cgalParamBegin{vertex_point_map} the property map with the points associated
 *    to the vertices of `tm`. Instance of a class model of `ReadWritePropertyMap`.
 *  \cgalParamEnd
 *  \cgalParamBegin{seeding_method} selection of seeding method.
 *  \cgalParamEnd
 *  \cgalParamBegin{max_nb_proxies} maximum number of proxies to approximate the geometry.
 *  \cgalParamEnd
 *  \cgalParamBegin{min_error_drop} minimum error drop of the approximation, expressed in ratio between two iterations of proxy addition.
 *  \cgalParamEnd
 *  \cgalParamBegin{nb_of_iterations} number of partitioning and fitting iterations after seeding.
 *  \cgalParamEnd
 *  \cgalParamBegin{nb_of_relaxations} number of relaxations interleaved within seeding.
 *  \cgalParamEnd
 *  \cgalParamBegin{mesh_chord_error} maximum chord approximation error used for meshing.
 *  \cgalParamEnd
 *  \cgalParamBegin{face_proxy_map} a ReadWritePropertyMap with
 * `boost::graph_traits<TriangleMesh>::%face_descriptor` as key and `std::size_t` as value type.
 * A proxy is a set of connected facets which are placed under the same proxy patch (see \cgalFigureRef{iterations}).
 * The proxy-id is contiguous in range [0, number_of_proxies - 1].
 *  \cgalParamEnd
 *  \cgalParamBegin{proxies} output iterator over proxies.
 *  \cgalParamEnd
 *  \cgalParamBegin{anchor_points} output iterator over anchor points. 
 *  \cgalParamEnd
 *  \cgalParamBegin{indexed_triangles} output iterator over indexed triangles.
 *  \cgalParamEnd
 * \cgalNamedParamsEnd
 */
template <typename TriangleMesh, typename NamedParameters>
bool mesh_approximation(const TriangleMesh &tm, const NamedParameters &np)
{
  using boost::get_param;
  using boost::choose_param;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;
  typedef typename Geom_traits::FT FT;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type Vertex_point_map;
  Vertex_point_map point_pmap = choose_param(get_param(np, internal_np::vertex_point),
    get_property_map(vertex_point, const_cast<TriangleMesh &>(tm)));

  typedef CGAL::Approximation_l21_traits<TriangleMesh, Vertex_point_map, false, Geom_traits> Approximation_traits;
  typedef CGAL::VSA_approximation<TriangleMesh, Vertex_point_map> L21_approx;

  L21_approx approx(tm, point_pmap);
  Approximation_traits l21_metric(tm, point_pmap);
  approx.set_metric(l21_metric);

  // default hierarchical seeding
  CGAL::Approximation_seeding_tag method = choose_param(
    get_param(np, internal_np::seeding_method), CGAL::Hierarchical);
  boost::optional<std::size_t> max_nb_proxies = choose_param(
    get_param(np, internal_np::max_nb_proxies), boost::optional<std::size_t>());
  boost::optional<FT> min_error_drop = choose_param(
    get_param(np, internal_np::min_error_drop), boost::optional<FT>());
  std::size_t nb_of_relaxations = choose_param(get_param(np, internal_np::nb_of_relaxations), 5);
  approx.seeding(method, max_nb_proxies, min_error_drop, nb_of_relaxations);

// reasonable default number of iterations
  std::size_t nb_of_iterations_default = max_nb_proxies ? num_faces(tm) / *max_nb_proxies : 30;
  nb_of_iterations_default = (std::min)((std::max)(
    nb_of_iterations_default, static_cast<std::size_t>(20)), static_cast<std::size_t>(60));
  const std::size_t nb_of_iterations = choose_param(
    get_param(np, internal_np::nb_of_iterations), nb_of_iterations_default);

  approx.run(nb_of_iterations);

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  std::cout << "#px = " << approx.proxies_size()
    << ", #itr = " << nb_of_iterations
    << ", #relx = " << nb_of_relaxations << std::endl;
#endif

  // get proxy map
  typedef typename boost::lookup_named_param_def<
    internal_np::facet_proxy_map_t,
    NamedParameters,
    internal_np::vsa_no_output_t>::type FPMap;
  FPMap fproxymap = choose_param(
    get_param(np, internal_np::facet_proxy_map), internal_np::vsa_no_output);
  get_proxy_map(approx, fproxymap);

  // get proxies
  typedef typename boost::lookup_named_param_def <
    internal_np::proxies_t,
    NamedParameters,
    internal_np::vsa_no_output_t>::type ProxiesOutItr;
  ProxiesOutItr pxies_out_itr = choose_param(
    get_param(np, internal_np::proxies), internal_np::vsa_no_output);
  get_proxies(approx, pxies_out_itr);

  // meshing
  const FT chord_error = choose_param(get_param(np, internal_np::mesh_chord_error), FT(5.0));
  const bool is_manifold = approx.extract_mesh(chord_error);

  // get anchor points
  typedef typename boost::lookup_named_param_def<
    internal_np::anchor_points_t,
    NamedParameters,
    internal_np::vsa_no_output_t>::type Anchor_point_output_iterator;
  Anchor_point_output_iterator apts_out_itr = choose_param(
    get_param(np, internal_np::anchor_points) , internal_np::vsa_no_output);
  get_anchor_points(approx, apts_out_itr);

  // get indexed triangles
  typedef typename boost::lookup_named_param_def<
    internal_np::indexed_triangles_t,
    NamedParameters,
    internal_np::vsa_no_output_t>::type Indexed_triangles_output_iterator;
  Indexed_triangles_output_iterator tris_out_itr = choose_param(
    get_param(np, internal_np::indexed_triangles) , internal_np::vsa_no_output);
  get_indexed_triangles(approx, tris_out_itr);

  return is_manifold;
}

} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
