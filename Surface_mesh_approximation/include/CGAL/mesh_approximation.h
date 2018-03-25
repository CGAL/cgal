#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H

#include <CGAL/license/Surface_mesh_approximation.h>


#include <CGAL/VSA_approximation.h>
#include <CGAL/internal/Surface_mesh_approximation/named_function_params.h>
#include <CGAL/internal/Surface_mesh_approximation/named_params_helper.h>

#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/optional.hpp>
#include <boost/type_traits/is_same.hpp>

#include <iostream>

namespace CGAL {

/// \ingroup PkgTSMA
/// @brief Verbose level enumeration for Variational Shape Approximation algorithm.
enum Approximation_verbose_level {
  /// Silent
  Silent,
  /// Main steps
  Main_steps,
  /// Verbose
  Verbose
};

/*!
 * \ingroup PkgTSMA
 * @brief Approximates the input mesh by fitting it with proxies.
 * This function uses the Variational Shape Approximation algorithm 
 * to approximate a triangle surface mesh, with indexed triangles as output.
 *
 * @tparam TriangleMesh model of `FaceListGraph`.
 *         If `TriangleMesh` has an internal property map for `CGAL::face_index_t`,
 *         and no `face_index_map` is given as a named parameter, then the internal one should be initialized.
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param tm triangle surface mesh to be approximated
 * @param np optional sequence of \ref namedparameters among the ones listed below
 * @return `true` if the indexed triangles represent a 2-manifold, oriented surface mesh, and `false` otherwise. 
 *
 * \cgalNamedParamsBegin
 *  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
 *    Exact constructions kernels are not supported by this function.
 *  \cgalParamEnd
 *  \cgalParamBegin{vertex_point_map} the property map with the points associated
 *    to the vertices of `tm`. Instance of a class model of `ReadablePropertyMap`.
 *  \cgalParamEnd
 *  \cgalParamBegin{seeding_method} selection of seeding method.
 *  \cgalParamEnd
 *  \cgalParamBegin{max_nb_proxies} maximum number of proxies to approximate the input mesh.
 *  \cgalParamEnd
 *  \cgalParamBegin{min_error_drop} minimum error drop of the approximation, expressed in ratio between two iterations of proxy addition.
 *  \cgalParamEnd
 *  \cgalParamBegin{nb_of_iterations} number of partitioning and fitting iterations after seeding.
 *  \cgalParamEnd
 *  \cgalParamBegin{nb_of_relaxations} number of relaxation iterations interleaved within seeding.
 *  \cgalParamEnd
 *  \cgalParamBegin{subdivision_ratio} chord subdivision ratio threshold to the chord length or average edge length.
 *  \cgalParamEnd
 *  \cgalParamBegin{relative_to_chord} set `true` if the subdivision_ratio is the ratio of the
 *    furthest vertex distance to the chord length, or to the average edge length otherwise.
 *  \cgalParamEnd
 *  \cgalParamBegin{with_dihedral_angle}  set `true` if subdivision_ratio is weighted by dihedral angle, `false` otherwise.
 *  \cgalParamEnd
 *  \cgalParamBegin{optimize_anchor_location}  set `true` if optimize the anchor locations, `false` otherwise.
 *  \cgalParamEnd
 *  \cgalParamBegin{pca_plane}  set `true` if use PCA plane fitting, otherwise use the default area averaged plane parameters.
 *  \cgalParamEnd
 *  \cgalParamBegin{facet_proxy_map} a ReadWritePropertyMap with
 * `boost::graph_traits<TriangleMesh>::%face_descriptor` as key and `std::size_t` as value type.
 * A proxy is a set of connected facets which are placed under the same proxy patch (see \cgalFigureRef{iterations}).
 * The proxy-ids are contiguous in range [0, number_of_proxies - 1].
 *  \cgalParamEnd
 *  \cgalParamBegin{proxies} output iterator over proxies.
 *  \cgalParamEnd
 *  \cgalParamBegin{anchors} output iterator over anchor points. 
 *  \cgalParamEnd
 *  \cgalParamBegin{triangles} output iterator over indexed triangles.
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

  typedef CGAL::VSA::L21_metric_vector_proxy_no_area_weighting<TriangleMesh, Vertex_point_map, Geom_traits> L21_metric;
  typedef CGAL::VSA_approximation<TriangleMesh, Vertex_point_map> L21_approx;
  typedef L21_approx::Error_metric L21_metric;

  const Approximation_verbose_level vl = choose_param(
    get_param(np, internal_np::verbose_level), CGAL::Main_steps);

  if (vl == CGAL::Main_steps || vl == CGAL::Verbose) {
    std::cout << "Variational shape approximation:"
      << "\n#f " << num_faces(tm)
      << "\n#v " << num_vertices(tm) << std::endl;
  }

  L21_approx approx(tm, point_pmap);
  L21_metric metric(tm, point_pmap);
  approx.set_metric(metric);

  // hierarchical seeding by default
  CGAL::Approximation_seeding_tag method = choose_param(
    get_param(np, internal_np::seeding_method), CGAL::Hierarchical);
  boost::optional<std::size_t> max_nb_proxies = choose_param(
    get_param(np, internal_np::max_nb_proxies), boost::optional<std::size_t>());
  boost::optional<FT> min_error_drop = choose_param(
    get_param(np, internal_np::min_error_drop), boost::optional<FT>());
  std::size_t nb_of_relaxations = choose_param(get_param(np, internal_np::nb_of_relaxations), 5);

  if (vl == CGAL::Verbose) {
    std::cout << (method == CGAL::Random ? "Random" :
      (method == CGAL::Incremental ? "Incremental" : "Hierarchical")) << " seeding.";
    std::cout << "\n#max_nb_proxies = " << *max_nb_proxies
      << "\n#min_error_drop = " << *min_error_drop
      << "\nnb_of_relaxations " << nb_of_relaxations << std::endl;
  }

  approx.seeding(method, max_nb_proxies, min_error_drop, nb_of_relaxations);

  if (vl == CGAL::Main_steps || vl == CGAL::Verbose)
    std::cout << "Seeding done." << std::endl;

  // default number of iterations
  std::size_t nb_of_iterations_default = max_nb_proxies ? num_faces(tm) / *max_nb_proxies : 30;
  nb_of_iterations_default = (std::min)((std::max)(
    nb_of_iterations_default, static_cast<std::size_t>(20)), static_cast<std::size_t>(60));
  const std::size_t nb_of_iterations = choose_param(
    get_param(np, internal_np::nb_of_iterations), nb_of_iterations_default);

  if (vl == CGAL::Verbose)
    std::cout << "\n#nb_of_iterations = " << nb_of_iterations << std::endl;

  approx.run(nb_of_iterations);

  if (vl == CGAL::Main_steps || vl == CGAL::Verbose) {
    std::cout << "Approximation done."
      << "\n#proxies = " << approx.proxies_size() << std::endl;
  }

  // get proxy map
  typedef typename boost::lookup_named_param_def<
    internal_np::facet_proxy_map_t,
    NamedParameters,
    internal_np::vsa_no_output_t>::type Face_proxy_map;
  Face_proxy_map fproxymap = choose_param(
    get_param(np, internal_np::facet_proxy_map), internal_np::vsa_no_output);
  facet_proxy_map(approx, fproxymap);

  if (!boost::is_same<Face_proxy_map, internal_np::vsa_no_output_t>::value
    && (vl == CGAL::Main_steps || vl == CGAL::Verbose))
    std::cout << "Filling facet proxy map done." << std::endl;

  // get proxies
  typedef typename boost::lookup_named_param_def<
    internal_np::proxies_t,
    NamedParameters,
    internal_np::vsa_no_output_t>::type Proxies_output_iterator;
  Proxies_output_iterator pxies_out_itr = choose_param(
    get_param(np, internal_np::proxies), internal_np::vsa_no_output);
  proxies(approx, pxies_out_itr);

  if (!boost::is_same<Proxies_output_iterator, internal_np::vsa_no_output_t>::value
    && (vl == CGAL::Main_steps || vl == CGAL::Verbose))
    std::cout << "Get proxies done." << std::endl;

  // meshing
  typedef typename boost::lookup_named_param_def<
    internal_np::anchors_t,
    NamedParameters,
    internal_np::vsa_no_output_t>::type Anchors_output_iterator;
  typedef typename boost::lookup_named_param_def<
    internal_np::triangles_t,
    NamedParameters,
    internal_np::vsa_no_output_t>::type Triangles_output_iterator;

  bool is_manifold = false;
  if (!boost::is_same<Anchors_output_iterator, internal_np::vsa_no_output_t>::value
    || !boost::is_same<Triangles_output_iterator, internal_np::vsa_no_output_t>::value) {
    if (vl == CGAL::Verbose) {
      const FT subdivision_ratio = choose_param(get_param(np, internal_np::subdivision_ratio), FT(5.0));
      const bool relative_to_chord = choose_param(get_param(np, internal_np::relative_to_chord), false);
      const bool with_dihedral_angle = choose_param(get_param(np, internal_np::with_dihedral_angle), false);
      const bool optimize_anchor_location = choose_param(get_param(np, internal_np::optimize_anchor_location), true);
      const bool pca_plane = choose_param(get_param(np, internal_np::pca_plane), false);
      std::cout << "Meshing: "
        << "\nchord_error = " << subdivision_ratio
        << "\nrelative_to_chord = " << relative_to_chord
        << "\nwith_dihedral_angle = " << with_dihedral_angle
        << "\noptimize_anchor_location = " << optimize_anchor_location
        << "\npca_plane = " << pca_plane << std::endl;
    }

    is_manifold = approx.extract_mesh(np);

    if (vl == CGAL::Main_steps || vl == CGAL::Verbose)
      std::cout << "Meshing done.\n"
        << (is_manifold ? "Can" : "Cannot") << " be built into 2-manifold surface." << std::endl;
  }

  // get anchor points
  Anchors_output_iterator apts_out_itr = choose_param(
    get_param(np, internal_np::anchors) , internal_np::vsa_no_output);
  anchors(approx, apts_out_itr);

  if (!boost::is_same<Anchors_output_iterator, internal_np::vsa_no_output_t>::value
    && (vl == CGAL::Main_steps || vl == CGAL::Verbose))
    std::cout << "Get anchors done." << std::endl;

  // get indexed triangles
  Triangles_output_iterator tris_out_itr = choose_param(
    get_param(np, internal_np::triangles) , internal_np::vsa_no_output);
  triangles(approx, tris_out_itr);

  if (!boost::is_same<Triangles_output_iterator, internal_np::vsa_no_output_t>::value
    && (vl == CGAL::Main_steps || vl == CGAL::Verbose))
    std::cout << "Get indexed triangles done." << std::endl;

  return is_manifold;
}

} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
