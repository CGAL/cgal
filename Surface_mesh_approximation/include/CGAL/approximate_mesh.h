#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H

#include <CGAL/license/Surface_mesh_approximation.h>


#include <CGAL/Variational_shape_approximation.h>
#include <CGAL/internal/Surface_mesh_approximation/named_function_params.h>
#include <CGAL/internal/Surface_mesh_approximation/named_params_helper.h>

#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/optional.hpp>
#include <boost/type_traits/is_same.hpp>

#include <iostream>

namespace CGAL {
namespace VSA {

/// \ingroup PkgTSMA
/// @brief Verbose level enumeration.
enum Verbose_level {
  /// Silent
  Silent,
  /// Main steps
  Main_steps,
  /// Verbose
  Verbose
};

// the named parameter header being not documented the doc is put here for now
#ifdef DOXYGEN_RUNNING
namespace parameters {

/*! \ingroup vsa_namedparameters
 * This function is used when default parameters are just fine for approximation or meshing.
 */
unspecified_type all_default();

} // namespace parameters
#endif

/*!
 * \ingroup PkgTSMA
 * @brief Approximates the input mesh with plane proxies.
 * This function uses the Variational Shape Approximation algorithm described in \cgalCite{cgal:cad-vsa-04}
 * to approximate a triangle surface mesh, with indexed triangles as output.
 *
 * @tparam TriangleMesh model of `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref vsa_namedparameters
 *
 * @param tm triangle surface mesh to be approximated
 * @param np optional sequence of \ref vsa_namedparameters among the ones listed below
 * @return `true` if the indexed triangles represent a 2-manifold, oriented surface mesh, and `false` otherwise. 
 *
 * \cgalNamedParamsBegin
 *  \cgalParamBegin{geom_traits} a geometric traits class instance, model of Kernel.
 *    Exact constructions kernels are not supported by this function.
 *  \cgalParamEnd
 *  \cgalParamBegin{vertex_point_map} the property map with the points associated
 *    to the vertices of `tm`. Instance of a class model of `ReadablePropertyMap`.
 *  \cgalParamEnd
 *  \cgalParamBegin{verbose_level} set verbose level.
 *  \cgalParamEnd
 *  \cgalParamBegin{seeding_method} selection of seeding method.
 *  \cgalParamEnd
 *  \cgalParamBegin{max_nb_of_proxies} maximum number of proxies to approximate the input mesh.
 *  \cgalParamEnd
 *  \cgalParamBegin{min_error_drop} minimum error drop of the approximation, expressed in ratio between two iterations of proxy addition.
 *  \cgalParamEnd
 *  \cgalParamBegin{nb_of_iterations} number of partitioning and fitting iterations after seeding.
 *  \cgalParamEnd
 *  \cgalParamBegin{nb_of_relaxations} number of relaxation iterations interleaved within seeding.
 *  \cgalParamEnd
 *  \cgalParamBegin{subdivision_ratio} chord subdivision ratio threshold to the chord length or average edge length.
 *  \cgalParamEnd
 *  \cgalParamBegin{relative_to_chord} if `true` the `subdivision_ratio` is the ratio of the
 *    furthest vertex distance to the chord length, otherwise is the average edge length.
 *  \cgalParamEnd
 *  \cgalParamBegin{with_dihedral_angle} if `true` the `subdivision_ratio` is weighted by dihedral angle, `false` otherwise.
 *  \cgalParamEnd
 *  \cgalParamBegin{optimize_anchor_location} if `true` optimize the anchor locations, `false` otherwise.
 *  \cgalParamEnd
 *  \cgalParamBegin{pca_plane} set `true` if use PCA plane fitting, otherwise use the default area averaged plane parameters.
 *  \cgalParamEnd
 *  \cgalParamBegin{face_proxy_map} a ReadWritePropertyMap with
 * `boost::graph_traits<TriangleMesh>::%face_descriptor` as key and `std::size_t` as value type.
 * A proxy is a set of connected faces which are placed under the same proxy patch (see \cgalFigureRef{iterations}).
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
bool approximate_mesh(const TriangleMesh &tm, const NamedParameters &np)
{
  using boost::get_param;
  using boost::choose_param;
  namespace vsa_np = CGAL::VSA::internal_np;

  typedef typename vsa_np::GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;
  typedef typename Geom_traits::FT FT;

  typedef typename vsa_np::GetVertexPointMap<TriangleMesh, NamedParameters>::type Vertex_point_map;
  Vertex_point_map point_pmap = choose_param(get_param(np, vsa_np::vertex_point),
    get_property_map(vertex_point, const_cast<TriangleMesh &>(tm)));

  typedef CGAL::Variational_shape_approximation<TriangleMesh, Vertex_point_map> L21_approx;
  typedef typename L21_approx::Error_metric L21_metric;

  const Verbose_level vl = choose_param(
    get_param(np, vsa_np::verbose_level), CGAL::VSA::Silent);

  if (vl == CGAL::VSA::Main_steps || vl == CGAL::VSA::Verbose) {
    std::cout << "Variational shape approximation:"
      << "\n#f " << num_faces(tm)
      << "\n#v " << num_vertices(tm) << std::endl;
  }

  L21_approx approx(tm, point_pmap);
  L21_metric metric(tm, point_pmap);
  approx.set_metric(metric);

  // hierarchical seeding by default
  CGAL::VSA::Seeding_method method = choose_param(
    get_param(np, vsa_np::seeding_method), CGAL::VSA::Hierarchical);
  boost::optional<std::size_t> max_nb_of_proxies = choose_param(
    get_param(np, vsa_np::max_nb_of_proxies), boost::optional<std::size_t>());
  boost::optional<FT> min_error_drop = choose_param(
    get_param(np, vsa_np::min_error_drop), boost::optional<FT>());
  std::size_t nb_of_relaxations = choose_param(get_param(np, vsa_np::nb_of_relaxations), 5);

  if (vl == CGAL::VSA::Verbose) {
    std::cout << (method == CGAL::VSA::Random ? "Random" :
      (method == CGAL::VSA::Incremental ? "Incremental" : "Hierarchical")) << " seeding.";
    std::cout << "\n#max_nb_of_proxies = " << *max_nb_of_proxies
      << "\n#min_error_drop = " << *min_error_drop
      << "\nnb_of_relaxations " << nb_of_relaxations << std::endl;
  }

  approx.initialize_seeds(method, max_nb_of_proxies, min_error_drop, nb_of_relaxations);

  if (vl == CGAL::VSA::Main_steps || vl == CGAL::VSA::Verbose)
    std::cout << "Seeding done." << std::endl;

  // default number of iterations
  std::size_t nb_of_iterations_default = max_nb_of_proxies ? num_faces(tm) / *max_nb_of_proxies : 30;
  nb_of_iterations_default = (std::min)((std::max)(
    nb_of_iterations_default, static_cast<std::size_t>(20)), static_cast<std::size_t>(60));
  const std::size_t nb_of_iterations = choose_param(
    get_param(np, vsa_np::nb_of_iterations), nb_of_iterations_default);

  if (vl == CGAL::VSA::Verbose)
    std::cout << "\n#nb_of_iterations = " << nb_of_iterations << std::endl;

  approx.run(nb_of_iterations);

  if (vl == CGAL::VSA::Main_steps || vl == CGAL::VSA::Verbose) {
    std::cout << "Approximation done."
      << "\n#proxies = " << approx.proxies_size() << std::endl;
  }

  // get proxy map
  typedef typename boost::lookup_named_param_def<
    vsa_np::face_proxy_map_t,
    NamedParameters,
    vsa_np::dummy_output_t>::type Face_proxy_map;
  Face_proxy_map fproxymap = choose_param(
    get_param(np, vsa_np::face_proxy_map), vsa_np::dummy_output);
  vsa_np::face_proxy_map_helper(approx, fproxymap);

  if (!boost::is_same<Face_proxy_map, vsa_np::dummy_output_t>::value
    && (vl == CGAL::VSA::Main_steps || vl == CGAL::VSA::Verbose))
    std::cout << "Filling face proxy map done." << std::endl;

  // get proxies
  typedef typename boost::lookup_named_param_def<
    vsa_np::proxies_t,
    NamedParameters,
    vsa_np::dummy_output_t>::type Proxies_output_iterator;
  Proxies_output_iterator pxies_out_itr = choose_param(
    get_param(np, vsa_np::proxies), vsa_np::dummy_output);
  vsa_np::proxies_helper(approx, pxies_out_itr);

  if (!boost::is_same<Proxies_output_iterator, vsa_np::dummy_output_t>::value
    && (vl == CGAL::VSA::Main_steps || vl == CGAL::VSA::Verbose))
    std::cout << "Get proxies done." << std::endl;

  // meshing
  typedef typename boost::lookup_named_param_def<
    vsa_np::anchors_t,
    NamedParameters,
    vsa_np::dummy_output_t>::type Anchors_output_iterator;
  typedef typename boost::lookup_named_param_def<
    vsa_np::triangles_t,
    NamedParameters,
    vsa_np::dummy_output_t>::type Triangles_output_iterator;

  bool is_manifold = false;
  if (!boost::is_same<Anchors_output_iterator, vsa_np::dummy_output_t>::value
    || !boost::is_same<Triangles_output_iterator, vsa_np::dummy_output_t>::value) {
    if (vl == CGAL::VSA::Verbose) {
      const FT subdivision_ratio = choose_param(get_param(np, vsa_np::subdivision_ratio), FT(5.0));
      const bool relative_to_chord = choose_param(get_param(np, vsa_np::relative_to_chord), false);
      const bool with_dihedral_angle = choose_param(get_param(np, vsa_np::with_dihedral_angle), false);
      const bool optimize_anchor_location = choose_param(get_param(np, vsa_np::optimize_anchor_location), true);
      const bool pca_plane = choose_param(get_param(np, vsa_np::pca_plane), false);
      std::cout << "Meshing: "
        << "\nchord_error = " << subdivision_ratio
        << "\nrelative_to_chord = " << relative_to_chord
        << "\nwith_dihedral_angle = " << with_dihedral_angle
        << "\noptimize_anchor_location = " << optimize_anchor_location
        << "\npca_plane = " << pca_plane << std::endl;
    }

    is_manifold = approx.extract_mesh(np);

    if (vl == CGAL::VSA::Main_steps || vl == CGAL::VSA::Verbose)
      std::cout << "Meshing done.\n"
        << (is_manifold ? "Can" : "Cannot") << " be built into 2-manifold surface." << std::endl;
  }

  // get anchor points
  Anchors_output_iterator apts_out_itr = choose_param(
    get_param(np, vsa_np::anchors) , vsa_np::dummy_output);
  vsa_np::anchors_helper(approx, apts_out_itr);

  if (!boost::is_same<Anchors_output_iterator, vsa_np::dummy_output_t>::value
    && (vl == CGAL::VSA::Main_steps || vl == CGAL::VSA::Verbose))
    std::cout << "Get anchors done." << std::endl;

  // get indexed triangles
  Triangles_output_iterator tris_out_itr = choose_param(
    get_param(np, vsa_np::triangles) , vsa_np::dummy_output);
  vsa_np::triangles_helper(approx, tris_out_itr);

  if (!boost::is_same<Triangles_output_iterator, vsa_np::dummy_output_t>::value
    && (vl == CGAL::VSA::Main_steps || vl == CGAL::VSA::Verbose))
    std::cout << "Get indexed triangles done." << std::endl;

  return is_manifold;
}

} // namespace VSA
} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
