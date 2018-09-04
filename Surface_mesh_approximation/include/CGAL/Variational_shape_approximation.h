#ifndef CGAL_VARIATIONAL_SHAPE_APPROXIMATION_H
#define CGAL_VARIATIONAL_SHAPE_APPROXIMATION_H

#include <CGAL/license/Surface_mesh_approximation.h>


#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/array.h>

#include <CGAL/VSA/L21_metric_plane_proxy.h>
#include <CGAL/Default.h>
#include <CGAL/tags.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>

#include <CGAL/internal/Surface_mesh_approximation/named_function_params.h>
#include <CGAL/internal/Surface_mesh_approximation/named_params_helper.h>

#include <vector>
#include <stack>
#include <queue>
#include <iterator>
#include <cmath>
#include <cstdlib>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#endif // CGAL_LINKED_WITH_TBB

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
#include <iostream>
#endif

#define CGAL_VSA_INVALID_TAG std::numeric_limits<std::size_t>::max()

namespace CGAL {

namespace VSA {
/// \ingroup PkgTSMA
/// @brief Seeding method enumeration for Variational Shape Approximation algorithm.
enum Seeding_method {
  /// Random seeding
  Random,
  /// Incremental seeding
  Incremental,
  /// Hierarchical seeding
  Hierarchical
};
} // namespace VSA

/// \ingroup PkgTSMA
/// @brief Main class for Variational Shape Approximation algorithm described in \cgalCite{cgal:cad-vsa-04}.
/// @tparam TriangleMesh a model of `FaceListGraph`
/// @tparam VertexPointMap a `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key and `GeomTraits::Point_3` as value type
/// @tparam ErrorMetricProxy a model of `ErrorMetricProxy`
/// @tparam GeomTraits a model of Kernel
/// @tparam Concurrency_tag concurrency tag.
template <typename TriangleMesh,
  typename VertexPointMap,
  typename ErrorMetricProxy = CGAL::Default,
  typename GeomTraits = CGAL::Default,
  typename Concurrency_tag = CGAL::Sequential_tag>
class Variational_shape_approximation {
// public typedefs
public:

  /// \name Types
  /// @{
#ifndef DOXYGEN_RUNNING
  // GeomTraits type
  typedef typename CGAL::Default::Get<
    GeomTraits,
    typename Kernel_traits<
      typename boost::property_traits<VertexPointMap>::value_type
    >::Kernel >::type Geom_traits;
#else
  typedef GeomTraits Geom_traits;
#endif
  // ErrorMetricProxy type
#ifndef DOXYGEN_RUNNING
  typedef typename CGAL::Default::Get<ErrorMetricProxy,
    CGAL::VSA::L21_metric_plane_proxy<TriangleMesh, VertexPointMap, Geom_traits> >::type Error_metric;
#else
  typedef ErrorMetricProxy Error_metric;
#endif
  // Proxy type
  typedef typename Error_metric::Proxy Proxy;

  // indexed triangle type
  typedef CGAL::cpp11::array<std::size_t, 3> Indexed_triangle;
  /// @}

// private typedefs and data member
private:
  // Geom_traits typedefs
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point_3;
  typedef typename Geom_traits::Vector_3 Vector_3;
  typedef typename Geom_traits::Plane_3 Plane_3;
  typedef typename Geom_traits::Construct_vector_3 Construct_vector_3;
  typedef typename Geom_traits::Construct_point_3 Construct_point_3;
  typedef typename Geom_traits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename Geom_traits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename Geom_traits::Compute_scalar_product_3 Compute_scalar_product_3;
  typedef typename Geom_traits::Construct_translated_point_3 Construct_translated_point_3;

  // graph_traits typedefs
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  // internal typedefs
  typedef CGAL::internal::vertex_property_t<std::size_t> Vertex_anchor_tag;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Vertex_anchor_tag >::type Vertex_anchor_map;

  typedef CGAL::internal::face_property_t<std::size_t> Face_proxy_tag;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_proxy_tag >::type Face_proxy_map;

  typedef std::vector<halfedge_descriptor> Boundary_chord;
  typedef typename Boundary_chord::iterator Boundary_chord_iterator;

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
public:
#endif
  // The proxy wrapper for approximation.
  struct Proxy_wrapper {
    Proxy_wrapper(const Proxy &p, const std::size_t &i, const face_descriptor s, const FT &e)
      : px(p), idx(i), seed(s), err(e) {}

    Proxy px; // parameterized proxy
    std::size_t idx; // proxy index, maintained to be the same as its position in proxies vector
    face_descriptor seed; // proxy seed
    FT err; // proxy fitting error
  };
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
private:
#endif

  // The proxy fitting plane for meshing.
  struct Proxy_plane {
    Proxy_plane(const Plane_3 &p, const Vector_3 &n, const FT &a)
      : plane(p), normal(n), area(a) {}

    Plane_3 plane;
    Vector_3 normal;
    FT area;
  };

  // The face candidate to be queued.
  struct Face_to_integrate {
    Face_to_integrate(const face_descriptor f_, const std::size_t &px_, const FT &err_)
      : f(f_), px(px_), err(err_) {}

    bool operator<(const Face_to_integrate &rhs) const {
      return err > rhs.err;
    }

    face_descriptor f; // face
    std::size_t px; // proxy index
    FT err; // fitting error
  };

  // Proxy error with its index.
  struct Proxy_error {
    Proxy_error(const std::size_t &px_, const FT &err_)
      : px(px_), err(err_) {}

    // in ascending order
    bool operator<(const Proxy_error &rhs) const {
      return err < rhs.err;
    }

    std::size_t px;
    FT err;
  };

  // The anchor attached to a vertex.
  struct Anchor {
    Anchor(const vertex_descriptor &vtx_, const Point_3 pos_)
      : vtx(vtx_), pos(pos_) {}

    vertex_descriptor vtx; // The associated vertex.
    Point_3 pos; // The position of the anchor.
  };

  // The boundary cycle of a region.
  // One region may have multiple boundary cycles.
  struct Boundary_cycle {
    Boundary_cycle(const halfedge_descriptor &h)
      : he_head(h), num_anchors(0) {}

    halfedge_descriptor he_head; // Heading halfedge of the boundary cycle.
    std::size_t num_anchors; // Number of anchors on the boundary cycle.
  };

  // Triangle polyhedron builder.
  template <typename HDS>
  class Triangle_polyhedron_builder : public CGAL::Modifier_base<HDS> {
    const std::vector<Point_3> &vtxs;
    const std::vector<Indexed_triangle> &tris;
  public:
    bool is_manifold;
    Triangle_polyhedron_builder(const std::vector<Point_3> &vtxs_,
      const std::vector<Indexed_triangle> &tris_)
      : vtxs(vtxs_), tris(tris_), is_manifold(true) {}

    void operator()(HDS &hds) {
      CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);
      typedef typename HDS::Vertex Vertex;
      typedef typename Vertex::Point Point;
      builder.begin_surface(vtxs.size(), tris.size());
      BOOST_FOREACH(const Point_3 &v, vtxs)
        builder.add_vertex(Point(v));
      BOOST_FOREACH(const Indexed_triangle &t, tris) {
        if (builder.test_facet(t.begin(), t.end())) {
          builder.begin_facet();
          builder.add_vertex_to_facet(t[0]);
          builder.add_vertex_to_facet(t[1]);
          builder.add_vertex_to_facet(t[2]);
          builder.end_facet();
        }
        else
          is_manifold = false;
      }
      builder.end_surface();
    }
  };

  // member variables
  // The triangle mesh.
  const TriangleMesh *m_ptm;
  // The mesh vertex point map.
  VertexPointMap m_vpoint_map;
  // The approximation object.
  const Error_metric *m_metric;

  Construct_vector_3 vector_functor;
  Construct_point_3 point_functor;
  Construct_scaled_vector_3 scale_functor;
  Construct_sum_of_vectors_3 sum_functor;
  Compute_scalar_product_3 scalar_product_functor;
  Construct_translated_point_3 translate_point_functor;

  // The face proxy index map.
  Face_proxy_map m_fproxy_map;
  // The attached anchor index of a vertex.
  Vertex_anchor_map m_vanchor_map;

  // Proxies.
  std::vector<Proxy_wrapper> m_proxies;
  // Proxy planes
  std::vector<Proxy_plane> m_px_planes;

  // All anchors.
  std::vector<Anchor> m_anchors;
  // All boundary cycles.
  std::vector<Boundary_cycle> m_bcycles;
  // The indexed triangle approximation.
  std::vector<Indexed_triangle> m_tris;

  // meshing parameters
  FT m_average_edge_length;

//member functions
public:
  /// \name Construction
  /// @{
  /*!
   * @brief Initializes internal data for the approximation.
   * @param tm `CGAL TriangleMesh` on which approximation operates
   * @param vpoint_map vertex point map of the mesh
   * @param error_metric an `ErrorMetricProxy` object
   */
  Variational_shape_approximation(const TriangleMesh &tm,
    const VertexPointMap &vpoint_map,
    const Error_metric &error_metric) :
    m_ptm(&tm),
    m_vpoint_map(vpoint_map),
    m_metric(&error_metric),
    m_average_edge_length(0.0) {

    Geom_traits traits;
    vector_functor = traits.construct_vector_3_object();
    point_functor = traits.construct_point_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scalar_product_functor = traits.compute_scalar_product_3_object();
    translate_point_functor = traits.construct_translated_point_3_object();

    m_vanchor_map = CGAL::internal::add_property(
      Vertex_anchor_tag("VSA-vertex_anchor"), *(const_cast<TriangleMesh *>(m_ptm)));

    m_fproxy_map = CGAL::internal::add_property(
      Face_proxy_tag("VSA-face_proxy"), *(const_cast<TriangleMesh *>(m_ptm)));
  }
  /// @}
  
  ~Variational_shape_approximation() {
    if (m_ptm) {
      CGAL::internal::remove_property(m_vanchor_map, *(const_cast<TriangleMesh *>(m_ptm)));
      CGAL::internal::remove_property(m_fproxy_map, *(const_cast<TriangleMesh *>(m_ptm)));
    }
  }

  /// \name Approximation
  /// @{
  /*!
   * @brief Initializes the seeds with both maximum number of proxies and minimum error drop stop criteria.
   * The first criterion met stops the seeding.
   * Parameters out of range are ignored.
   * @tparam NamedParameters a sequence of \ref vsa_namedparameters

   * @param np optional sequence of \ref vsa_namedparameters among the ones listed below
   * @return number of proxies initialized

   * \cgalNamedParamsBegin{Seeding Named Parameters}
   *  \cgalParamBegin{seeding_method} selection of seeding method.
   *  \cgalParamEnd
   *  \cgalParamBegin{max_nb_of_proxies} maximum number of proxies to approximate the input mesh.
   *  \cgalParamEnd
   *  \cgalParamBegin{min_error_drop} minimum error drop of the approximation, expressed in ratio between two iterations of proxy addition.
   *  \cgalParamEnd
   *  \cgalParamBegin{nb_of_relaxations} number of relaxation iterations interleaved within seeding.
   *  \cgalParamEnd
   * \cgalNamedParamsEnd
   */
  template <typename NamedParameters>
  std::size_t initialize_seeds(const NamedParameters &np) {
    using boost::get_param;
    using boost::choose_param;
    namespace vsa_np = CGAL::VSA::internal_np;

    const CGAL::VSA::Seeding_method method = choose_param(
      get_param(np, vsa_np::seeding_method), CGAL::VSA::Hierarchical);
    const boost::optional<std::size_t> max_nb_of_proxies = choose_param(
      get_param(np, vsa_np::max_nb_of_proxies), boost::optional<std::size_t>());
    const boost::optional<FT> min_error_drop = choose_param(
      get_param(np, vsa_np::min_error_drop), boost::optional<FT>());
    const std::size_t nb_of_relaxations = choose_param(get_param(np, vsa_np::nb_of_relaxations), 5);

    // maximum number of proxies internally, maybe better choice?
    const std::size_t nb_px = num_faces(*m_ptm) / 3;

    // initialize proxies and the proxy map to prepare for insertion
    bootstrap_from_connected_components();

    if (min_error_drop && *min_error_drop > FT(0.0) && *min_error_drop < FT(1.0)) {
      // as long as minimum error is specified and valid
      // maximum number of proxies always exists, no matter specified or not or out of range
      // there is always a maximum number of proxies explicitly (max_nb_of_proxies) or implicitly (nb_px)
      std::size_t max_nb_px_adjusted = nb_px;
      if (max_nb_of_proxies && *max_nb_of_proxies < nb_px && *max_nb_of_proxies > 0)
        max_nb_px_adjusted = *max_nb_of_proxies;
      switch (method) {
        case VSA::Random:
          return init_random_error(max_nb_px_adjusted, *min_error_drop, nb_of_relaxations);
        case VSA::Incremental:
          return init_incremental_error(max_nb_px_adjusted, *min_error_drop, nb_of_relaxations);
        case VSA::Hierarchical:
          return init_hierarchical_error(max_nb_px_adjusted, *min_error_drop, nb_of_relaxations);
        default:
          return 0;
      }
    }
    else if (max_nb_of_proxies && *max_nb_of_proxies < nb_px && *max_nb_of_proxies > 0) {
      // no valid min_error_drop provided, only max_nb_of_proxies
      switch (method) {
        case VSA::Random:
          return init_random(*max_nb_of_proxies, nb_of_relaxations);
        case VSA::Incremental:
          return init_incremental(*max_nb_of_proxies, nb_of_relaxations);
        case VSA::Hierarchical:
          return init_hierarchical(*max_nb_of_proxies, nb_of_relaxations);
        default:
          return 0;
      }
    }
    else {
      // both parameters are unspecified or out of range
      const FT e(0.1);
      switch (method) {
        case VSA::Random:
          return init_random_error(nb_px, e, nb_of_relaxations);
        case VSA::Incremental:
          return init_incremental_error(nb_px, e, nb_of_relaxations);
        case VSA::Hierarchical:
          return init_hierarchical_error(nb_px, e, nb_of_relaxations);
        default:
          return 0;
      }
    }
  }

  /*!
   * @brief Runs the partitioning and fitting processes on the whole surface.
   * @param nb_of_iterations number of iterations.
   * @return total fitting error
   */
  FT run(std::size_t nb_of_iterations = 1) {
    for (std::size_t i = 0; i < nb_of_iterations; ++i) {
      // tag the whole surface
      BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
        put(m_fproxy_map, f, CGAL_VSA_INVALID_TAG);

      partition(m_proxies.begin(), m_proxies.end());
      fit(m_proxies.begin(), m_proxies.end(), Concurrency_tag());
    }

    return compute_total_error();
  }

  /*!
   * @brief Calls `run` while error decrease is greater than the threshold.
   * @param cvg_threshold the percentage of error change between two successive runs,
   * should be in range (0, 1).
   * @param max_iterations maximum number of iterations allowed
   * @param avg_interval size of error average interval to have smoother convergence curve,
   * if 0 is assigned, 1 is used instead.
   * @return `true` if converged before hitting the maximum iterations, `false` otherwise
   */
  bool run_to_convergence(const FT cvg_threshold,
    const std::size_t max_iterations = 100,
    std::size_t avg_interval = 3) {
    if (avg_interval == 0)
      avg_interval = 1;
    FT drop_pct(0.0);
    FT pre_err = compute_total_error();
    for (std::size_t itr_count = 0; itr_count < max_iterations; itr_count += avg_interval) {
      if (pre_err == FT(0.0))
        return true;

      FT avg_err(0.0);
      for (std::size_t i = 0; i < avg_interval; ++i)
        avg_err += run();
      avg_err /= static_cast<FT>(avg_interval);

      drop_pct = (pre_err - avg_err) / pre_err;
      // the error may fluctuate
      if (drop_pct < FT(0.0))
        drop_pct = -drop_pct;
      if (drop_pct < cvg_threshold)
        return true;

      pre_err = avg_err;
    }

    return false;
  }

  /*!
   * @brief Computes fitting error of current partition to the proxies.
   * @return total fitting error
   */
  FT compute_total_error() {
    FT sum_error(0.0);
    BOOST_FOREACH(const Proxy_wrapper &pxw, m_proxies)
      sum_error += pxw.err;

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
    static std::size_t count = 0;
    std::cerr << '#' << count++ << ": " << sum_error << std::endl;
#endif

    return sum_error;
  }

  /*!
   * @brief Adds proxies to the worst regions one by one.
   * The re-fitting is performed after each proxy is inserted.
   * @param num_proxies number of proxies to be added
   * @param nb_of_iterations number of re-fitting iterations
   * @return number of proxies added
   */
  std::size_t add_to_furthest_proxies(const std::size_t num_proxies,
    const std::size_t nb_of_iterations = 5) {
    std::size_t num_added = 0;
    while (num_added < num_proxies) {
      if (!add_to_furthest_proxy())
        break;
      ++num_added;
      run(nb_of_iterations);
    }
    return num_added;
  }

  /*!
   * @brief Adds proxies by diffusing fitting error into current partition.
   * Each partition is added with the number of proxies in proportion to its fitting error.
   * @param num_proxies number of proxies to be added
   * @return number of proxies successfully added
   */
  std::size_t add_proxies_error_diffusion(const std::size_t num_proxies) {
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
    std::cerr << "#px " << m_proxies.size() << std::endl;
#endif

    const double sum_error = CGAL::to_double(compute_total_error());
    const double avg_error = sum_error / static_cast<double>(num_proxies);

    // number of proxies to be added to each region
    std::vector<std::size_t> num_to_add(m_proxies.size(), 0);
    if (avg_error <= 0.0) {
      // rare case on extremely regular geometry like a cube
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
      std::cerr << "zero error, diffuse w.r.t. number of faces" << std::endl;
#endif
      const double avg_face =
        static_cast<double>(num_faces(*m_ptm)) / static_cast<double>(num_proxies);
      std::vector<double> px_size(m_proxies.size(), 0.0);
      BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
        px_size[get(m_fproxy_map, f)] += 1.0;
      double residual = 0.0;
      for (std::size_t i = 0; i < m_proxies.size(); ++i) {
        const double to_add = (residual + px_size[i]) / avg_face;
        const double to_add_round_up = std::floor(to_add + 0.5);
        residual = (to_add - to_add_round_up) * avg_face;
        num_to_add[i] = static_cast<std::size_t>(to_add_round_up);
      }
    }
    else {
      std::vector<Proxy_error> px_error;
      for (std::size_t i = 0; i < m_proxies.size(); ++i)
        px_error.push_back(Proxy_error(i, m_proxies[i].err));
      // sort partition by error
      std::sort(px_error.begin(), px_error.end());

      // residual from previous proxy in range (-0.5, 0.5] * avg_error
      double residual = 0.0;
      for (std::size_t i = 0; i < m_proxies.size(); ++i) {
        // add error residual from previous proxy
        // to_add maybe negative but greater than -0.5
        const double to_add = (residual + CGAL::to_double(px_error[i].err)) / avg_error;
        const double to_add_round_up = std::floor(to_add + 0.5);
        residual = (to_add - to_add_round_up) * avg_error;
        num_to_add[px_error[i].px] = static_cast<std::size_t>(to_add_round_up);
      }

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
      for (std::size_t i = 0; i < px_error.size(); ++i)
        std::cerr << "#px " << px_error[i].px
          << ", #error " << px_error[i].err
          << ", #num_to_add " << num_to_add[px_error[i].px] << std::endl;
#endif
    }

    std::size_t num_added = 0;
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm)) {
      const std::size_t px_id = get(m_fproxy_map, f);
      if (m_proxies[px_id].seed == f)
        continue;

      if (num_to_add[px_id] > 0) {
        add_one_proxy_at(f);
        --num_to_add[px_id];
        ++num_added;
      }
    }

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
    std::cerr << "#requested/added "
      << num_proxies << '/' << num_added << std::endl;
#endif

    return num_added;
  }
  /// @}

  /// \name Refinement Operations
  /// @{
  /*!
   * @brief Teleports the local minimum to the worst region by combining the merging and adding processes.
   * The re-fitting is performed after each teleportation.
   * Here if we specify more than one proxy this means we teleport in a naive iterative fashion.
   * @param num_proxies number of proxies requested to teleport.
   * @param nb_of_iterations number of re-fitting iterations.
   * @param if_force set `true` to force the teleportation (no merge test).
   * @return number of proxies teleported.
   */
  std::size_t teleport_proxies(const std::size_t num_proxies,
    const std::size_t nb_of_iterations = 5,
    const bool if_force = false) {
    std::size_t num_teleported = 0;
    while (num_teleported < num_proxies) {
      // find worst proxy
      std::size_t px_worst = 0;
      FT max_error = m_proxies.front().err;
      for (std::size_t i = 0; i < m_proxies.size(); ++i) {
        if (max_error < m_proxies[i].err) {
          max_error = m_proxies[i].err;
          px_worst = i;
        }
      }
      bool found = false;
      face_descriptor tele_to;
      BOOST_FOREACH(face_descriptor f, faces(*m_ptm)) {
        if (get(m_fproxy_map, f) == px_worst && f != m_proxies[px_worst].seed) {
          // teleport to anywhere but the seed
          tele_to = f;
          found = true;
          break;
        }
      }
      if (!found)
        return num_teleported;

      // find the best merge pair
      std::size_t px_enlarged = 0, px_merged = 0;
      if (!find_best_merge(px_enlarged, px_merged, !if_force))
        return num_teleported;
      if (px_worst == px_enlarged || px_worst == px_merged)
        return num_teleported;

      // teleport to a face of the worst region
      // update merged proxies
      std::list<face_descriptor> merged_patch;
      BOOST_FOREACH(face_descriptor f, faces(*m_ptm)) {
        std::size_t px_idx = get(m_fproxy_map, f);
        if (px_idx == px_enlarged || px_idx == px_merged) {
          put(m_fproxy_map, f, px_enlarged);
          merged_patch.push_back(f);
        }
      }
      m_proxies[px_enlarged] = fit_proxy_from_patch(merged_patch, px_enlarged);
      // replace the merged proxy position to the newly teleported proxy
      m_proxies[px_merged] = fit_proxy_from_face(tele_to, px_merged);

      num_teleported++;
      // coarse re-fitting
      run(nb_of_iterations);

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
      std::cerr << "teleported" << std::endl;
#endif
    }

    return num_teleported;
  }

  /*!
   * @brief Merges two specified adjacent regions.
   * The overall re-fitting is not performed and the proxy map is maintained.
   * @pre two proxies must be adjacent, and px0 < px1 < proxies.size()
   * @param px0 the enlarged proxy
   * @param px1 the merged proxy
   * @return change of error
   */
  FT merge(const std::size_t px0, const std::size_t px1) {
    // ensure px0 < px1
    if (px0 >= px1 || px1 >= m_proxies.size())
      return FT(0.0);

    const FT pre_err = m_proxies[px0].err + m_proxies[px1].err;
    // merge px1 to px0
    std::list<face_descriptor> merged_patch;
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm)) {
      std::size_t px_idx = get(m_fproxy_map, f);
      if (px_idx == px1 || px_idx == px0) {
        put(m_fproxy_map, f, px0);
        merged_patch.push_back(f);
      }
    }
    m_proxies[px0] = fit_proxy_from_patch(merged_patch, px0);

    // erase px1 and maintain proxy index
    m_proxies.erase(m_proxies.begin() + px1);
    for (std::size_t i = 0; i < m_proxies.size(); ++i)
      m_proxies[i].idx = i;
    // keep face proxy map valid
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm)) {
      if (get(m_fproxy_map, f) > px1)
        put(m_fproxy_map, f, get(m_fproxy_map, f) - 1);
    }

    return m_proxies[px0].err - pre_err;
  }

  /*!
   * @brief Simulates merging and local re-fitting of all pairs of adjacent proxies
   * and finds the best pair to merge.
   * @note The <b>best</b> is defined as the minimum merged sum error
   * <b>change</b> (increase or decrease) among all pairs.
   * @param[out] px_tobe_enlarged proxy index to be enlarged (FIX TERM ENLARGED, NOT DEFINED)
   * @param[out] px_tobe_merged proxy index to be merged,
   * guaranteed to be greater than <em>px_tobe_enlarged</em>.
   * @param if_test set `true` to activate the merge test.
   * The merge test is considered successful if the merged error change
   * is lower than half of the maximum proxy error.
   * @return `true` if best merge pair found, `false` otherwise.
   */
  bool find_best_merge(std::size_t &px_tobe_enlarged,
    std::size_t &px_tobe_merged,
    const bool if_test) {
    typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
    typedef std::pair<std::size_t, std::size_t> ProxyPair;
    typedef std::set<ProxyPair> MergedPair;

    std::vector<std::list<face_descriptor> > px_faces(m_proxies.size());
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
      px_faces[get(m_fproxy_map, f)].push_back(f);

    // find best merge
    MergedPair merged_set;
    FT min_error_change = FT(0.0);
    bool first_merge = true;
    BOOST_FOREACH(edge_descriptor e, edges(*m_ptm)) {
      if (CGAL::is_border(e, *m_ptm))
        continue;
      std::size_t pxi = get(m_fproxy_map, face(halfedge(e, *m_ptm), *m_ptm));
      std::size_t pxj = get(m_fproxy_map, face(opposite(halfedge(e, *m_ptm), *m_ptm), *m_ptm));
      if (pxi == pxj)
        continue;
      if (pxi > pxj)
        std::swap(pxi, pxj);
      if (merged_set.find(ProxyPair(pxi, pxj)) != merged_set.end())
        continue;

      merged_set.insert(ProxyPair(pxi, pxj));
      // simulated merge
      std::list<face_descriptor> merged_patch(px_faces[pxi]);
      BOOST_FOREACH(face_descriptor f, px_faces[pxj])
        merged_patch.push_back(f);
      Proxy_wrapper pxw_tmp = fit_proxy_from_patch(merged_patch, CGAL_VSA_INVALID_TAG);

      const FT error_change = pxw_tmp.err - (m_proxies[pxi].err + m_proxies[pxj].err);
      if (first_merge || error_change < min_error_change) {
        first_merge = false;
        min_error_change = error_change;
        px_tobe_enlarged = pxi;
        px_tobe_merged = pxj;
      }
    }

    if (merged_set.empty())
      return false;

    // test if merge worth it
    if (if_test) {
      FT max_error = m_proxies.front().err;
      for (std::size_t i = 0; i < m_proxies.size(); ++i) {
        if (max_error < m_proxies[i].err)
          max_error = m_proxies[i].err;
      }
      if (min_error_change > max_error / FT(2.0))
        return false;
    }

    return true;
  }

  /*!
   * @brief Splits one proxy area via N-section (by default bisection).
   * @param px_idx proxy index.
   * @param n number of split sections.
   * @param nb_of_relaxations number of relaxation on the confined proxy area. (DEFINE CONFINED)
   * @return `true` if split succeeds, `false` otherwise.
   */
  bool split(const std::size_t px_idx,
    const std::size_t n = 2,
    const std::size_t nb_of_relaxations = 10) {
    if (px_idx >= m_proxies.size())
      return false;

    // collect confined proxy area
    std::list<face_descriptor> confined_area;
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
      if (get(m_fproxy_map, f) == px_idx)
        confined_area.push_back(f);
    // not enough faces to split
    if (n > confined_area.size())
      return false;

    // a copy of confined proxies
    std::vector<Proxy_wrapper> confined_proxies;
    confined_proxies.push_back(m_proxies[px_idx]);

    // select seed faces in the confined area
    std::size_t count = 1;
    BOOST_FOREACH(face_descriptor f, confined_area) {
      if (count >= n)
        break;

      if (f != m_proxies[px_idx].seed) {
        add_one_proxy_at(f);
        ++count;
        // copy
        confined_proxies.push_back(m_proxies.back());
      }
    }

    // relaxation on confined area and proxies
    for (std::size_t i = 0; i < nb_of_relaxations; ++i) {
      BOOST_FOREACH(face_descriptor f, confined_area)
        put(m_fproxy_map, f, CGAL_VSA_INVALID_TAG);

      partition(confined_proxies.begin(), confined_proxies.end());
      fit(confined_proxies.begin(), confined_proxies.end(), Concurrency_tag());
    }

    // copy back
    BOOST_FOREACH(const Proxy_wrapper &pxw, confined_proxies)
      m_proxies[pxw.idx] = pxw;

    return true;
  }
  /// @}

  /// \name Meshing
  /// @{
  /*!
   * @brief Extracts the output mesh in the form of an indexed triangle set.
   * @tparam NamedParameters a sequence of \ref vsa_namedparameters
   *
   * @param np optional sequence of \ref vsa_namedparameters among the ones listed below
   * @return `true` if the extracted surface mesh is manifold, `false` otherwise.
   *
   * \cgalNamedParamsBegin{Meshing Named Parameters}
   *   \cgalParamBegin{subdivision_ratio} chord subdivision ratio threshold to the chord length or average edge length.
   *   \cgalParamEnd
   *   \cgalParamBegin{relative_to_chord} set `true` if the subdivision_ratio is the ratio of the
   *     furthest vertex distance to the chord length, or to the average edge length otherwise.
   *   \cgalParamEnd
   *   \cgalParamBegin{with_dihedral_angle}  set `true` if subdivision_ratio is weighted by dihedral angle, `false` otherwise.
   *   \cgalParamEnd
   *   \cgalParamBegin{optimize_anchor_location}  set `true` if optimize the anchor locations, `false` otherwise.
   *   \cgalParamEnd
   *   \cgalParamBegin{pca_plane}  set `true` if use PCA plane fitting, otherwise use the default area averaged plane parameters.
   *   \cgalParamEnd
   * \cgalNamedParamsEnd
   */
  template <typename NamedParameters>
  bool extract_mesh(const NamedParameters &np) {
    using boost::get_param;
    using boost::choose_param;
    namespace vsa_np = CGAL::VSA::internal_np;

    const FT subdivision_ratio = choose_param(get_param(np, vsa_np::subdivision_ratio), FT(5.0));
    const bool relative_to_chord = choose_param(get_param(np, vsa_np::relative_to_chord), false);
    const bool with_dihedral_angle = choose_param(get_param(np, vsa_np::with_dihedral_angle), false);
    const bool optimize_anchor_location = choose_param(get_param(np, vsa_np::optimize_anchor_location), true);
    const bool pca_plane = choose_param(get_param(np, vsa_np::pca_plane), false);

    // compute averaged edge length, used in chord subdivision
    m_average_edge_length = compute_averaged_edge_length(*m_ptm, m_vpoint_map);

    // initialize all vertex anchor status
    BOOST_FOREACH(vertex_descriptor v, vertices(*m_ptm))
      put(m_vanchor_map, v, CGAL_VSA_INVALID_TAG);
    m_anchors.clear();
    m_bcycles.clear();
    m_tris.clear();
    m_px_planes.clear();

    // compute proxy planes, used for subdivision and anchor location
    compute_proxy_planes(pca_plane);

    // generate anchors
    find_anchors();
    find_edges(subdivision_ratio, relative_to_chord, with_dihedral_angle);
    add_anchors();

    // discrete constrained Delaunay triangulation
    pseudo_cdt();

    if (optimize_anchor_location)
      this->optimize_anchor_location();

    // check manifold-oriented
    CGAL::Polyhedron_3<Geom_traits> tm_test;
    return build_polyhedron_surface(tm_test);
  }
  /// @}

  /// \name Output
  /// @{
  /*!
   * @brief Outputs approximation results.
   * @tparam NamedParameters a sequence of \ref vsa_namedparameters

   * @param np optional sequence of \ref vsa_namedparameters among the ones listed below

   * \cgalNamedParamsBegin{Output Named Parameters}
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
  template <typename NamedParameters>
  void output(const NamedParameters &np) const {
    using boost::get_param;
    using boost::choose_param;
    namespace vsa_np = CGAL::VSA::internal_np;

    // get proxy map
    typedef typename boost::lookup_named_param_def<
      vsa_np::face_proxy_map_t,
      NamedParameters,
      vsa_np::dummy_output_t>::type Output_face_proxy_map;
    Output_face_proxy_map fproxymap = choose_param(
      get_param(np, vsa_np::face_proxy_map), vsa_np::dummy_output);
    vsa_np::face_proxy_map_helper(*this, fproxymap);

    // get proxies
    typedef typename boost::lookup_named_param_def<
      vsa_np::proxies_t,
      NamedParameters,
      vsa_np::dummy_output_t>::type Proxies_output_iterator;
    Proxies_output_iterator pxies_out_itr = choose_param(
      get_param(np, vsa_np::proxies), vsa_np::dummy_output);
    vsa_np::proxies_helper(*this, pxies_out_itr);

    // get anchor points
    typedef typename boost::lookup_named_param_def<
      vsa_np::anchors_t,
      NamedParameters,
      vsa_np::dummy_output_t>::type Anchors_output_iterator;
    Anchors_output_iterator apts_out_itr = choose_param(
      get_param(np, vsa_np::anchors) , vsa_np::dummy_output);
    vsa_np::anchors_helper(*this, apts_out_itr);

    // get indexed triangles
    typedef typename boost::lookup_named_param_def<
      vsa_np::triangles_t,
      NamedParameters,
      vsa_np::dummy_output_t>::type Triangles_output_iterator;
    Triangles_output_iterator tris_out_itr = choose_param(
      get_param(np, vsa_np::triangles) , vsa_np::dummy_output);
    vsa_np::triangles_helper(*this, tris_out_itr);
  }

  /// @cond CGAL_DOCUMENT_INTERNAL
  /*!
   * @brief Gets the face-proxy index map.
   * @tparam FaceProxyMap `WritablePropertyMap` with
   * `boost::graph_traits<TriangleMesh>::%face_descriptor` as key and `std::size_t` as value type
   * @param[out] face_proxy_map face proxy index map
   */
  template <typename FaceProxyMap>
  void proxy_map(FaceProxyMap &face_proxy_map) const {
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
      face_proxy_map[f] = get(m_fproxy_map, f);
  }

  /*!
   * @brief Gets the face region of the specified proxy.
   * @tparam OutputIterator output iterator with `boost::graph_traits<TriangleMesh>::%face_descriptor` as value type
   * @param px_idx proxy index
   * @param out output iterator
   */
  template <typename OutputIterator>
  void proxy_region(const std::size_t px_idx, OutputIterator out) const {
    if (px_idx >= m_proxies.size())
      return;

    BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
      if (get(m_fproxy_map, f) == px_idx)
        *out++ = f;
  }

  /*!
   * @brief Gets the proxies size.
   * @return number of proxies
   */
  std::size_t proxies_size() const { return m_proxies.size(); }

  /*!
   * @brief Gets the proxies.
   * @tparam OutputIterator output iterator with Proxy as value type
   * @param out output iterator
   */
  template <typename OutputIterator>
  void proxies(OutputIterator out) const {
    BOOST_FOREACH(const Proxy_wrapper &pxw, m_proxies)
      *out++ = pxw.px;
  }

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  /*!
   * @brief Gets the wrapped proxies.
   * @tparam OutputIterator output iterator with Proxy_wrapper as value type
   * @param out output iterator
   */
  template <typename OutputIterator>
  void wrapped_proxies(OutputIterator out) const {
    BOOST_FOREACH(const Proxy_wrapper &pxw, m_proxies)
      *out++ = pxw;
  }
#endif

  /*!
   * @brief Gets the anchor points, which have the area-averaged position of the projected anchor vertex points on the incident proxies.
   * @tparam OutputIterator output iterator with Point_3 as value type
   * @param out output iterator
   */
  template <typename OutputIterator>
  void anchor_points(OutputIterator out) const {
    BOOST_FOREACH(const Anchor &a, m_anchors)
      *out++ = a.pos;
  }

  /*!
   * @brief Gets the anchor vertices.
   * @tparam OutputIterator output iterator with vertex_descriptor as value type
   * @param out output iterator
   */
  template <typename OutputIterator>
  void anchor_vertices(OutputIterator out) const {
    BOOST_FOREACH(const Anchor &a, m_anchors)
      *out++ = a.vtx;
  }

  /*!
   * @brief Gets the indexed triangles, with
   * one triplet of integers per triangles, which refers to the anchor point indices.
   * @tparam OutputIterator output iterator with Indexed_triangle as value type
   * @param out output iterator
   */
  template <typename OutputIterator>
  void indexed_triangles(OutputIterator out) const {
    BOOST_FOREACH(const Indexed_triangle &t, m_tris)
      *out++ = t;
  }

  /*!
   * @brief Gets the indexed boundary polygon approximation.
   * @tparam OutputIterator output iterator with std::vector<std::size_t> as value type
   * @param out output iterator
   */
  template <typename OutputIterator>
  void indexed_boundary_polygons(OutputIterator out) const {
    BOOST_FOREACH(const Boundary_cycle &bcycle, m_bcycles) {
      std::vector<std::size_t> plg;
      halfedge_descriptor he = bcycle.he_head;
      do {
        Boundary_chord chord;
        walk_to_next_anchor(he, chord);
        plg.push_back(get(m_vanchor_map, target(he, *m_ptm)));
      } while (he != bcycle.he_head);
      *out++ = plg;
    }
  }
  /// @endcond
  /// @}

// private member functions
private:
  /*!
   * @brief Randomly initializes proxies to target number of proxies.
   * @note To ensure the randomness, call `std::srand()` beforehand.
   * @param max_nb_of_proxies maximum number of proxies, 
   * should be in range (nb_connected_components, num_faces(*m_ptm))
   * @param nb_of_iterations number of re-fitting iterations 
   * @return number of proxies initialized
   */
  std::size_t init_random(const std::size_t max_nb_of_proxies,
    const std::size_t nb_of_iterations) {
    // pick from current non seed faces randomly
    std::vector<face_descriptor> picked_seeds;
    if (random_pick_non_seed_faces(max_nb_of_proxies - m_proxies.size(), picked_seeds)) {
      BOOST_FOREACH(face_descriptor f, picked_seeds)
        add_one_proxy_at(f);
      run(nb_of_iterations);
    }

    return m_proxies.size();
  }

  /*!
   * @brief Incrementally initializes proxies to target number of proxies.
   * @param max_nb_of_proxies maximum number of proxies, 
   * should be in range (nb_connected_components, num_faces(*m_ptm))
   * @param nb_of_iterations number of re-fitting iterations 
   * before each incremental proxy insertion
   * @return number of proxies initialized
   */
  std::size_t init_incremental(const std::size_t max_nb_of_proxies,
    const std::size_t nb_of_iterations) {
    if (m_proxies.size() < max_nb_of_proxies)
      add_to_furthest_proxies(max_nb_of_proxies - m_proxies.size(), nb_of_iterations);

    return m_proxies.size();
  }

  /*!
   * @brief Hierarchically initializes proxies to target number of proxies.
   * @param max_nb_of_proxies maximum number of proxies, 
   * should be in range (nb_connected_components, num_faces(*m_ptm))
   * @param nb_of_iterations number of re-fitting iterations
   * before each hierarchical proxy insertion
   * @return number of proxies initialized
   */
  std::size_t init_hierarchical(const std::size_t max_nb_of_proxies,
    const std::size_t nb_of_iterations) {
    while (m_proxies.size() < max_nb_of_proxies) {
      // try to double current number of proxies each time
      std::size_t target_px = m_proxies.size();
      if (target_px * 2 > max_nb_of_proxies)
        target_px = max_nb_of_proxies;
      else
        target_px *= 2;
      // add proxies by error diffusion
      add_proxies_error_diffusion(target_px - m_proxies.size());
      run(nb_of_iterations);
    }

    return m_proxies.size();
  }

  /*!
   * @brief Randomly initializes proxies
   * with both maximum number of proxies and minimum error drop stop criteria,
   * where the first criterion met stops the seeding.
   * @note To ensure the randomness, call `std::srand()` beforehand.
   * @param max_nb_of_proxies maximum number of proxies, should be in range (nb_connected_components, num_faces(tm) / 3)
   * @param min_error_drop minimum error drop, should be in range (0.0, 1.0)
   * @param nb_of_iterations number of re-fitting iterations 
   * @return number of proxies initialized
   */
  std::size_t init_random_error(const std::size_t max_nb_of_proxies,
    const FT min_error_drop,
    const std::size_t nb_of_iterations) {

    const FT initial_err = compute_total_error();
    FT error_drop = min_error_drop * FT(2.0);
    while (m_proxies.size() < max_nb_of_proxies && error_drop > min_error_drop) {
      // try to double current number of proxies each time
      const std::size_t nb_px = m_proxies.size();
      const std::size_t nb_to_add =
        (nb_px * 2 > max_nb_of_proxies) ? max_nb_of_proxies - nb_px : nb_px;

      // pick from current non seed faces randomly
      std::vector<face_descriptor> picked_seeds;
      if (!random_pick_non_seed_faces(nb_to_add, picked_seeds))
        return m_proxies.size();

      BOOST_FOREACH(face_descriptor f, picked_seeds)
        add_one_proxy_at(f);
      const FT err = run(nb_of_iterations);
      error_drop = err / initial_err;
    }

    return m_proxies.size();
  }

  /*!
   * @brief Incrementally initializes proxies
   * with both maximum number of proxies and minimum error drop stop criteria,
   * The first criterion met stops the seeding.
   * @param max_nb_of_proxies maximum number of proxies, should be in range (nb_connected_components, num_faces(tm) / 3)
   * @param min_error_drop minimum error drop, should be in range (0.0, 1.0)
   * @param nb_of_iterations number of re-fitting iterations 
   * @return number of proxies initialized
   */
  std::size_t init_incremental_error(const std::size_t max_nb_of_proxies,
    const FT min_error_drop,
    const std::size_t nb_of_iterations) {
    const FT initial_err = compute_total_error();
    FT error_drop = min_error_drop * FT(2.0);
    while (m_proxies.size() < max_nb_of_proxies && error_drop > min_error_drop) {
      add_to_furthest_proxy();
      const FT err = run(nb_of_iterations);
      error_drop = err / initial_err;
    }

    return m_proxies.size();
  }

  /*!
   * @brief Hierarchically initializes proxies
   * with both maximum number of proxies and minimum error drop stop criteria,
   * where the first criterion met stops the seeding.
   * @param max_nb_of_proxies maximum number of proxies, should be in range (nb_connected_components, num_faces(tm) / 3)
   * @param min_error_drop minimum error drop, should be in range (0.0, 1.0)
   * @param nb_of_iterations number of re-fitting iterations 
   * @return number of proxies initialized
   */
  std::size_t init_hierarchical_error(const std::size_t max_nb_of_proxies,
    const FT min_error_drop,
    const std::size_t nb_of_iterations) {
    const FT initial_err = compute_total_error();
    FT error_drop = min_error_drop * FT(2.0);
    while (m_proxies.size() < max_nb_of_proxies && error_drop > min_error_drop) {
      // try to double current number of proxies each time
      std::size_t target_px = m_proxies.size();
      if (target_px * 2 > max_nb_of_proxies)
        target_px = max_nb_of_proxies;
      else
        target_px *= 2;
      add_proxies_error_diffusion(target_px - m_proxies.size());
      const FT err = run(nb_of_iterations);
      error_drop = err / initial_err;
    }

    return m_proxies.size();
  }

  /*!
   * @brief Partitions the area tagged with CGAL_VSA_INVALID_TAG with proxies, global face proxy map is updated.
   * Propagates the proxy seed faces and floods the tagged area to minimize the fitting error.
   * @tparam ProxyWrapperIterator forward iterator with Proxy_wrapper as value type
   * @param beg iterator point to the first element
   * @param end iterator point to the one past the last element
   */
  template<typename ProxyWrapperIterator>
  void partition(const ProxyWrapperIterator beg, const ProxyWrapperIterator end) {
    std::priority_queue<Face_to_integrate> face_pqueue;
    for (ProxyWrapperIterator pxw_itr = beg; pxw_itr != end; ++pxw_itr) {
      face_descriptor f = pxw_itr->seed;
      put(m_fproxy_map, f, pxw_itr->idx);

      BOOST_FOREACH(face_descriptor fadj, faces_around_face(halfedge(f, *m_ptm), *m_ptm)) {
        if (fadj != boost::graph_traits<TriangleMesh>::null_face()
            && get(m_fproxy_map, fadj) == CGAL_VSA_INVALID_TAG) {
          face_pqueue.push(Face_to_integrate(
            fadj, pxw_itr->idx, m_metric->compute_error(*m_ptm, fadj, pxw_itr->px)));
        }
      }
    }

    while (!face_pqueue.empty()) {
      const Face_to_integrate c = face_pqueue.top();
      face_pqueue.pop();
      if (get(m_fproxy_map, c.f) == CGAL_VSA_INVALID_TAG) {
        put(m_fproxy_map, c.f, c.px);
        BOOST_FOREACH(face_descriptor fadj, faces_around_face(halfedge(c.f, *m_ptm), *m_ptm)) {
          if (fadj != boost::graph_traits<TriangleMesh>::null_face()
            && get(m_fproxy_map, fadj) == CGAL_VSA_INVALID_TAG) {
            face_pqueue.push(Face_to_integrate(
              fadj, c.px, m_metric->compute_error(*m_ptm, fadj, m_proxies[c.px].px)));
          }
        }
      }
    }
  }

  /*!
   * @brief Refits and updates input range of proxies, sequential.
   * @tparam ProxyWrapperIterator forward iterator with Proxy_wrapper as value type
   * @param beg iterator point to the first element
   * @param end iterator point to the one past the last element
   * @param t concurrency tag
   */
  template<typename ProxyWrapperIterator>
  void fit(const ProxyWrapperIterator beg, const ProxyWrapperIterator end, const CGAL::Sequential_tag &) {
    std::vector<std::list<face_descriptor> > px_faces(m_proxies.size());
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
      px_faces[get(m_fproxy_map, f)].push_back(f);

    // update proxy parameters and seed
    for (ProxyWrapperIterator pxw_itr = beg; pxw_itr != end; ++pxw_itr) {
      const std::size_t px_idx = pxw_itr->idx;
      *pxw_itr = fit_proxy_from_patch(px_faces[px_idx], px_idx);
    }
  }

#ifdef CGAL_LINKED_WITH_TBB
  /*!
   * @brief Refits and updates input range of proxies, parallel.
   * @tparam ProxyWrapperIterator forward iterator with Proxy_wrapper as value type
   * @param beg iterator point to the first element
   * @param end iterator point to the one past the last element
   * @param t concurrency tag
   */
  template<typename ProxyWrapperIterator>
  void fit(const ProxyWrapperIterator beg, const ProxyWrapperIterator end, const CGAL::Parallel_tag &) {
    std::vector<std::list<face_descriptor> > px_faces(m_proxies.size());
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
      px_faces[get(m_fproxy_map, f)].push_back(f);

    // update proxy parameters and seed
    tbb::parallel_for(tbb::blocked_range<ProxyWrapperIterator>(beg, end),
      [&](tbb::blocked_range<ProxyWrapperIterator> &r) {
        for (ProxyWrapperIterator pxw_itr = r.begin(); pxw_itr != r.end(); ++pxw_itr) {
          const std::size_t px_idx = pxw_itr->idx;
          *pxw_itr = fit_proxy_from_patch(px_faces[px_idx], px_idx);
        }
      });
  }
#endif // CGAL_LINKED_WITH_TBB

  /*!
   * @brief Adds a proxy seed at the face with the maximum fitting error.
   * @return `true` add successfully, `false` otherwise
   */
  bool add_to_furthest_proxy() {
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
    std::cerr << "add furthest " << m_proxies.size() << std::endl;
#endif
    FT max_error = m_proxies.front().err;
    std::size_t px_worst = 0;
    for (std::size_t i = 0; i < m_proxies.size(); ++i) {
      if (max_error < m_proxies[i].err) {
        max_error = m_proxies[i].err;
        px_worst = i;
      }
    }

    face_descriptor fworst;
    bool first = true;
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm)) {
      std::size_t px_idx = get(m_fproxy_map, f);
      if (px_idx != px_worst || f == m_proxies[px_idx].seed)
        continue;

      FT err = m_metric->compute_error(*m_ptm, f, m_proxies[px_idx].px);
      if (first || max_error < err) {
        first = false;
        max_error = err;
        fworst = f;
      }
    }

    if (first)
      return false;

    add_one_proxy_at(fworst);

    return true;
  }

  /*!
   * @brief Fits a new (wrapped) proxy from a region patch.
   * 1. Compute proxy parameters from a list of faces.
   * 2. Find proxy seed face.
   * 3. Sum the proxy error.
   * @tparam FacePatch container with `face_descriptor` as data type
   * @param px_patch proxy patch container
   * @param px_idx proxy index
   * @return fitted wrapped proxy
   */
  template<typename FacePatch>
  Proxy_wrapper fit_proxy_from_patch(const FacePatch &px_patch, const std::size_t px_idx) {
    CGAL_assertion(!px_patch.empty());

    // use Proxy_fitting functor to fit proxy parameters
    const Proxy px = m_metric->fit_proxy(px_patch, *m_ptm);

    // find proxy seed and sum error
    face_descriptor seed = *px_patch.begin();
    FT err_min = m_metric->compute_error(*m_ptm, seed, px);
    FT sum_error(0.0);
    BOOST_FOREACH(face_descriptor f, px_patch) {
      const FT err = m_metric->compute_error(*m_ptm, f, px);
      sum_error += err;
      if (err < err_min) {
        err_min = err;
        seed = f;
      }
    }

    return Proxy_wrapper(px, px_idx, seed, sum_error);
  }

  /*!
   * @brief Adds a proxy at face f.
   * @param f where to the proxy is initialized from
   */
  void add_one_proxy_at(const face_descriptor f) {
    m_proxies.push_back(fit_proxy_from_face(f, m_proxies.size()));
  }

  /*!
   * @brief Fits a new (wrapped) proxy from a face.
   * 1. Compute proxy parameters from the face.
   * 2. Set seed to this face.
   * 3. Update the proxy error.
   * 4. Update proxy map.
   * @pre current face proxy map is valid
   * @param face_descriptor face
   * @param px_idx proxy index
   * @return fitted wrapped proxy
   */
  Proxy_wrapper fit_proxy_from_face(const face_descriptor f, const std::size_t px_idx) {
    // fit proxy parameters
    std::vector<face_descriptor> fvec(1, f);
    const Proxy px = m_metric->fit_proxy(fvec, *m_ptm);
    const FT err = m_metric->compute_error(*m_ptm, f, px);

    // original proxy map should always be falid
    const std::size_t prev_px_idx = get(m_fproxy_map, f);
    CGAL_assertion(prev_px_idx != CGAL_VSA_INVALID_TAG);
    // update the proxy error and proxy map
    m_proxies[prev_px_idx].err -= m_metric->compute_error(*m_ptm, f, m_proxies[prev_px_idx].px);
    put(m_fproxy_map, f, px_idx);

    return Proxy_wrapper(px, px_idx, f, err);
  }

  /*!
   * @brief Picks a number of non-seed faces into an empty vector randomly.
   * @param nb_requested requested number of faces
   * @param[out] picked_faces shuffled faces vector
   * @return `true` if requested number of faces are selected, `false` otherwise
   */
  bool random_pick_non_seed_faces(const std::size_t nb_requested,
    std::vector<face_descriptor> &picked_faces) {
    if (nb_requested + m_proxies.size() >= num_faces(*m_ptm))
      return false;

    std::set<face_descriptor> seed_faces_set;
    BOOST_FOREACH(const Proxy_wrapper &pxw, m_proxies)
      seed_faces_set.insert(pxw.seed);

    const std::size_t nb_nsf = num_faces(*m_ptm) - m_proxies.size();
    std::vector<face_descriptor> non_seed_faces;
    non_seed_faces.reserve(nb_nsf);
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm)) {
      if (seed_faces_set.find(f) != seed_faces_set.end())
        continue;
      non_seed_faces.push_back(f);
    }

    // random shuffle first few faces
    for (std::size_t i = 0; i < nb_requested; ++i) {
      // swap ith element with a random one
      std::size_t r = static_cast<std::size_t>(
        static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX) * 
        static_cast<double>(nb_nsf - 1));
      std::swap(non_seed_faces[i], non_seed_faces[r]);
    }

    for (std::size_t i = 0; i < nb_requested; ++i)
      picked_faces.push_back(non_seed_faces[i]);

    return true;
  }

  /*!
   * @brief Initializes proxies from each connected component of the input mesh.
   * @note This function clears proxy vector and sets face proxy map to initial state,
   * intended only for bootstrapping initialization.
   * Coarse approximation iteration is not performed, because it is inaccurate anyway
   * and may yield degenerate cases (e.g. a standard cube model).
   */
  void bootstrap_from_connected_components() {
    // set all faces invalid to mark as unvisited / untagged
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
      put(m_fproxy_map, f, CGAL_VSA_INVALID_TAG);

    // prepare for connected components visiting
    std::vector<std::list<face_descriptor> > cc_patches;
    bool if_all_visited = false;
    std::size_t cc_idx = 0;
    face_descriptor seed_face = *(faces(*m_ptm).first);
    while (!if_all_visited) {
      // use current seed face to traverse the conneceted componnets
      std::list<face_descriptor> cc_patch;
      cc_patch.push_back(seed_face);
      std::stack<face_descriptor> fstack;
      fstack.push(seed_face);
      put(m_fproxy_map, seed_face, cc_idx);
      while (!fstack.empty()) {
        face_descriptor active_face = fstack.top();
        fstack.pop();
        BOOST_FOREACH(face_descriptor fadj,
          faces_around_face(halfedge(active_face, *m_ptm), *m_ptm)) {
          if (fadj != boost::graph_traits<TriangleMesh>::null_face()
            && get(m_fproxy_map, fadj) == CGAL_VSA_INVALID_TAG) {
            cc_patch.push_back(fadj);
            fstack.push(fadj);
            put(m_fproxy_map, fadj, cc_idx);
          }
        }
      }
      cc_patches.push_back(cc_patch);
      // check if all visited
      if_all_visited = true;
      BOOST_FOREACH(face_descriptor f, faces(*m_ptm)) {
        if (get(m_fproxy_map, f) == CGAL_VSA_INVALID_TAG) {
          if_all_visited = false;
          ++cc_idx;
          seed_face = f;
          break;
        }
      }
    }

    m_proxies.clear();
    BOOST_FOREACH(const std::list<face_descriptor> &cc_patch, cc_patches)
      m_proxies.push_back(fit_proxy_from_patch(cc_patch, m_proxies.size()));
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
    std::cerr << "#cc " << m_proxies.size() << std::endl;
#endif
  }

  /*!
   * @brief Computes proxy planes.
   * The proxy may not contain the plane related properties, so we need these internal planes,
   * used in the chord subdivision and anchor location.
   * @param if_pca_plane set `true` to use the PCA plane fitting
   */
  void compute_proxy_planes(const bool if_pca_plane) {
    // fit proxy planes, areas, normals
    std::vector<std::list<face_descriptor> > px_faces(m_proxies.size());
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
      px_faces[get(m_fproxy_map, f)].push_back(f);

    BOOST_FOREACH(const std::list<face_descriptor> &px_patch, px_faces) {
      Plane_3 fit_plane = if_pca_plane ? 
        fit_plane_pca(px_patch.begin(), px_patch.end()) :
          fit_plane_area_averaged(px_patch.begin(), px_patch.end());

      Vector_3 norm = CGAL::NULL_VECTOR;
      FT area(0.0);
      BOOST_FOREACH(face_descriptor f, px_patch) {
        halfedge_descriptor he = halfedge(f, *m_ptm);
        const Point_3 &p0 = m_vpoint_map[source(he, *m_ptm)];
        const Point_3 &p1 = m_vpoint_map[target(he, *m_ptm)];
        const Point_3 &p2 = m_vpoint_map[target(next(he, *m_ptm), *m_ptm)];
        const FT farea(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
        const Vector_3 fnorm = CGAL::unit_normal(p0, p1, p2);

        norm = sum_functor(norm, scale_functor(fnorm, farea));
        area += farea;
      }
      norm = scale_functor(norm, FT(1.0 / std::sqrt(CGAL::to_double(norm.squared_length()))));

      m_px_planes.push_back(Proxy_plane(fit_plane, norm, area));
    }
  }

  /*!
   * @brief Finds the anchors.
   */
  void find_anchors() {
    BOOST_FOREACH(vertex_descriptor vtx, vertices(*m_ptm)) {
      std::size_t border_count = 0;

      BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(vtx, *m_ptm)) {
        if (CGAL::is_border_edge(h, *m_ptm))
          ++border_count;
        else if (get(m_fproxy_map, face(h, *m_ptm)) != get(m_fproxy_map, face(opposite(h, *m_ptm), *m_ptm)))
          ++border_count;
      }
      if (border_count >= 3)
        attach_anchor(vtx);
    }
  }

  /*!
   * @brief Finds and approximates the chord connecting the anchors.
   * @param subdivision_ratio boundary chord approximation recursive split creterion
   * @param relative_to_chord set `true` if the subdivision_ratio is relative to the chord length (relative sense),
   * otherwise it's relative to the average edge length (absolute sense).
   * @param with_dihedral_angle set `true` if add dihedral angle weight to the distance, `false` otherwise
   */
  void find_edges(const FT subdivision_ratio,
    const bool relative_to_chord,
    const bool with_dihedral_angle) {
    // collect candidate halfedges in a set
    std::set<halfedge_descriptor> he_candidates;
    BOOST_FOREACH(halfedge_descriptor h, halfedges(*m_ptm)) {
      if (!CGAL::is_border(h, *m_ptm)
        && (CGAL::is_border(opposite(h, *m_ptm), *m_ptm)
          || get(m_fproxy_map, face(h, *m_ptm)) != get(m_fproxy_map, face(opposite(h, *m_ptm), *m_ptm))))
        he_candidates.insert(h);
    }

    // pick up one candidate halfedge each time and traverse the connected boundary cycle
    while (!he_candidates.empty()) {
      halfedge_descriptor he_start = *he_candidates.begin();
      walk_to_first_anchor(he_start);
      // no anchor in this connected boundary cycle, make a new anchor
      if (!is_anchor_attached(he_start))
        attach_anchor(he_start);

      // a new connected boundary cycle
      m_bcycles.push_back(Boundary_cycle(he_start));

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
      std::cerr << "#bcycle " << m_bcycles.size() << std::endl;
#endif

      const halfedge_descriptor he_mark = he_start;
      do {
        Boundary_chord chord;
        walk_to_next_anchor(he_start, chord);
        m_bcycles.back().num_anchors += subdivide_chord(chord.begin(), chord.end(),
          subdivision_ratio, relative_to_chord, with_dihedral_angle);

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
        std::cerr << "#chord_anchor " << m_bcycles.back().num_anchors << std::endl;
#endif

        BOOST_FOREACH(const halfedge_descriptor &he, chord)
          he_candidates.erase(he);
      } while (he_start != he_mark);
    }
  }

  /*!
   * @brief Adds anchors to the boundary cycles with only 2 anchors.
   */
  void add_anchors() {
    BOOST_FOREACH(const Boundary_cycle &bcycle, m_bcycles) {
      if (bcycle.num_anchors > 2)
        continue;

      // 2 initial anchors at least
      CGAL_assertion(bcycle.num_anchors == 2);
      // borders with only 2 initial anchors
      Point_3 pt_begin = m_vpoint_map[target(bcycle.he_head, *m_ptm)];
      Point_3 pt_end = pt_begin;

      halfedge_descriptor he = bcycle.he_head;
      Boundary_chord chord;
      std::size_t count = 0;
      do {
        walk_to_next_border_halfedge(he);
        if (!is_anchor_attached(he))
          chord.push_back(he);
        else {
          if (count == 0)
            pt_end = m_vpoint_map[target(he, *m_ptm)];
          ++count;
        }
      } while (he != bcycle.he_head);

      // anchor count may be increased to more than 2 afterwards
      // due to the new anchors added by the neighboring boundary cycle (< 2 anchors)
      if (count > 2) {
        const_cast<Boundary_cycle &>(bcycle).num_anchors = count;
        continue;
      }

      FT dist_max(0.0);
      halfedge_descriptor he_max;
      Vector_3 chord_vec = vector_functor(pt_begin, pt_end);
      chord_vec = scale_functor(chord_vec,
        FT(1.0 / std::sqrt(CGAL::to_double(chord_vec.squared_length()))));
      BOOST_FOREACH(const halfedge_descriptor &he, chord) {
        Vector_3 vec = vector_functor(pt_begin, m_vpoint_map[target(he, *m_ptm)]);
        vec = CGAL::cross_product(chord_vec, vec);
        FT dist(std::sqrt(CGAL::to_double(vec.squared_length())));
        if (dist > dist_max) {
          dist_max = dist;
          he_max = he;
        }
      }
      // add one anchors to this boundary cycle
      attach_anchor(he_max);
      const_cast<Boundary_cycle &>(bcycle).num_anchors++;
    }
  }

  /*!
   * @brief Runs the pseudo Constrained Delaunay Triangulation at each proxy region,
   * and stores the extracted indexed triangles in @a tris.
   * @pre all anchors are found, i.e. all boundary cycles have been visited
   * and attached with at least 3 anchors.
   */
  void pseudo_cdt() {
    // subgraph attached with vertex anchor status and edge weight
    typedef boost::property<boost::vertex_index1_t, std::size_t,
      boost::property<boost::vertex_index2_t, std::size_t> > VertexProperty;
    typedef boost::property<boost::edge_weight_t, FT,
      boost::property<boost::edge_index_t, std::size_t> > EdgeProperty;
    typedef boost::subgraph<boost::adjacency_list<
      boost::listS, boost::vecS,
      boost::undirectedS,
      VertexProperty, EdgeProperty> > SubGraph;
    typedef typename boost::property_map<SubGraph, boost::vertex_index1_t>::type VertexIndex1Map;
    typedef typename boost::property_map<SubGraph, boost::vertex_index2_t>::type VertexIndex2Map;
    typedef typename boost::property_map<SubGraph, boost::edge_weight_t>::type EdgeWeightMap;
    typedef typename SubGraph::vertex_descriptor sg_vertex_descriptor;
    typedef std::vector<sg_vertex_descriptor> VertexVector;

    typedef boost::unordered_map<vertex_descriptor, sg_vertex_descriptor> VertexMap;
    typedef boost::associative_property_map<VertexMap> ToSGVertexMap;
    VertexMap vmap;
    ToSGVertexMap to_sgv_map(vmap);

    // mapping the TriangleMesh mesh into a SubGraph
    SubGraph gmain;
    VertexIndex1Map global_vanchor_map = get(boost::vertex_index1, gmain);
    VertexIndex2Map global_vtag_map = get(boost::vertex_index2, gmain);
    EdgeWeightMap global_eweight_map = get(boost::edge_weight, gmain);
    BOOST_FOREACH(vertex_descriptor v, vertices(*m_ptm)) {
      sg_vertex_descriptor sgv = add_vertex(gmain);
      global_vanchor_map[sgv] = get(m_vanchor_map, v);
      global_vtag_map[sgv] = get(m_vanchor_map, v);
      vmap.insert(std::pair<vertex_descriptor, sg_vertex_descriptor>(v, sgv));
    }
    BOOST_FOREACH(edge_descriptor e, edges(*m_ptm)) {
      vertex_descriptor vs = source(e, *m_ptm);
      vertex_descriptor vt = target(e, *m_ptm);
      FT len(std::sqrt(CGAL::to_double(
        CGAL::squared_distance(m_vpoint_map[vs], m_vpoint_map[vt]))));
      add_edge(to_sgv_map[vs], to_sgv_map[vt], len, gmain);
    }

    std::vector<VertexVector> vertex_patches(m_proxies.size());
    BOOST_FOREACH(vertex_descriptor v, vertices(*m_ptm)) {
      std::set<std::size_t> px_set;
      BOOST_FOREACH(face_descriptor f, faces_around_target(halfedge(v, *m_ptm), *m_ptm)) {
        if (f != boost::graph_traits<TriangleMesh>::null_face())
          px_set.insert(get(m_fproxy_map, f));
      }
      BOOST_FOREACH(std::size_t p, px_set)
        vertex_patches[p].push_back(to_sgv_map[v]);
    }
    BOOST_FOREACH(VertexVector &vpatch, vertex_patches) {
      // add a super vertex connecting to its boundary anchors in each patch
      const sg_vertex_descriptor superv = add_vertex(gmain);
      global_vanchor_map[superv] = CGAL_VSA_INVALID_TAG;
      global_vtag_map[superv] = CGAL_VSA_INVALID_TAG;
      BOOST_FOREACH(sg_vertex_descriptor v, vpatch) {
        if (is_anchor_attached(v, global_vanchor_map))
          add_edge(superv, v, FT(0.0), gmain);
      }
      vpatch.push_back(superv);
    }

    // multi-source Dijkstra's shortest path algorithm applied to each proxy patch
    BOOST_FOREACH(VertexVector &vpatch, vertex_patches) {
      // construct subgraph
      SubGraph &glocal = gmain.create_subgraph();
      BOOST_FOREACH(sg_vertex_descriptor v, vpatch)
        add_vertex(v, glocal);

      // most subgraph functions work with local descriptors
      VertexIndex1Map local_vanchor_map = get(boost::vertex_index1, glocal);
      VertexIndex2Map local_vtag_map = get(boost::vertex_index2, glocal);
      EdgeWeightMap local_eweight_map = get(boost::edge_weight, glocal);

      const sg_vertex_descriptor source = glocal.global_to_local(vpatch.back());
      VertexVector pred(num_vertices(glocal));
      boost::dijkstra_shortest_paths(glocal, source,
        boost::predecessor_map(&pred[0]).weight_map(local_eweight_map));

      // backtrack to the anchor and tag each vertex in the local patch graph
      BOOST_FOREACH(sg_vertex_descriptor v, vertices(glocal)) {
        // skip the added super source vertex in the patch
        if (v == source)
          continue;
        sg_vertex_descriptor curr = v;
        while (!is_anchor_attached(curr, local_vanchor_map))
          curr = pred[curr];
        local_vtag_map[v] = local_vanchor_map[curr];
      }
    }

    // tag all boundary chords
    BOOST_FOREACH(const Boundary_cycle &bcycle, m_bcycles) {
      halfedge_descriptor he = bcycle.he_head;
      do {
        Boundary_chord chord;
        walk_to_next_anchor(he, chord);

        std::vector<FT> vdist;
        vdist.push_back(FT(0.0));
        BOOST_FOREACH(halfedge_descriptor h, chord) {
          FT elen = global_eweight_map[edge(
            to_sgv_map[source(h, *m_ptm)],
            to_sgv_map[target(h, *m_ptm)],
            gmain).first];
          vdist.push_back(vdist.back() + elen);
        }

        FT half_chord_len = vdist.back() / FT(2.0);
        const std::size_t anchorleft = get(m_vanchor_map, source(chord.front(), *m_ptm));
        const std::size_t anchorright = get(m_vanchor_map, target(chord.back(), *m_ptm));
        typename std::vector<FT>::iterator ditr = vdist.begin() + 1;
        for (Boundary_chord_iterator citr = chord.begin(); citr != chord.end() - 1; ++citr, ++ditr) {
          if (*ditr < half_chord_len)
            global_vtag_map[to_sgv_map[target(*citr, *m_ptm)]] = anchorleft;
          else
            global_vtag_map[to_sgv_map[target(*citr, *m_ptm)]] = anchorright;
        }
      } while (he != bcycle.he_head);
    }

    // collect triangles
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm)) {
      halfedge_descriptor he = halfedge(f, *m_ptm);
      std::size_t i = global_vtag_map[to_sgv_map[source(he, *m_ptm)]];
      std::size_t j = global_vtag_map[to_sgv_map[target(he, *m_ptm)]];
      std::size_t k = global_vtag_map[to_sgv_map[target(next(he, *m_ptm), *m_ptm)]];
      if (i != j && i != k && j != k) {
        CGAL::cpp11::array<std::size_t, 3> t;
        t[0] = i;
        t[1] = j;
        t[2] = k;
        m_tris.push_back(t);
      }
    }
  }

  /*!
   * @brief Walks along the region boundary cycle to the first halfedge
   * pointing to a vertex associated with an anchor.
   * @param[in/out] he_start region boundary halfedge
   */
  void walk_to_first_anchor(halfedge_descriptor &he_start) {
    const halfedge_descriptor start_mark = he_start;
    while (!is_anchor_attached(he_start)) {
      // no anchor attached to the halfedge target
      walk_to_next_border_halfedge(he_start);
      if (he_start == start_mark) // back to where started, a boundary cycle without anchor
        return;
    }
  }

  /*!
   * @brief Walks along the region boundary cycle to the next anchor
   * and records the path as a `Boundary_chord`.
   * @param[in/out] he_start starting region boundary halfedge
   * pointing to a vertex associated with an anchor
   * @param[out] chord recorded path chord
   */
  void walk_to_next_anchor(halfedge_descriptor &he_start, Boundary_chord &chord) const {
    do {
      walk_to_next_border_halfedge(he_start);
      chord.push_back(he_start);
    } while (!is_anchor_attached(he_start));
  }

  /*!
   * @brief Walks to the next boundary cycle halfedge.
   * @param[in/out] he_start region boundary halfedge
   */
  void walk_to_next_border_halfedge(halfedge_descriptor &he_start) const {
    const std::size_t px_idx = get(m_fproxy_map, face(he_start, *m_ptm));
    BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(he_start, *m_ptm)) {
      if (CGAL::is_border(h, *m_ptm) || get(m_fproxy_map, face(h, *m_ptm)) != px_idx) {
        he_start = opposite(h, *m_ptm);
        return;
      }
    }
  }

  /*!
   * @brief Subdivides a chord recursively in range [@a chord_begin, @a chord_end).
   * @param chord_begin begin iterator of the chord
   * @param chord_end end iterator of the chord
   * @param subdivision_ratio the chord recursive split error threshold
   * @param relative_to_chord set `true` if the subdivision_ratio is relative to the the chord length (relative sense),
   * otherwise it's relative to the average edge length (absolute sense).
   * @param with_dihedral_angle set `true` if add dihedral angle weight to the distance, `false` otherwise
   * @return the number of anchors of the chord apart from the first one
   */
  std::size_t subdivide_chord(
    const Boundary_chord_iterator &chord_begin,
    const Boundary_chord_iterator &chord_end,
    const FT subdivision_ratio,
    const bool relative_to_chord,
    const bool with_dihedral_angle) {
    const std::size_t chord_size = std::distance(chord_begin, chord_end);
    const halfedge_descriptor he_first = *chord_begin;
    const halfedge_descriptor he_last = *(chord_end - 1);
    const std::size_t anchor_first = get(m_vanchor_map, source(he_first, *m_ptm));
    const std::size_t anchor_last = get(m_vanchor_map, target(he_last, *m_ptm));

    // do not subdivide trivial non-circular chord
    if ((anchor_first != anchor_last) && (chord_size < 4))
      return 1;

    bool if_subdivide = false;
    Boundary_chord_iterator chord_max;
    const Point_3 &pt_begin = m_vpoint_map[source(he_first, *m_ptm)];
    const Point_3 &pt_end = m_vpoint_map[target(he_last, *m_ptm)];
    if (anchor_first == anchor_last) {
      // circular chord
      CGAL_assertion(chord_size > 2);

      FT dist_max(0.0);
      for (Boundary_chord_iterator citr = chord_begin; citr != chord_end; ++citr) {
        FT dist = CGAL::squared_distance(pt_begin, m_vpoint_map[target(*citr, *m_ptm)]);
        dist = FT(std::sqrt(CGAL::to_double(dist)));
        if (dist > dist_max) {
          chord_max = citr;
          dist_max = dist;
        }
      }

      if_subdivide = true;
    }
    else {
      FT dist_max(0.0);
      Vector_3 chord_vec = vector_functor(pt_begin, pt_end);
      FT chord_len(std::sqrt(CGAL::to_double(chord_vec.squared_length())));
      chord_vec = scale_functor(chord_vec, FT(1.0) / chord_len);

      for (Boundary_chord_iterator citr = chord_begin; citr != chord_end; ++citr) {
        Vector_3 vec = vector_functor(pt_begin, m_vpoint_map[target(*citr, *m_ptm)]);
        vec = CGAL::cross_product(chord_vec, vec);
        FT dist(std::sqrt(CGAL::to_double(vec.squared_length())));
        if (dist > dist_max) {
          chord_max = citr;
          dist_max = dist;
        }
      }

      FT criterion = dist_max;
      if (relative_to_chord)
        criterion /= chord_len;
      else
        criterion /= m_average_edge_length;

      if (with_dihedral_angle) {
        // suppose the proxy normal angle is acute
        std::size_t px_left = get(m_fproxy_map, face(he_first, *m_ptm));
        std::size_t px_right = px_left;
        if (!CGAL::is_border(opposite(he_first, *m_ptm), *m_ptm))
          px_right = get(m_fproxy_map, face(opposite(he_first, *m_ptm), *m_ptm));
        FT norm_sin(1.0);
        if (!CGAL::is_border(opposite(he_first, *m_ptm), *m_ptm)) {
          Vector_3 vec = CGAL::cross_product(
            m_px_planes[px_left].normal, m_px_planes[px_right].normal);
          norm_sin = FT(std::sqrt(CGAL::to_double(scalar_product_functor(vec, vec))));
        }
        criterion *= norm_sin;
      }

      if (criterion > subdivision_ratio)
        if_subdivide = true;
    }

    if (if_subdivide) {
      // subdivide at the most remote vertex
      attach_anchor(*chord_max);

      const std::size_t num_left = subdivide_chord(chord_begin, chord_max + 1,
        subdivision_ratio, relative_to_chord, with_dihedral_angle);
      const std::size_t num_right = subdivide_chord(chord_max + 1, chord_end,
        subdivision_ratio, relative_to_chord, with_dihedral_angle);

      return num_left + num_right;
    }

    return 1;
  }

  /*!
   * @brief Tests if the target vertex of a halfedge is attached with an anchor.
   * @param he a halfedge descriptor
   * @return `true` is attached with an anchor, and `false` otherwise.
   */
  bool is_anchor_attached(const halfedge_descriptor &he) const {
    return is_anchor_attached(target(he, *m_ptm), m_vanchor_map);
  }

  /*!
   * @brief Checks if a vertex is attached with an anchor.
   * @tparam VertexAnchorIndexMap `WritablePropertyMap`
   * with `boost::graph_traights<TriangleMesh>::vertex_descriptor` as key and `std::size_t` as value type
   * @param vtx a vertex descriptor
   * @param vanchor_map vertex anchor index map
   */
  template<typename VertexAnchorIndexMap>
  bool is_anchor_attached(
    const typename boost::property_traits<VertexAnchorIndexMap>::key_type &vtx,
    const VertexAnchorIndexMap &vanchor_map) const {
    return get(vanchor_map, vtx) != CGAL_VSA_INVALID_TAG;
  }

  /*!
   * @brief Attaches an anchor to the vertex.
   * @param vtx vertex
   */
  void attach_anchor(const vertex_descriptor &vtx) {
    put(m_vanchor_map, vtx, m_anchors.size());
    // default anchor location is the vertex point
    m_anchors.push_back(Anchor(vtx, m_vpoint_map[vtx]));
  }

  /*!
   * @brief Attaches an anchor to the target vertex of the halfedge.
   * @param he halfedge
   */
  void attach_anchor(const halfedge_descriptor &he) {
    attach_anchor(target(he, *m_ptm));
  }

  /*!
   * @brief Optimizes the anchor location by averaging the projection points of
   * the anchor vertex to the incident proxy plane.
   */
  void optimize_anchor_location() {
    BOOST_FOREACH(Anchor &a, m_anchors) {
      const vertex_descriptor v = a.vtx;
      // incident proxy set
      std::set<std::size_t> px_set;
      BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(v, *m_ptm)) {
        if (!CGAL::is_border(h, *m_ptm))
          px_set.insert(get(m_fproxy_map, face(h, *m_ptm)));
      }

      // projection
      FT sum_area(0.0);
      Vector_3 vec = vector_functor(CGAL::ORIGIN, point_functor(CGAL::ORIGIN));
      const Point_3 vtx_pt = m_vpoint_map[v];
      BOOST_FOREACH(const std::size_t px_idx, px_set) {
        const Vector_3 proj = vector_functor(
          CGAL::ORIGIN, m_px_planes[px_idx].plane.projection(vtx_pt));
        const FT area = m_px_planes[px_idx].area;
        vec = sum_functor(vec, scale_functor(proj, area));
        sum_area += area;
      }
      vec = scale_functor(vec, FT(1.0) / sum_area);

      a.pos = translate_point_functor(CGAL::ORIGIN, vec);
    }
  }

  /*!
   * @brief Calculates the averaged edge length of a triangle mesh.
   * @param tm the input triangle mesh
   * @param vpoint_map vertex point map
   * @return averaged edge length
   */
  FT compute_averaged_edge_length(const TriangleMesh &tm, const VertexPointMap &vpoint_map) const {
    // compute average edge length
    FT sum(0.0);
    BOOST_FOREACH(edge_descriptor e, edges(tm)) {
      const vertex_descriptor vs = source(e, tm);
      const vertex_descriptor vt = target(e, tm);
      sum += FT(std::sqrt(CGAL::to_double(
        CGAL::squared_distance(vpoint_map[vs], vpoint_map[vt]))));
    }
    return sum / num_edges(tm);
  }

  /*!
   * @brief Uses an incremental builder to build and tests
   * if the indexed triangle surface is manifold.
   * @tparam PolyhedronSurface should be `CGAL::Polyhedron_3`
   * @param[out] poly input polyhedorn mesh
   * @return `true` if build manifold surface successfully, `false` otherwise
   */
  template <typename PolyhedronSurface>
  bool build_polyhedron_surface(PolyhedronSurface &poly) {
    std::vector<Point_3> vtx;
    BOOST_FOREACH(const Anchor &a, m_anchors)
      vtx.push_back(a.pos);

    typedef typename PolyhedronSurface::HalfedgeDS HDS;
    Triangle_polyhedron_builder<HDS> tpbuilder(vtx, m_tris);
    poly.delegate(tpbuilder);

    return tpbuilder.is_manifold;
  }

  /*!
   * @brief Fits an area averaged plane from a range of faces.
   * @tparam FaceIterator face_descriptor container iterator
   * @param beg container begin
   * @param end container end
   * @return fitted plane
   */
  template <typename FaceIterator>
  Plane_3 fit_plane_area_averaged(const FaceIterator &beg, const FaceIterator &end) {
    CGAL_assertion(beg != end);
    // area average normal and centroid
    Vector_3 norm = CGAL::NULL_VECTOR;
    Vector_3 cent = CGAL::NULL_VECTOR;
    FT sum_area(0.0);
    for (FaceIterator fitr = beg; fitr != end; ++fitr) {
      const halfedge_descriptor he = halfedge(*fitr, *m_ptm);
      const Point_3 &p0 = m_vpoint_map[source(he, *m_ptm)];
      const Point_3 &p1 = m_vpoint_map[target(he, *m_ptm)];
      const Point_3 &p2 = m_vpoint_map[target(next(he, *m_ptm), *m_ptm)];

      Vector_3 vec = vector_functor(CGAL::ORIGIN, CGAL::centroid(p0, p1, p2));
      FT farea(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
      Vector_3 fnorm = CGAL::unit_normal(p0, p1, p2);

      norm = sum_functor(norm, scale_functor(fnorm, farea));
      cent = sum_functor(cent, scale_functor(vec, farea));
      sum_area += farea;
    }
    norm = scale_functor(norm,
      FT(1.0 / std::sqrt(CGAL::to_double(norm.squared_length()))));
    cent = scale_functor(cent, FT(1.0) / sum_area);

    return Plane_3(CGAL::ORIGIN + cent, norm);
  }

  /*!
   * @brief Fits a plane from a range of faces with PCA algorithm.
   * @tparam FaceIterator face_descriptor container iterator
   * @param beg container begin
   * @param end container end
   * @return fitted plane
   */
  template <typename FaceIterator>
  Plane_3 fit_plane_pca(const FaceIterator &beg, const FaceIterator &end) {
    CGAL_assertion(beg != end);

    typedef typename Geom_traits::Triangle_3 Triangle_3;
    std::list<Triangle_3> tri_list;
    for (FaceIterator fitr = beg; fitr != end; ++fitr) {
      halfedge_descriptor he = halfedge(*fitr, *m_ptm);
      const Point_3 &p0 = m_vpoint_map[source(he, *m_ptm)];
      const Point_3 &p1 = m_vpoint_map[target(he, *m_ptm)];
      const Point_3 &p2 = m_vpoint_map[target(next(he, *m_ptm), *m_ptm)];
      tri_list.push_back(Triangle_3(p0, p1, p2));
    }

    // construct and fit proxy plane
    Plane_3 fit_plane;
    CGAL::linear_least_squares_fitting_3(
      tri_list.begin(),
      tri_list.end(),
      fit_plane,
      CGAL::Dimension_tag<2>());
    
    return fit_plane;
  }
};

} // end namespace CGAL

#undef CGAL_VSA_INVALID_TAG

#endif // CGAL_VARIATIONAL_SHAPE_APPROXIMATION_H
