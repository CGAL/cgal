#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_APPROXIMATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_APPROXIMATION_H

#include <CGAL/license/Surface_mesh_approximation.h>


#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <CGAL/vsa_metrics.h>
#include <CGAL/Default.h>

#include <boost/unordered_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>

#include <vector>
#include <queue>
#include <iterator>
#include <cmath>
#include <cstdlib>

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
#include <iostream>
#endif

namespace CGAL {
namespace VSA {

/*!
 * \ingroup PkgTSMA
 * @brief Seeding method enumeration for Variational Shape Approximation algorithm.
 */
  enum Seeding {
    /// Random seeding
    Random,
    /// Incremental seeding
    Incremental,
    /// Hierarchical seeding
    Hierarchical
  };

/*!
 * \ingroup PkgTSMA
 * @brief Main class for Variational Shape Approximation algorithm.
 * @tparam TriangleMesh a CGAL TriangleMesh
 * @tparam VertexPointMap vertex point map
 * @tparam ErrorMetric error metric type
 * @tparam ProxyFitting proxy fitting type
 * @tparam GeomTraits geometric traits type
 */
template <typename TriangleMesh,
  typename VertexPointMap,
  typename ErrorMetric = CGAL::Default,
  typename ProxyFitting = CGAL::Default,
  typename GeomTraits = CGAL::Default>
class Mesh_approximation {
// public typedefs
public:
  // Default typedefs
  /// Geometric trait type
  typedef typename CGAL::Default::Get<
    GeomTraits,
    typename Kernel_traits<
      typename boost::property_traits<VertexPointMap>::value_type
    >::Kernel >::type Geom_traits;
  /// ErrorMetric type
  typedef typename CGAL::Default::Get<ErrorMetric,
    CGAL::VSA::L21_metric<TriangleMesh, VertexPointMap, false, Geom_traits> >::type Error_metric;
  /// ProxyFitting type
  typedef typename CGAL::Default::Get<ProxyFitting,
    CGAL::VSA::L21_proxy_fitting<TriangleMesh, VertexPointMap, Geom_traits> >::type Proxy_fitting;
  /// Proxy type
  typedef typename Error_metric::Proxy Proxy;

// private typedefs and data member
private:
  // Geom_traits typedefs
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point_3;
  typedef typename Geom_traits::Vector_3 Vector_3;
  typedef typename Geom_traits::Plane_3 Plane_3;
  typedef typename Geom_traits::Construct_vector_3 Construct_vector_3;
  typedef typename Geom_traits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename Geom_traits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename Geom_traits::Compute_scalar_product_3 Compute_scalar_product_3;

  // graph_traits typedefs
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  // internal typedefs
  typedef boost::associative_property_map<boost::unordered_map<vertex_descriptor, int> > Vertex_anchor_map;

  typedef std::vector<halfedge_descriptor> Chord_vector;
  typedef typename Chord_vector::iterator Chord_vector_iterator;

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
public:
#endif
  // The proxy wrapper for approximation.
  struct Proxy_wrapper {
    Proxy_wrapper(const Proxy &_p, const face_descriptor &_s)
      : px(_p), seed(_s), err(0.0) {}

    Proxy px; // parameterized proxy
    face_descriptor seed; // proxy seed
    FT err; // proxy fitting error
  };
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
private:
#endif

  // The proxy fitting plane for meshing.
  struct Proxy_plane {
    Proxy_plane(const Plane_3 &_p, const Vector_3 &_n, const FT &_a)
      : plane(_p), normal(_n), area(_a) {}

    Plane_3 plane;
    Vector_3 normal;
    FT area;
  };

  // The facet candidate to be queued.
  struct Facet_to_integrate {
    Facet_to_integrate(const face_descriptor &_f, const std::size_t &_px, const FT &_err)
      : f(_f), px(_px), err(_err) {}

    bool operator<(const Facet_to_integrate &rhs) const {
      return err > rhs.err;
    }

    face_descriptor f; // facet
    std::size_t px; // proxy index
    FT err; // fitting error
  };

  // Proxy error with its index.
  struct Proxy_error {
    Proxy_error(const std::size_t &_px, const FT &_err)
      : px(_px), err(_err) {}

    // in ascending order
    bool operator<(const Proxy_error &rhs) const {
      return err < rhs.err;
    }

    std::size_t px;
    FT err;
  };

  // The average positioned anchor attached to a vertex.
  struct Anchor {
    Anchor(const vertex_descriptor &_vtx, const Point_3 _pos)
      : vtx(_vtx), pos(_pos) {}

    vertex_descriptor vtx; // The associated vertex.
    Point_3 pos; // The position of the anchor.
  };

  // The border cycle of a region.
  // One region may have multiple border cycles.
  struct Border {
    Border(const halfedge_descriptor &h)
      : he_head(h), num_anchors(0) {}

    halfedge_descriptor he_head; // The heading halfedge of the border cylce.
    std::size_t num_anchors; // The number of anchors on the border.
  };

  // Triangle polyhedron builder.
  template <typename HDS>
  class Triangle_polyhedron_builder : public CGAL::Modifier_base<HDS> {
    const std::vector<Point_3> &vtxs;
    const std::vector<std::vector<std::size_t> > &tris;
  public:
    bool is_manifold;
    Triangle_polyhedron_builder(const std::vector<Point_3> &_vtxs,
      const std::vector<std::vector<std::size_t> > &_tris)
      : vtxs(_vtxs), tris(_tris), is_manifold(true) {}

    void operator()(HDS &hds) {
      CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);
      typedef typename HDS::Vertex Vertex;
      typedef typename Vertex::Point Point;
      builder.begin_surface(vtxs.size(), tris.size());
      BOOST_FOREACH(const Point_3 &v, vtxs)
        builder.add_vertex(Point(v));
      BOOST_FOREACH(const std::vector<std::size_t> &t, tris) {
        std::vector<std::size_t>::const_iterator itr = t.begin();
        if (builder.test_facet(itr, itr + 3)) {
          builder.begin_facet();
          builder.add_vertex_to_facet(*itr);
          builder.add_vertex_to_facet(*(itr + 1));
          builder.add_vertex_to_facet(*(itr + 2));
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
  const TriangleMesh *m_pmesh;
  // The mesh vertex point map.
  VertexPointMap point_pmap;
  // The error metric.
  const Error_metric *fit_error;
  // The proxy fitting functor.
  const Proxy_fitting *proxy_fitting;

  Construct_vector_3 vector_functor;
  Construct_scaled_vector_3 scale_functor;
  Construct_sum_of_vectors_3 sum_functor;
  Compute_scalar_product_3 scalar_product_functor;

  // The facet proxy index map.
  boost::unordered_map<face_descriptor, std::size_t> internal_fidx_map;
  boost::associative_property_map<boost::unordered_map<face_descriptor, std::size_t> > fproxy_map;
  // The attached anchor index of a vertex.
  boost::unordered_map<vertex_descriptor, int> internal_vidx_map;
  Vertex_anchor_map vanchor_map;

  // Proxies.
  std::vector<Proxy_wrapper> proxies;
  // Proxy planes
  std::vector<Proxy_plane> px_planes;

  // All anchors.
  std::vector<Anchor> anchors;
  // All borders cycles.
  std::vector<Border> borders;
  // The indexed triangle approximation.
  std::vector<std::vector<std::size_t> > tris;

//member functions
public:
  /*!
   * %Default constructor.
   */
  Mesh_approximation() :
    m_pmesh(NULL),
    fit_error(NULL),
    proxy_fitting(NULL),
    fproxy_map(internal_fidx_map),
    vanchor_map(internal_vidx_map) {
    Geom_traits traits;
    vector_functor = traits.construct_vector_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scalar_product_functor = traits.compute_scalar_product_3_object();
  }

  /*!
   * Initialize and prepare for the approximation.
   * @param _mesh `CGAL TriangleMesh` on which approximation operates.
   * @param _point_pmap vertex point map of the mesh
   */
  Mesh_approximation(const TriangleMesh &_mesh,
    const VertexPointMap &_point_pmap) :
    m_pmesh(&_mesh),
    point_pmap(_point_pmap),
    fit_error(NULL),
    proxy_fitting(NULL),
    fproxy_map(internal_fidx_map),
    vanchor_map(internal_vidx_map) {
    Geom_traits traits;
    vector_functor = traits.construct_vector_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scalar_product_functor = traits.compute_scalar_product_3_object();
  }

  /*!
   * Set the mesh for approximation and rebuild the internal data structure.
   * @pre @a _mesh.is_pure_triangle()
   * @param _mesh `CGAL TriangleMesh` on which approximation operates.
   * @param _point_pmap vertex point map of the mesh
   */
  void set_mesh(const TriangleMesh &_mesh, const VertexPointMap &_point_pmap) {
    m_pmesh = &_mesh;
    point_pmap = _point_pmap;
    rebuild();
  }

  /*!
   * Set the error and fitting functor.
   * @param _error_metric a `ErrorMetric` functor.
   * @param _proxy_fitting a `ProxyFitting` functor.
   */
  void set_metric(const Error_metric &_error_metric,
    const Proxy_fitting &_proxy_fitting) {
    fit_error = &_error_metric;
    proxy_fitting = &_proxy_fitting;
  }

  /*!
   * Rebuild the internal data structure.
   */
  void rebuild() {

    // cleanup
    proxies.clear();
    px_planes.clear();
    anchors.clear();
    borders.clear();
    tris.clear();

    if (!m_pmesh)
      return;

    // rebuild internal data structure
    internal_fidx_map.clear();
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      internal_fidx_map[f] = 0;

    internal_vidx_map.clear();
    BOOST_FOREACH(vertex_descriptor v, vertices(*m_pmesh))
      internal_vidx_map.insert(std::pair<vertex_descriptor, int>(v, 0));
  }

  /*!
   * @brief Initialize the seeds with both maximum number of proxies
   * and minmum error drop stop criteria.
   * The first criterion met stops the seeding.
   * Parameters out of range are ignored.
   * @param method seeding method
   * @param max_nb_proxies maximum target number of proxies,
   * should be in range (0, num_faces(_mesh) / 3)
   * @param min_error_drop minimum error drop,
   * should be in range (0.0, 1.0)
   * @param nb_relaxations number of interleaved refitting relaxations
   * in incremental and hierarchical seeding
   * @return number of proxies initialized
   */
  std::size_t seeding(const Seeding method = Hierarchical,
    const boost::optional<std::size_t> max_nb_proxies = boost::optional<std::size_t>(),
    const boost::optional<FT> min_error_drop = boost::optional<FT>(),
    const std::size_t nb_relaxations = 5) {
    // maximum number of proxies internally, maybe better choice?
    const std::size_t nb_px = num_faces(*m_pmesh) / 3;

    if (min_error_drop && *min_error_drop > FT(0.0) && *min_error_drop < FT(1.0)) {
      // as long as minimum error is specified and valid
      // maximum number of proxies always exist, no matter specified or not or out of range
      // there is always a maximum number of proxies explicitly (max_nb_proxies) or implicitly (nb_px)
      std::size_t max_nb_px_adjusted = nb_px;
      if (max_nb_proxies && *max_nb_proxies < nb_px && *max_nb_proxies > 0)
        max_nb_px_adjusted = *max_nb_proxies;
      switch (method) {
        case Random:
          return init_random_error(max_nb_px_adjusted, *min_error_drop, nb_relaxations);
        case Incremental:
          return init_incremental_error(max_nb_px_adjusted, *min_error_drop, nb_relaxations);
        case Hierarchical:
          return init_hierarchical_error(max_nb_px_adjusted, *min_error_drop, nb_relaxations);
        default:
          return 0;
      }
    }
    else if (max_nb_proxies && *max_nb_proxies < nb_px && *max_nb_proxies > 0) {
      // no valid min_error_drop provided, only max_nb_proxies
      switch (method) {
        case Random:
          return init_random(*max_nb_proxies, nb_relaxations);
        case Incremental:
          return init_incremental(*max_nb_proxies, nb_relaxations);
        case Hierarchical:
          return init_hierarchical(*max_nb_proxies, nb_relaxations);
        default:
          return 0;
      }
    }
    else {
      // both parameters are unspecified or out of range
      const FT e(0.1);
      switch (method) {
        case Random:
          return init_random_error(nb_px, e, nb_relaxations);
        case Incremental:
          return init_incremental_error(nb_px, e, nb_relaxations);
        case Hierarchical:
          return init_hierarchical_error(nb_px, e, nb_relaxations);
        default:
          return 0;
      }
    }
  }

  /*!
   * @brief Run the partitioning and fitting processes.
   * @param nb_iterations number of iterations.
   */
  void run(std::size_t nb_iterations = 1) {
    for (std::size_t i = 0; i < nb_iterations; ++i) {
      partition();
      fit();
    }
  }

  /*!
   * @brief Run the partitioning and fitting process until no significant error change
   * @param cvg_threshold the percentage of error change between two successive runs,
   * should be in range (0, 1).
   * @param max_iterations maximum number of iterations allowed
   * @param avg_interval size of error average interval to have smoother convergence curve,
   * if 0 is assigned, 1 is used instead.
   * @return true if converged before hitting the maximum iterations, false otherwise
   */
  bool run_to_converge(const FT cvg_threshold,
    const std::size_t max_iterations = 100,
    std::size_t avg_interval = 3) {
    if (avg_interval == 0)
      avg_interval = 1;
    FT drop_pct(0.0);
    FT pre_err = compute_fitting_error();
    for (std::size_t itr_count = 0; itr_count < max_iterations; itr_count += avg_interval) {
      if (pre_err == FT(0.0))
        return true;

      FT avg_err(0.0);
      for (std::size_t i = 0; i < avg_interval; ++i) {
        run();
        avg_err += compute_fitting_error();
      }
      avg_err /= static_cast<FT>(avg_interval);

      drop_pct = (pre_err - avg_err) / pre_err;
      // the error may fluctuates
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
  FT compute_fitting_error() {
    BOOST_FOREACH(Proxy_wrapper &pxw, proxies)
      pxw.err = FT(0.0);
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
      std::size_t pxidx = fproxy_map[f];
      proxies[pxidx].err += (*fit_error)(f, proxies[pxidx].px);
    }

    FT sum_error(0.0);
    BOOST_FOREACH(const Proxy_wrapper &pxw, proxies)
      sum_error += pxw.err;

    return sum_error;
  }

  /*!
   * @brief Partition the geometry with current proxies.
   * Propagates the proxy seed facets and floods the whole mesh to minimize the fitting error.
   */
  void partition() {
#define CGAL_NOT_TAGGED_ID std::numeric_limits<std::size_t>::max()
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      fproxy_map[f] = CGAL_NOT_TAGGED_ID;

    std::priority_queue<Facet_to_integrate> facet_pqueue;
    for (std::size_t i = 0; i < proxies.size(); ++i) {
      face_descriptor f = proxies[i].seed;
      fproxy_map[f] = i;

      BOOST_FOREACH(face_descriptor fadj, faces_around_face(halfedge(f, *m_pmesh), *m_pmesh)) {
        if (fadj != boost::graph_traits<TriangleMesh>::null_face()
            && fproxy_map[fadj] == CGAL_NOT_TAGGED_ID) {
          facet_pqueue.push(Facet_to_integrate(
            fadj, i, (*fit_error)(fadj, proxies[i].px)));
        }
      }
    }

    while (!facet_pqueue.empty()) {
      const Facet_to_integrate c = facet_pqueue.top();
      facet_pqueue.pop();
      if (fproxy_map[c.f] == CGAL_NOT_TAGGED_ID) {
        fproxy_map[c.f] = c.px;
        BOOST_FOREACH(face_descriptor fadj, faces_around_face(halfedge(c.f, *m_pmesh), *m_pmesh)) {
          if (fadj != boost::graph_traits<TriangleMesh>::null_face()
            && fproxy_map[fadj] == CGAL_NOT_TAGGED_ID) {
            facet_pqueue.push(Facet_to_integrate(
            fadj, c.px, (*fit_error)(fadj, proxies[c.px].px)));
          }
        }
      }
    }
#undef CGAL_NOT_TAGGED_ID
  }

  /*!
   * @brief Refitting of current partitioning, update proxy parameters.
   * Calculates and updates the fitting proxies of current partition.
   */
  void fit() {
    std::vector<std::list<face_descriptor> > px_facets(proxies.size());
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      px_facets[fproxy_map[f]].push_back(f);

    // update proxy parameters and seed
    for (std::size_t i = 0; i < proxies.size(); ++i)
      proxies[i] = fit_new_proxy(px_facets[i].begin(), px_facets[i].end());
  }

  /*!
   * @brief Add proxies to the worst regions one by one.
   * The re-fitting is performed after each proxy is inserted.
   * @pre current facet proxy map is valid
   * @note after the addition, the facet proxy map remains valid
   * @param num_proxies number of proxies to be added
   * @param num_iterations number of re-fitting iterations
   * @return number of proxies added
   */
  std::size_t add_proxies_furthest(const std::size_t num_proxies,
    const std::size_t num_iterations = 5) {
    std::size_t num_added = 0;
    for (; num_added < num_proxies; ++num_added) {
      if (!add_proxy_furthest())
        break;
      for (std::size_t i = 0; i < num_iterations; ++i) {
        partition();
        fit();
      }
    }
    return num_added;
  }

  /*!
   * @brief Add proxies by diffusing fitting error into current partitions.
   * Each partition is added with the number of proxies in proportional to its fitting error.
   * @pre current facet proxy map is valid
   * @note after the addition, the facet proxy map is invalid
   * @param num_proxies number of proxies to be added
   * @return number of proxies successfully added
   */
  std::size_t add_proxies_error_diffusion(const std::size_t num_proxies) {
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
    std::cerr << "#px " << proxies.size() << std::endl;
#endif

    const FT sum_error = compute_fitting_error();
    const FT avg_error = sum_error / FT(static_cast<double>(num_proxies));

    std::vector<Proxy_error> px_error;
    for (std::size_t i = 0; i < proxies.size(); ++i)
      px_error.push_back(Proxy_error(i, proxies[i].err));
    // sort partition by error
    std::sort(px_error.begin(), px_error.end());

    // number of proxies to be added to each region
    std::vector<std::size_t> num_to_add(proxies.size(), 0);
    if (avg_error == FT(0.0)) {
      // rare case on extremely regular geometry like a cube
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
      std::cerr << "zero error, diffuse w.r.t. number of facets" << std::endl;
#endif
      const FT avg_facet = FT(
        static_cast<double>(num_faces(*m_pmesh)) / static_cast<double>(num_proxies));
      std::vector<FT> px_size(proxies.size(), FT(0.0));
      BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
        px_size[fproxy_map[f]] += FT(1.0);
      FT residual(0.0);
      for (std::size_t i = 0; i < proxies.size(); ++i) {
        FT to_add = (residual + px_size[i]) / avg_facet;
        FT floor_to_add = FT(std::floor(CGAL::to_double(to_add)));
        FT q_to_add = FT(CGAL::to_double(
          ((to_add - floor_to_add) > FT(0.5)) ? (floor_to_add + FT(1.0)) : floor_to_add));
        residual = (to_add - q_to_add) * avg_facet;
        num_to_add[i] = static_cast<std::size_t>(CGAL::to_double(q_to_add));
      }
    }
    else {
      // residual from previous proxy in range (-0.5, 0.5] * avg_error
      FT residual(0.0);
      BOOST_FOREACH(const Proxy_error &pxe, px_error) {
        // add error residual from previous proxy
        // to_add maybe negative but greater than -0.5
        FT to_add = (residual + pxe.err) / avg_error;
        // floor_to_add maybe negative but no less than -1
        FT floor_to_add = FT(std::floor(CGAL::to_double(to_add)));
        FT q_to_add = FT(CGAL::to_double(
          ((to_add - floor_to_add) > FT(0.5)) ? (floor_to_add + FT(1.0)) : floor_to_add));
        residual = (to_add - q_to_add) * avg_error;
        num_to_add[pxe.px] = static_cast<std::size_t>(CGAL::to_double(q_to_add));
      }
    }

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
    for (std::size_t i = 0; i < px_error.size(); ++i)
      std::cerr << "#px " << px_error[i].px
        << ", #error " << px_error[i].err
        << ", #num_to_add " << num_to_add[px_error[i].px] << std::endl;
#endif

    std::size_t num_added = 0;
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
      const std::size_t px_id = fproxy_map[f];
      if (proxies[px_id].seed == f)
        continue;

      if (num_to_add[px_id] > 0) {
        proxies.push_back(fit_new_proxy(f));
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

  /*!
   * @brief Teleport the local minima to the worst region, this combines the merging and adding processes.
   * The re-fitting is performed after each teleportation.
   * Here if we specify more than one proxy this means we teleport in a naive iterative fashion.
   * @param num_proxies number of proxies request to teleport
   * @param num_iterations number of re-fitting iterations
   * @param if_test true if do the merge test before the teleportation (attempt to escape from local minima).
   * @return number of proxies teleported.
   */
  std::size_t teleport_proxies(const std::size_t num_proxies,
    const std::size_t num_iterations = 5,
    const bool if_test = true) {
    std::size_t num_teleported = 0;
    while (num_teleported < num_proxies) {
      // find worst proxy
      compute_fitting_error();
      std::size_t px_worst = 0;
      FT max_error = proxies.front().err;
      for (std::size_t i = 0; i < proxies.size(); ++i) {
        if (max_error < proxies[i].err) {
          max_error = proxies[i].err;
          px_worst = i;
        }
      }
      bool found = false;
      face_descriptor tele_to;
      BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
        if (fproxy_map[f] == px_worst && f != proxies[px_worst].seed) {
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
      if (!find_best_merge(px_enlarged, px_merged, if_test))
        return num_teleported;
      if (px_worst == px_enlarged || px_worst == px_merged)
        return num_teleported;

      // teleport to a facet of the worst region
      // update merged proxies
      std::list<face_descriptor> merged_patch;
      BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
        std::size_t &px_idx = fproxy_map[f];
        if (px_idx == px_enlarged || px_idx == px_merged) {
          px_idx = px_enlarged;
          merged_patch.push_back(f);
        }
      }
      proxies[px_enlarged] = fit_new_proxy(merged_patch.begin(), merged_patch.end());
      // replace the merged proxy position to the newly teleported proxy
      proxies[px_merged] = fit_new_proxy(tele_to);
      fproxy_map[tele_to] = px_merged;

      num_teleported++;
      // coarse re-fitting
      for (std::size_t i = 0; i < num_iterations; ++i) {
        partition();
        fit();
      }

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
      std::cerr << "teleported" << std::endl;
#endif
    }

    return num_teleported;
  }

  /*!
   * @brief Merge two specified adjacent regions.
   * The overall re-fitting is not performed and the proxy map is maintained.
   * @pre two proxies must be adjacent
   * @param px0 the enlarged proxy
   * @param px1 the merged proxy
   * @return change of error
   */
  FT merge(std::size_t px0, std::size_t px1) {
    if (px0 >= proxies.size() || px1 >= proxies.size() || px0 == px1)
      return FT(0.0);

    // ensure px0 < px1
    if (px0 > px1)
      std::swap(px0, px1);

    // merge px1 to px0
    FT err_sum(0.0);
    std::list<face_descriptor> merged_patch;
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
      std::size_t px_idx = fproxy_map[f];
      if (px_idx == px1 || px_idx == px0) {
        err_sum += (*fit_error)(f, proxies[px_idx].px);
        fproxy_map[f] = px0;
        merged_patch.push_back(f);
      }
    }
    proxies[px0] = fit_new_proxy(merged_patch.begin(), merged_patch.end());

    proxies.erase(proxies.begin() + px1);
    // update facet proxy map
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
      if (fproxy_map[f] > px1)
        --fproxy_map[f];
    }

    FT err_merged(0.0);
    BOOST_FOREACH(face_descriptor f, merged_patch)
      err_merged += (*fit_error)(f, proxies[px0].px);

    return err_merged - err_sum;
  }

  /*!
   * @brief Find the best two regions to merge.
   * TODO: define 'best', it is minimum merged sum error now
   * @param px_enlarged the proxy to be enlarged
   * @param px_merged the proxy to be merged
   * @param if_test set true to activate the merge test
   * @return true if found, false otherwise
   */
  bool find_best_merge(std::size_t &px_enlarged, std::size_t &px_merged, const bool if_test) {
    typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
    typedef std::pair<std::size_t, std::size_t> ProxyPair;
    typedef std::set<ProxyPair> MergedPair;

    std::vector<std::list<face_descriptor> > px_facets(proxies.size());
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      px_facets[fproxy_map[f]].push_back(f);

    // find best merge
    MergedPair merged_set;
    FT min_merged_error = FT(0.0);
    bool first_merge = true;
    BOOST_FOREACH(edge_descriptor e, edges(*m_pmesh)) {
      if (CGAL::is_border(e, *m_pmesh))
        continue;
      std::size_t pxi = fproxy_map[face(halfedge(e, *m_pmesh), *m_pmesh)];
      std::size_t pxj = fproxy_map[face(opposite(halfedge(e, *m_pmesh), *m_pmesh), *m_pmesh)];
      if (pxi == pxj)
        continue;
      if (pxi > pxj)
        std::swap(pxi, pxj);
      if (merged_set.find(ProxyPair(pxi, pxj)) != merged_set.end())
        continue;

      std::list<face_descriptor> merged_patch(px_facets[pxi]);
      BOOST_FOREACH(face_descriptor f, px_facets[pxj])
        merged_patch.push_back(f);

      Proxy_wrapper pxw = fit_new_proxy(merged_patch.begin(), merged_patch.end());
      FT sum_error(0.0);
      BOOST_FOREACH(face_descriptor f, merged_patch)
        sum_error += (*fit_error)(f, pxw.px);
      merged_set.insert(ProxyPair(pxi, pxj));

      if (first_merge || sum_error < min_merged_error) {
        first_merge = false;
        min_merged_error = sum_error;
        px_enlarged = pxi;
        px_merged = pxj;
      }
    }

    // test if merge worth it
    if (if_test) {
      compute_fitting_error();
      FT max_error = proxies.front().err;
      for (std::size_t i = 0; i < proxies.size(); ++i) {
        if (max_error < proxies[i].err)
          max_error = proxies[i].err;
      }
      const FT merge_thre = max_error / FT(2.0);
      const FT increase = min_merged_error - (proxies[px_enlarged].err + proxies[px_merged].err);
      if (increase > merge_thre)
        return false;
    }

    return true;
  }

  /*!
   * @brief Split one proxy by default bisection, but N-section is also possible
   * No re-fitting performed and the proxy map is maintained.
   * @param px proxy index
   * @param n split section
   * @return change of error
   */
  FT split(const std::size_t px, const std::size_t n = 2) {
    if (px >= proxies.size())
      return FT(0.0);

    std::size_t count = 1;
    FT sum_err(0.0);
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
      if (count >= n)
        break;

      if (fproxy_map[f] == px && f != proxies[px].seed) {
        sum_err += (*fit_error)(f, proxies[px].px);
        fproxy_map[f] = proxies.size();
        proxies.push_back(fit_new_proxy(f));
        ++count;
      }
    }

    return sum_err;
  }

  /*!
   * @brief Extract the approximated surface mesh.
   * @note If the extracted surface mesh contains non-manifold facets, 
   * they are not built into the output polyhedron.
   * @tparam PolyhedronSurface should be `CGAL::Polyhedron_3`
   * @param[out] tm_out output triangle mesh
   * @param chord_error boundary approximation recursively split criterion
   * @param pca_plane true if use PCA plane fitting, otherwise use the default area averaged plane parameters
   * @return true if the extracted surface mesh is manifold, false otherwise.
   */
  template <typename PolyhedronSurface>
  bool extract_mesh(PolyhedronSurface &tm_out, const FT chord_error = FT(0.2), bool pca_plane = false) {
    // initialize all vertex anchor status
    enum Vertex_status { NO_ANCHOR = -1 };
    BOOST_FOREACH(vertex_descriptor v, vertices(*m_pmesh))
      internal_vidx_map[v] = static_cast<int>(NO_ANCHOR);
    anchors.clear();
    borders.clear();
    tris.clear();

    px_planes.clear();
    init_proxy_planes(pca_plane);

    find_anchors();
    find_edges(chord_error);
    add_anchors();
    pseudo_CDT();

    return build_polyhedron_surface(tm_out);
  }

  /*!
   * @brief Get the facet-proxy index map.
   * @tparam FacetProxyMap `WritablePropertyMap` with
   * `boost::graph_traits<TriangleMesh>::%face_descriptor` as key and `std::size_t` as value type
   * @param[out] facet_proxy_map facet proxy index map
   */
  template <typename FacetProxyMap>
  void get_proxy_map(FacetProxyMap &facet_proxy_map) const {
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      facet_proxy_map[f] = fproxy_map[f];
  }

  /*!
   * @brief Get the facet region of the specified proxy.
   * @tparam OutputIterator output iterator with `boost::graph_traits<TriangleMesh>::%face_descriptor` as value type
   * @param px_idx proxy index
   * @param out_itr output iterator
   */
  template <typename OutputIterator>
  void get_proxy_region(const std::size_t px_idx, OutputIterator out_itr) const {
    if (px_idx >= proxies.size())
      return;

    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      if (fproxy_map[f] == px_idx)
        *out_itr++ = f;
  }

  /*!
   * @brief Get the proxies.
   * @tparam OutputIterator output iterator with Proxy as value type
   * @param out_itr output iterator
   */
  template <typename OutputIterator>
  void get_proxies(OutputIterator out_itr) const {
    BOOST_FOREACH(const Proxy_wrapper &pxw, proxies)
      *out_itr++ = pxw.px;
  }

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  /*!
   * @brief Get the wrapped proxies.
   * @tparam OutputIterator output iterator with Proxy_wrapper as value type
   * @param out_itr output iterator
   */
  template <typename OutputIterator>
  void get_wrapped_proxies(OutputIterator out_itr) const {
    BOOST_FOREACH(const Proxy_wrapper &pxw, proxies)
      *out_itr++ = pxw;
  }
#endif

  /*!
   * @brief Get the proxies size.
   * @return number of proxies
   */
  std::size_t get_proxies_size() const { return proxies.size(); }

  /*!
   * @brief Get the anchor points, which have the area-averaged position of the projected anchor vertex points on the incident proxies.
   * @tparam OutputIterator output iterator with Point_3 as value type
   * @param out_itr output iterator
   */
  template <typename OutputIterator>
  void get_anchor_points(OutputIterator out_itr) const {
    BOOST_FOREACH(const Anchor &a, anchors)
      *out_itr++ = a.pos;
  }

  /*!
   * @brief Get the anchor vertices.
   * @tparam OutputIterator output iterator with vertex_descriptor as value type
   * @param out_itr output iterator
   */
  template <typename OutputIterator>
  void get_anchor_vertices(OutputIterator out_itr) const {
    BOOST_FOREACH(const Anchor &a, anchors)
      *out_itr++ = a.vtx;
  }

  /*!
   * @brief Get the indexed triangles,
   * one triplet of integers per triangles, which refers to the anchor point indexes.
   * @tparam OutputIterator output iterator with std::vector<size_t> as value type
   * @param out_itr output iterator
   */
  template <typename OutputIterator>
  void get_indexed_triangles(OutputIterator out_itr) const {
    BOOST_FOREACH(const std::vector<std::size_t> &t, tris)
      *out_itr++ = t;
  }

  /*!
   * @brief Get the indexed boundary polygon approximation.
   * @tparam OutputIterator output iterator with std::vector<std::size_t> as value type
   */
  template <typename OutputIterator>
  void get_indexed_boundary_polygons(OutputIterator out_itr) const {
    for (typename std::vector<Border>::const_iterator bitr = borders.begin();
      bitr != borders.end(); ++bitr) {
      std::vector<std::size_t> bdr;
      const halfedge_descriptor he_mark = bitr->he_head;
      halfedge_descriptor he = he_mark;
      do {
        Chord_vector chord;
        walk_to_next_anchor(he, chord);
        bdr.push_back(vanchor_map[target(he, *m_pmesh)]);
      } while (he != he_mark);
      *out_itr++ = bdr;
    }
  }

// private member functions
private:
  /*!
   * @brief Random initialize proxies to target number of proxies.
   * @note To ensure the randomness, call `std::srand()` beforehand.
   * @param max_nb_proxies maximum number of proxies, 
   * should be in range (0, num_faces(*m_pmesh))
   * @param num_iterations number of re-fitting iterations 
   * @return number of proxies initialized
   */
  std::size_t init_random(const std::size_t max_nb_proxies,
    const std::size_t num_iterations) {
    // fill a temporary vector of facets
    std::vector<face_descriptor> facets;
    random_shuffle_facets(facets);

    proxies.clear();
    // reach to the number of proxies
    for (std::size_t i = 0; i < max_nb_proxies; ++i)
      proxies.push_back(fit_new_proxy(facets[i]));
    run(num_iterations);

    return proxies.size();
  }

  /*!
   * @brief Incremental initialize proxies to target number of proxies.
   * @param max_nb_proxies maximum number of proxies, 
   * should be in range (0, num_faces(*m_pmesh))
   * @param num_iterations number of re-fitting iterations 
   * before each incremental proxy insertion
   * @return number of proxies initialized
   */
  std::size_t init_incremental(const std::size_t max_nb_proxies,
    const std::size_t num_iterations) {
    // initialize a proxy and the proxy map to prepare for the insertion
    bootstrap_from_first_facet();

    add_proxies_furthest(max_nb_proxies - 1, num_iterations);

    return proxies.size();
  }

  /*!
   * @brief Hierarchical initialize proxies to target number of proxies.
   * @param max_nb_proxies maximum number of proxies, 
   * should be in range (0, num_faces(*m_pmesh))
   * @param num_iterations number of re-fitting iterations
   * before each hierarchical proxy insertion
   * @return number of proxies initialized
   */
  std::size_t init_hierarchical(const std::size_t max_nb_proxies,
    const std::size_t num_iterations) {
    // initialize a proxy and the proxy map to prepare for the insertion
    bootstrap_from_first_facet();

    while (proxies.size() < max_nb_proxies) {
      // try to double current number of proxies each time
      std::size_t target_px = proxies.size();
      if (target_px * 2 > max_nb_proxies)
        target_px = max_nb_proxies;
      else
        target_px *= 2;
      // add proxies by error diffusion
      add_proxies_error_diffusion(target_px - proxies.size());
      run(num_iterations);
    }

    return proxies.size();
  }

  /*!
   * @brief Randomly initialize proxies
   * with both maximum number of proxies and minimum error drop stop criteria,
   * The first criterion met stops the seeding.
   * @note To ensure the randomness, call `std::srand()` beforehand.
   * @param max_nb_proxies maximum number of proxies, should be in range (0, num_faces(_mesh) / 3)
   * @param min_error_drop minimum error drop, should be in range (0.0, 1.0)
   * @param num_iterations number of re-fitting iterations 
   * @return number of proxies initialized
   */
  std::size_t init_random_error(const std::size_t max_nb_proxies,
    const FT min_error_drop,
    const std::size_t num_iterations) {
    // fill a temporary vector of facets
    std::vector<face_descriptor> facets;
    random_shuffle_facets(facets);

    bootstrap_from_first_facet();
    const FT initial_err = compute_fitting_error();
    FT error_drop = min_error_drop * FT(2.0);
    while (proxies.size() < max_nb_proxies && error_drop > min_error_drop) {
      // try to double current number of proxies each time
      std::size_t target_px = proxies.size();
      if (target_px * 2 > max_nb_proxies)
        target_px = max_nb_proxies;
      else
        target_px *= 2;
      proxies.clear();
      for (std::size_t j = 0; j < target_px; ++j)
        proxies.push_back(fit_new_proxy(facets[j]));
      run(num_iterations);
      const FT err = compute_fitting_error();
      error_drop = err / initial_err;
    }

    return proxies.size();
  }

  /*!
   * @brief Incrementally initialize proxies
   * with both maximum number of proxies and minimum error drop stop criteria,
   * The first criterion met stops the seeding.
   * @param max_nb_proxies maximum number of proxies, should be in range (0, num_faces(_mesh) / 3)
   * @param min_error_drop minimum error drop, should be in range (0.0, 1.0)
   * @param num_iterations number of re-fitting iterations 
   * @return number of proxies initialized
   */
  std::size_t init_incremental_error(const std::size_t max_nb_proxies,
    const FT min_error_drop,
    const std::size_t num_iterations) {
    // initialize a proxy and the proxy map to prepare for the insertion
    bootstrap_from_first_facet();
    const FT initial_err = compute_fitting_error();
    FT error_drop = min_error_drop * FT(2.0);
    while (proxies.size() < max_nb_proxies && error_drop > min_error_drop) {
      add_proxy_furthest();
      run(num_iterations);
      const FT err = compute_fitting_error();
      error_drop = err / initial_err;
    }

    return proxies.size();
  }

  /*!
   * @brief Hierarchically initialize proxies
   * with both maximum number of proxies and minimum error drop stop criteria,
   * The first criterion met stops the seeding.
   * @param max_nb_proxies maximum number of proxies, should be in range (0, num_faces(_mesh) / 3)
   * @param min_error_drop minimum error drop, should be in range (0.0, 1.0)
   * @param num_iterations number of re-fitting iterations 
   * @return number of proxies initialized
   */
  std::size_t init_hierarchical_error(const std::size_t max_nb_proxies,
    const FT min_error_drop,
    const std::size_t num_iterations) {
    // initialize a proxy and the proxy map to prepare for the insertion
    bootstrap_from_first_facet();
    const FT initial_err = compute_fitting_error();
    FT error_drop = min_error_drop * FT(2.0);
    while (proxies.size() < max_nb_proxies && error_drop > min_error_drop) {
      // try to double current number of proxies each time
      std::size_t target_px = proxies.size();
      if (target_px * 2 > max_nb_proxies)
        target_px = max_nb_proxies;
      else
        target_px *= 2;
      add_proxies_error_diffusion(target_px - proxies.size());
      run(num_iterations);
      const FT err = compute_fitting_error();
      error_drop = err / initial_err;
    }

    return proxies.size();
  }

  /*!
   * @brief Add a proxy seed at the facet with the maximum fitting error.
   * @pre current facet proxy map is valid, proxy error is computed
   * @note No re-fitting is performed. After the operation, the facet proxy map remains valid.
   * @return true add successfully, false otherwise
   */
  bool add_proxy_furthest() {
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
    std::cerr << "add furthest " << proxies.size() << std::endl;
#endif
    compute_fitting_error();
    FT max_error = proxies.front().err;
    std::size_t px_worst = 0;
    for (std::size_t i = 0; i < proxies.size(); ++i) {
      if (max_error < proxies[i].err) {
        max_error = proxies[i].err;
        px_worst = i;
      }
    }

    face_descriptor fworst;
    bool first = true;
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
      std::size_t px_idx = fproxy_map[f];
      if (px_idx != px_worst || f == proxies[px_idx].seed)
        continue;

      FT err = (*fit_error)(f, proxies[px_idx].px);
      if (first || max_error < err) {
        first = false;
        max_error = err;
        fworst = f;
      }
    }

    if (first)
      return false;

    fproxy_map[fworst] = proxies.size();
    proxies.push_back(fit_new_proxy(fworst));

    return true;
  }

  /*!
   * @brief Fitting a new proxy.
   * 1. Fit proxy parameters from a list of facets.
   * 2. Set seed.
   * @tparam FacetIterator face_descriptor container iterator
   * @param beg container begin
   * @param end container end
   */
  template<typename FacetIterator>
  Proxy_wrapper fit_new_proxy(const FacetIterator &beg, const FacetIterator &end) {
    CGAL_assertion(beg != end);

    // use proxy_fitting functor to fit proxy parameters
    Proxy px = (*proxy_fitting)(beg, end);

    // find proxy seed
    face_descriptor seed = *beg;
    FT err_min = (*fit_error)(*beg, px);
    std::pair<FacetIterator, FacetIterator> facets(beg, end);
    BOOST_FOREACH(face_descriptor f, facets) {
      FT err = (*fit_error)(f, px);
      if (err < err_min) {
        err_min = err;
        seed = f;
      }
    }

    return Proxy_wrapper(px, seed);
  }

  /*!
   * @brief Fitting a new proxy from a single facet.
   * 1. Fit proxy parameters from one facet.
   * 2. Set seed.
   * @param face_descriptor facet
   */
  Proxy_wrapper fit_new_proxy(const face_descriptor &f) {
    std::vector<face_descriptor> fvec(1, f);
    // fit proxy parameters
    Proxy px = (*proxy_fitting)(fvec.begin(), fvec.end());

    return Proxy_wrapper(px, f);
  }

  /*!
   * @brief Random shuffle the surface facets into an empty vector.
   */
  void random_shuffle_facets(std::vector<face_descriptor> &facets) {
    const std::size_t nbf = num_faces(*m_pmesh);
    facets.reserve(nbf);
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      facets.push_back(f);
    // random shuffle
    for (std::size_t i = 0; i < nbf; ++i) {
      // swap ith element with a random one
      std::size_t r = static_cast<std::size_t>(
        static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX) * 
        static_cast<double>(nbf - 1));
      face_descriptor tmp = facets[r];
      facets[r] = facets[i];
      facets[i] = tmp;
    }
  }

  /*!
   * @brief Initialize a proxy from the first facet of the surface.
   * @note This function clears proxy vector and set facet proxy map to initial state,
   * intended only for bootstrapping initialization.
   * Coarse approximation iteration is not performed, because it's inaccurate anyway
   * and may cause serious degenerate cases(e.g. a standard cube mode).
   */
  void bootstrap_from_first_facet() {
    proxies.clear();
    proxies.push_back(fit_new_proxy(*(faces(*m_pmesh).first)));
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      fproxy_map[f] = 0;
  }

  /*!
   * @brief Initialize proxy planes.
   * @param if_pca_plane true to use the PCA plane fitting
   */
  void init_proxy_planes(const bool if_pca_plane) {
    // fit proxy planes, areas, normals
    std::vector<std::list<face_descriptor> > px_facets(proxies.size());
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      px_facets[fproxy_map[f]].push_back(f);

    BOOST_FOREACH(const std::list<face_descriptor> &px_patch, px_facets) {
      Plane_3 fit_plane = if_pca_plane ? 
        fit_plane_pca(px_patch.begin(), px_patch.end()) :
          fit_plane_area_averaged(px_patch.begin(), px_patch.end());

      Vector_3 norm = CGAL::NULL_VECTOR;
      FT area(0.0);
      BOOST_FOREACH(face_descriptor f, px_patch) {
        halfedge_descriptor he = halfedge(f, *m_pmesh);
        const Point_3 &p0 = point_pmap[source(he, *m_pmesh)];
        const Point_3 &p1 = point_pmap[target(he, *m_pmesh)];
        const Point_3 &p2 = point_pmap[target(next(he, *m_pmesh), *m_pmesh)];
        const FT farea(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
		const Vector_3 fnorm = CGAL::unit_normal(p0, p1, p2);

        norm = sum_functor(norm, scale_functor(fnorm, farea));
        area += farea;
      }
      norm = scale_functor(norm, FT(1.0 / std::sqrt(CGAL::to_double(norm.squared_length()))));

      px_planes.push_back(Proxy_plane(fit_plane, norm, area));
    }
  }

  /*!
   * @brief Finds the anchors.
   */
  void find_anchors() {
    BOOST_FOREACH(vertex_descriptor vtx, vertices(*m_pmesh)) {
      std::size_t border_count = 0;

      BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(vtx, *m_pmesh)) {
        if (CGAL::is_border_edge(h, *m_pmesh))
          ++border_count;
        else if (fproxy_map[face(h, *m_pmesh)] != fproxy_map[face(opposite(h, *m_pmesh), *m_pmesh)])
          ++border_count;
      }
      if (border_count >= 3)
        attach_anchor(vtx);
    }
  }

  /*!
   * @brief Finds and approximates the chord connecting the anchors.
   * @param chord_error boundary chord approximation recursive split creterion
   */
  void find_edges(const FT chord_error) {
    // collect candidate halfedges in a set
    std::set<halfedge_descriptor> he_candidates;
    BOOST_FOREACH(halfedge_descriptor h, halfedges(*m_pmesh)) {
      if (!CGAL::is_border(h, *m_pmesh)
        && (CGAL::is_border(opposite(h, *m_pmesh), *m_pmesh)
          || fproxy_map[face(h, *m_pmesh)] != fproxy_map[face(opposite(h, *m_pmesh), *m_pmesh)]))
        he_candidates.insert(h);
    }

    // pick up one candidate halfedge each time and traverse the connected border
    while (!he_candidates.empty()) {
      halfedge_descriptor he_start = *he_candidates.begin();
      walk_to_first_anchor(he_start);
      // no anchor in this connected border, make a new anchor
      if (!is_anchor_attached(he_start))
        attach_anchor(he_start);

      // a new connected border
      borders.push_back(Border(he_start));

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
      std::cerr << "#border " << borders.size() << std::endl;
#endif

      const halfedge_descriptor he_mark = he_start;
      do {
        Chord_vector chord;
        walk_to_next_anchor(he_start, chord);
        borders.back().num_anchors += subdivide_chord(chord.begin(), chord.end(), chord_error);

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
        std::cerr << "#chord_anchor " << borders.back().num_anchors << std::endl;
#endif

        for (Chord_vector_iterator citr = chord.begin(); citr != chord.end(); ++citr)
          he_candidates.erase(*citr);
      } while (he_start != he_mark);
    }
  }

  /*!
   * @brief Adds anchors to the border cycles with only 2 anchors.
   */
  void add_anchors() {
    typedef typename std::vector<Border>::iterator BorderIterator;
    for (BorderIterator bitr = borders.begin(); bitr != borders.end(); ++bitr) {
      if (bitr->num_anchors > 2)
        continue;

      // 2 initial anchors at least
      CGAL_assertion(bitr->num_anchors == 2);
      // borders with only 2 initial anchors
      const halfedge_descriptor he_mark = bitr->he_head;
      Point_3 pt_begin = point_pmap[target(he_mark, *m_pmesh)];
      Point_3 pt_end = pt_begin;

      halfedge_descriptor he = he_mark;
      Chord_vector chord;
      std::size_t count = 0;
      do {
        walk_to_next_border_halfedge(he);
        if (!is_anchor_attached(he))
          chord.push_back(he);
        else {
          if (count == 0)
            pt_end = point_pmap[target(he, *m_pmesh)];
          ++count;
        }
      } while (he != he_mark);

      // anchor count may be increased to more than 2 afterwards
      // due to the new anchors added by the neighboring border (< 2 anchors)
      if (count > 2) {
        bitr->num_anchors = count;
        continue;
      }

      FT dist_max(0.0);
      halfedge_descriptor he_max;
      Vector_3 chord_vec = vector_functor(pt_begin, pt_end);
      chord_vec = scale_functor(chord_vec,
        FT(1.0 / std::sqrt(CGAL::to_double(chord_vec.squared_length()))));
      for (Chord_vector_iterator citr = chord.begin(); citr != chord.end(); ++citr) {
        Vector_3 vec = vector_functor(pt_begin, point_pmap[target(*citr, *m_pmesh)]);
        vec = CGAL::cross_product(chord_vec, vec);
        FT dist(std::sqrt(CGAL::to_double(vec.squared_length())));
        if (dist > dist_max) {
          dist_max = dist;
          he_max = *citr;
        }
      }
      attach_anchor(he_max);
      // increase border anchors by one
      bitr->num_anchors++;
    }
  }

  /*!
   * @brief Runs the pseudo Constrained Delaunay Triangulation at each region, and stores the extracted indexed triangles in @a tris.
   */
  void pseudo_CDT() {
    // subgraph attached with vertex anchor status and edge weight
    typedef boost::property<boost::vertex_index1_t, int,
      boost::property<boost::vertex_index2_t, int> > VertexProperty;
    typedef boost::property<boost::edge_weight_t, FT,
      boost::property<boost::edge_index_t, int> > EdgeProperty;
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
    BOOST_FOREACH(vertex_descriptor v, vertices(*m_pmesh)) {
      sg_vertex_descriptor sgv = add_vertex(gmain);
      global_vanchor_map[sgv] = vanchor_map[v];
      global_vtag_map[sgv] = vanchor_map[v];
      vmap.insert(std::pair<vertex_descriptor, sg_vertex_descriptor>(v, sgv));
    }
    BOOST_FOREACH(edge_descriptor e, edges(*m_pmesh)) {
      vertex_descriptor vs = source(e, *m_pmesh);
      vertex_descriptor vt = target(e, *m_pmesh);
      FT len(std::sqrt(CGAL::to_double(
        CGAL::squared_distance(point_pmap[vs], point_pmap[vt]))));
      add_edge(to_sgv_map[vs], to_sgv_map[vt], len, gmain);
    }

    std::vector<VertexVector> vertex_patches(proxies.size());
    BOOST_FOREACH(vertex_descriptor v, vertices(*m_pmesh)) {
      std::set<std::size_t> px_set;
      BOOST_FOREACH(face_descriptor f, faces_around_target(halfedge(v, *m_pmesh), *m_pmesh)) {
        if (f != boost::graph_traits<TriangleMesh>::null_face())
          px_set.insert(fproxy_map[f]);
      }
      BOOST_FOREACH(std::size_t p, px_set)
        vertex_patches[p].push_back(to_sgv_map[v]);
    }
    BOOST_FOREACH(VertexVector &vpatch, vertex_patches) {
      // add a super vertex connecting to its boundary anchors into the main graph
      const sg_vertex_descriptor superv = add_vertex(gmain);
      global_vanchor_map[superv] = 0;
      global_vtag_map[superv] = 0;
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
        sg_vertex_descriptor curr = v;
        while (!is_anchor_attached(curr, local_vanchor_map))
          curr = pred[curr];
        local_vtag_map[v] = local_vanchor_map[curr];
      }
    }

    // tag all boundary chord
    BOOST_FOREACH(const Border &bdr, borders) {
      const halfedge_descriptor he_mark = bdr.he_head;
      halfedge_descriptor he = he_mark;
      do {
        Chord_vector chord;
        walk_to_next_anchor(he, chord);

        std::vector<FT> vdist;
        vdist.push_back(FT(0.0));
        BOOST_FOREACH(halfedge_descriptor h, chord) {
          FT elen = global_eweight_map[edge(
            to_sgv_map[source(h, *m_pmesh)],
            to_sgv_map[target(h, *m_pmesh)],
            gmain).first];
          vdist.push_back(vdist.back() + elen);
        }

        FT half_chord_len = vdist.back() / FT(2.0);
        const int anchorleft = vanchor_map[source(chord.front(), *m_pmesh)];
        const int anchorright = vanchor_map[target(chord.back(), *m_pmesh)];
        typename std::vector<FT>::iterator ditr = vdist.begin() + 1;
        for (typename Chord_vector::iterator hitr = chord.begin();
          hitr != chord.end() - 1; ++hitr, ++ditr) {
          if (*ditr < half_chord_len)
            global_vtag_map[to_sgv_map[target(*hitr, *m_pmesh)]] = anchorleft;
          else
            global_vtag_map[to_sgv_map[target(*hitr, *m_pmesh)]] = anchorright;
        }
      } while (he != he_mark);
    }

    // collect triangles
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
      halfedge_descriptor he = halfedge(f, *m_pmesh);
      int i = global_vtag_map[to_sgv_map[source(he, *m_pmesh)]];
      int j = global_vtag_map[to_sgv_map[target(he, *m_pmesh)]];
      int k = global_vtag_map[to_sgv_map[target(next(he, *m_pmesh), *m_pmesh)]];
      if (i != j && i != k && j != k) {
        std::vector<std::size_t> t;
        t.push_back(static_cast<std::size_t>(i));
        t.push_back(static_cast<std::size_t>(j));
        t.push_back(static_cast<std::size_t>(k));
        tris.push_back(t);
      }
    }
  }

  /*!
   * @brief Walks along the region border to the first halfedge pointing to a vertex associated with an anchor.
   * @param[in/out] he_start region border halfedge
   */
  void walk_to_first_anchor(halfedge_descriptor &he_start) {
    const halfedge_descriptor start_mark = he_start;
    while (!is_anchor_attached(he_start)) {
      // no anchor attached to the halfedge target
      walk_to_next_border_halfedge(he_start);
      if (he_start == start_mark) // back to where started, a circular border
        return;
    }
  }

  /*!
   * @brief Walks along the region border to the next anchor and records the path as @a chord.
   * @param[in/out] he_start starting region border halfedge pointing to a vertex associated with an anchor
   * @param[out] chord recorded path chord
   */
  void walk_to_next_anchor(halfedge_descriptor &he_start, Chord_vector &chord) const {
    do {
      walk_to_next_border_halfedge(he_start);
      chord.push_back(he_start);
    } while (!is_anchor_attached(he_start));
  }

  /*!
   * @brief Walks to next border halfedge.
   * @param[in/out] he_start region border halfedge
   */
  void walk_to_next_border_halfedge(halfedge_descriptor &he_start) const {
    const std::size_t px_idx = fproxy_map[face(he_start, *m_pmesh)];
    BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(he_start, *m_pmesh)) {
      if (CGAL::is_border(h, *m_pmesh) || fproxy_map[face(h, *m_pmesh)] != px_idx) {
        he_start = opposite(h, *m_pmesh);
        return;
      }
    }
  }

  /*!
   * @brief Subdivides a chord recursively in range [@a chord_begin, @a chord_end).
   * @param chord_begin begin iterator of the chord
   * @param chord_end end iterator of the chord
   * @param the recursive split threshold
   * @return the number of anchors of the chord apart from the first one
   */
  std::size_t subdivide_chord(
    const Chord_vector_iterator &chord_begin,
    const Chord_vector_iterator &chord_end,
    const FT thre) {
    const std::size_t chord_size = std::distance(chord_begin, chord_end);
    const halfedge_descriptor he_first = *chord_begin;
    const halfedge_descriptor he_last = *(chord_end - 1);
    const std::size_t anchor_first = vanchor_map[source(he_first, *m_pmesh)];
    const std::size_t anchor_last = vanchor_map[target(he_last, *m_pmesh)];

    // do not subdivide trivial non-circular chord
    if ((anchor_first != anchor_last) && (chord_size < 4))
      return 1;

    bool if_subdivide = false;
    Chord_vector_iterator chord_max;
    const Point_3 &pt_begin = point_pmap[source(he_first, *m_pmesh)];
    const Point_3 &pt_end = point_pmap[target(he_last, *m_pmesh)];
    if (anchor_first == anchor_last) {
      // circular chord
      CGAL_assertion(chord_size > 2);

      FT dist_max(0.0);
      for (Chord_vector_iterator citr = chord_begin; citr != chord_end; ++citr) {
        FT dist = CGAL::squared_distance(pt_begin, point_pmap[target(*citr, *m_pmesh)]);
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

      for (Chord_vector_iterator citr = chord_begin; citr != chord_end; ++citr) {
        Vector_3 vec = vector_functor(pt_begin, point_pmap[target(*citr, *m_pmesh)]);
        vec = CGAL::cross_product(chord_vec, vec);
        FT dist(std::sqrt(CGAL::to_double(vec.squared_length())));
        if (dist > dist_max) {
          chord_max = citr;
          dist_max = dist;
        }
      }

      // suppose the proxy normal angle is acute
      std::size_t px_left = fproxy_map[face(he_first, *m_pmesh)];
      std::size_t px_right = px_left;
      if (!CGAL::is_border(opposite(he_first, *m_pmesh), *m_pmesh))
        px_right = fproxy_map[face(opposite(he_first, *m_pmesh), *m_pmesh)];
      FT norm_sin(1.0);
      if (!CGAL::is_border(opposite(he_first, *m_pmesh), *m_pmesh)) {
        Vector_3 vec = CGAL::cross_product(
          px_planes[px_left].normal, px_planes[px_right].normal);
        norm_sin = FT(std::sqrt(CGAL::to_double(scalar_product_functor(vec, vec))));
      }
      FT criterion = dist_max * norm_sin / chord_len;
      if (criterion > thre)
        if_subdivide = true;
    }

    if (if_subdivide) {
      // subdivide at the most remote vertex
      attach_anchor(*chord_max);

      std::size_t num0 = subdivide_chord(chord_begin, chord_max + 1, thre);
      std::size_t num1 = subdivide_chord(chord_max + 1, chord_end, thre);

      return num0 + num1;
    }

    return 1;
  }

  /*!
   * @brief Return true if the target vertex of a halfedge is attached with an anchor, and false otherwise.
   * @param he halfedge
   */
  bool is_anchor_attached(const halfedge_descriptor &he) const {
    return is_anchor_attached(target(he, *m_pmesh), vanchor_map);
  }

  /*!
   * @brief Check if a vertex is attached with an anchor.
   * @tparam VertexAnchorIndexMap `WritablePropertyMap` with `boost::graph_traights<TriangleMesh>::vertex_descriptor` as key and `std::size_t` as value type
   * @param v vertex
   * @param vertex_anchor_map vertex anchor index map
   */
  template<typename VertexAnchorIndexMap>
  bool is_anchor_attached(
    const typename boost::property_traits<VertexAnchorIndexMap>::key_type &v,
    const VertexAnchorIndexMap &vertex_anchor_map) const {
    return vertex_anchor_map[v] >= 0;
  }

  /*!
   * @brief Attachs an anchor to the vertex.
   * @param vtx vertex
   */
  void attach_anchor(const vertex_descriptor &vtx) {
    vanchor_map[vtx] = static_cast<int>(anchors.size());
    anchors.push_back(Anchor(vtx, compute_anchor_position(vtx)));
  }

  /*!
   * @brief Attachs an anchor to the target vertex of the halfedge.
   * @param he halfedge
   */
  void attach_anchor(const halfedge_descriptor &he) {
    vertex_descriptor vtx = target(he, *m_pmesh);
    attach_anchor(vtx);
  }

  /*!
   * @brief Calculate the anchor positions from a vertex.
   * @param v the vertex descriptor
   * @return the anchor position
   */
  Point_3 compute_anchor_position(const vertex_descriptor &v) {
    // construct an anchor from vertex and the incident proxies
    std::set<std::size_t> px_set;
    BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(v, *m_pmesh)) {
      if (!CGAL::is_border(h, *m_pmesh))
        px_set.insert(fproxy_map[face(h, *m_pmesh)]);
    }

    // construct an anchor from vertex and the incident proxies
    FT avgx(0.0), avgy(0.0), avgz(0.0), sum_area(0.0);
    const Point_3 vtx_pt = point_pmap[v];
    for (std::set<std::size_t>::iterator pxitr = px_set.begin();
      pxitr != px_set.end(); ++pxitr) {
      std::size_t px_idx = *pxitr;
      Point_3 proj = px_planes[px_idx].plane.projection(vtx_pt);
      FT area = px_planes[px_idx].area; // TOFIX: use Vector and CGAL::ORIGIN + vector
      avgx += proj.x() * area;
      avgy += proj.y() * area;
      avgz += proj.z() * area;
      sum_area += area;
    }
    return Point_3(avgx / sum_area, avgy / sum_area, avgz / sum_area);
  }

  /*!
   * @brief Use an incremental builder to build and test if the indexed triangle surface is manifold
   * @tparam PolyhedronSurface should be `CGAL::Polyhedron_3`
   * @param[out] poly input polyhedorn mesh
   * @return true if build manifold surface successfully
   */
  template <typename PolyhedronSurface>
  bool build_polyhedron_surface(PolyhedronSurface &poly) {
    std::vector<Point_3> vtx;
    BOOST_FOREACH(const Anchor &a, anchors)
      vtx.push_back(a.pos);

    typedef typename PolyhedronSurface::HalfedgeDS HDS;
    Triangle_polyhedron_builder<HDS> tpbuilder(vtx, tris);
    poly.delegate(tpbuilder);

    return tpbuilder.is_manifold;
  }

  /*!
   * @brief Fit an area averaged plane from a range of facets.
   * @tparam FacetIterator face_descriptor container iterator
   * @param beg container begin
   * @param end container end
   * @return fitted plane
   */
  template <typename FacetIterator>
  Plane_3 fit_plane_area_averaged(const FacetIterator &beg, const FacetIterator &end) {
    CGAL_assertion(beg != end);
    // area average normal and centroid
    Vector_3 norm = CGAL::NULL_VECTOR;
    Vector_3 cent = CGAL::NULL_VECTOR;
    FT sum_area(0.0);
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      const halfedge_descriptor he = halfedge(*fitr, *m_pmesh);
      const Point_3 &p0 = point_pmap[source(he, *m_pmesh)];
      const Point_3 &p1 = point_pmap[target(he, *m_pmesh)];
      const Point_3 &p2 = point_pmap[target(next(he, *m_pmesh), *m_pmesh)];

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
   * @brief Fit a plane from a range of facets with PCA algorithm.
   * @tparam FacetIterator face_descriptor container iterator
   * @param beg container begin
   * @param end container end
   * @return fitted plane
   */
  template <typename FacetIterator>
  Plane_3 fit_plane_pca(const FacetIterator &beg, const FacetIterator &end) {
    CGAL_assertion(beg != end);

    typedef typename Geom_traits::Triangle_3 Triangle_3;
    std::list<Triangle_3> tri_list;
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      halfedge_descriptor he = halfedge(*fitr, *m_pmesh);
      const Point_3 &p0 = point_pmap[source(he, *m_pmesh)];
      const Point_3 &p1 = point_pmap[target(he, *m_pmesh)];
      const Point_3 &p2 = point_pmap[target(next(he, *m_pmesh), *m_pmesh)];
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

} // end namespace VSA
} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_APPROXIMATION_H
