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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>

#include <vector>
#include <stack>
#include <queue>
#include <iterator>
#include <cmath>
#include <cstdlib>

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
#include <iostream>
#endif

#define CGAL_VSA_INVALID_TAG std::numeric_limits<std::size_t>::max()

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
    Proxy_wrapper(const Proxy &p, const std::size_t &i, const face_descriptor &s, const FT &e)
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

  // The facet candidate to be queued.
  struct Facet_to_integrate {
    Facet_to_integrate(const face_descriptor &f_, const std::size_t &px_, const FT &err_)
      : f(f_), px(px_), err(err_) {}

    bool operator<(const Facet_to_integrate &rhs) const {
      return err > rhs.err;
    }

    face_descriptor f; // facet
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

    halfedge_descriptor he_head; // The heading halfedge of the boundary cylce.
    std::size_t num_anchors; // The number of anchors on the boundary cycle.
  };

  // Triangle polyhedron builder.
  template <typename HDS>
  class Triangle_polyhedron_builder : public CGAL::Modifier_base<HDS> {
    const std::vector<Point_3> &vtxs;
    const std::vector<std::vector<std::size_t> > &tris;
  public:
    bool is_manifold;
    Triangle_polyhedron_builder(const std::vector<Point_3> &vtxs_,
      const std::vector<std::vector<std::size_t> > &tris_)
      : vtxs(vtxs_), tris(tris_), is_manifold(true) {}

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
  const TriangleMesh *m_ptm;
  // The mesh vertex point map.
  VertexPointMap m_vpoint_map;
  // The error metric functor.
  const Error_metric *m_perror_metric;
  // The proxy fitting functor.
  const Proxy_fitting *m_pproxy_fitting;

  Construct_vector_3 vector_functor;
  Construct_point_3 point_functor;
  Construct_scaled_vector_3 scale_functor;
  Construct_sum_of_vectors_3 sum_functor;
  Compute_scalar_product_3 scalar_product_functor;
  Construct_translated_point_3 translate_point_functor;

  // The facet proxy index map.
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
  std::vector<std::vector<std::size_t> > m_tris;

  // meshing parameters
  FT m_average_edge_length;

//member functions
public:
  /*!
   * %Default constructor.
   */
  Mesh_approximation() :
    m_ptm(NULL),
    m_perror_metric(NULL),
    m_pproxy_fitting(NULL),
    m_average_edge_length(0.0) {

    Geom_traits traits;
    vector_functor = traits.construct_vector_3_object();
    point_functor = traits.construct_point_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scalar_product_functor = traits.compute_scalar_product_3_object();
    translate_point_functor = traits.construct_translated_point_3_object();
  }

  /*!
   * Initialize and prepare for the approximation.
   * @param tm `CGAL TriangleMesh` on which approximation operates.
   * @param vpoint_map vertex point map of the mesh
   */
  Mesh_approximation(const TriangleMesh &tm, const VertexPointMap &vpoint_map) :
    m_ptm(&tm),
    m_vpoint_map(vpoint_map),
    m_perror_metric(NULL),
    m_pproxy_fitting(NULL),
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

  
  ~Mesh_approximation() {
    if (m_ptm) {
      CGAL::internal::remove_property(m_vanchor_map, *(const_cast<TriangleMesh *>(m_ptm)));
      CGAL::internal::remove_property(m_fproxy_map, *(const_cast<TriangleMesh *>(m_ptm)));
    }
  }
  
  /*!
   * Set the mesh for approximation and rebuild the internal data structure.
   * @pre @a tm.is_pure_triangle()
   * @param tm `CGAL TriangleMesh` on which approximation operates.
   * @param vpoint_map vertex point map of the mesh
   */
  void set_mesh(const TriangleMesh &tm, const VertexPointMap &vpoint_map) {
    m_ptm = &tm;
    m_vpoint_map = vpoint_map;
    rebuild();
  }

  /*!
   * Set the error and fitting functor.
   * @param error_metric_ an `ErrorMetric` functor.
   * @param proxy_fitting_ a `ProxyFitting` functor.
   */
  void set_metric(const Error_metric &error_metric_,
    const Proxy_fitting &proxy_fitting_) {
    m_perror_metric = &error_metric_;
    m_pproxy_fitting = &proxy_fitting_;
  }

  /*!
   * Rebuild the internal data structure.
   */
  void rebuild() {
    // cleanup
    m_proxies.clear();
    m_px_planes.clear();
    m_anchors.clear();
    m_bcycles.clear();
    m_tris.clear();

    if (!m_ptm)
      return;

    // rebuild internal data structure
    CGAL::internal::remove_property(m_fproxy_map, *m_ptm);
    m_fproxy_map = CGAL::internal::add_property(
      Face_proxy_tag("VSA-face_proxy"), *(const_cast<TriangleMesh *>(m_ptm)));
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
      put(m_fproxy_map, f, CGAL_VSA_INVALID_TAG);

    CGAL::internal::remove_property(m_vanchor_map, *m_ptm);
    m_vanchor_map = CGAL::internal::add_property(
      Vertex_anchor_tag("VSA-vertex_anchor"), *(const_cast<TriangleMesh *>(m_ptm)));
    BOOST_FOREACH(vertex_descriptor v, vertices(*m_ptm))
      put(m_vanchor_map, v, CGAL_VSA_INVALID_TAG);
  }

  /*!
   * @brief Initialize the seeds with both maximum number of proxies
   * and minmum error drop stop criteria.
   * The first criterion met stops the seeding.
   * Parameters out of range are ignored.
   * @param method seeding method
   * @param max_nb_proxies maximum target number of proxies,
   * should be in range (nb_connected_components, num_faces(tm) / 3)
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
    const std::size_t nb_px = num_faces(*m_ptm) / 3;

    // initialize proxies and the proxy map to prepare for insertion
    bootstrap_from_connected_components();

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
   * @brief Run the partitioning and fitting processes on the whole surface.
   * @param nb_iterations number of iterations.
   */
  void run(std::size_t nb_iterations = 1) {
    for (std::size_t i = 0; i < nb_iterations; ++i) {
      // tag the whole surface
      BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
        put(m_fproxy_map, f, CGAL_VSA_INVALID_TAG);

      partition(m_proxies.begin(), m_proxies.end());
      fit(m_proxies.begin(), m_proxies.end());
    }
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
    static std::size_t count = 0;
    std::cerr << '#' << count++ << ": " << compute_fitting_error() << std::endl;
#endif
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
    BOOST_FOREACH(Proxy_wrapper &pxw, m_proxies)
      pxw.err = FT(0.0);
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm)) {
      std::size_t pxidx = get(m_fproxy_map, f);
      m_proxies[pxidx].err += (*m_perror_metric)(f, m_proxies[pxidx].px);
    }

    FT sum_error(0.0);
    BOOST_FOREACH(const Proxy_wrapper &pxw, m_proxies)
      sum_error += pxw.err;

    return sum_error;
  }

  /*!
   * @brief Add proxies to the worst regions one by one.
   * The re-fitting is performed after each proxy is inserted.
   * @pre current facet proxy map is valid
   * @note after the addition, the facet proxy map remains valid
   * @param num_proxies number of proxies to be added
   * @param nb_iterations number of re-fitting iterations
   * @return number of proxies added
   */
  std::size_t add_proxies_furthest(const std::size_t num_proxies,
    const std::size_t nb_iterations = 5) {
    std::size_t num_added = 0;
    for (; num_added < num_proxies; ++num_added) {
      if (!add_proxy_furthest())
        break;
      run(nb_iterations);
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
    std::cerr << "#px " << m_proxies.size() << std::endl;
#endif

    const FT sum_error = compute_fitting_error();
    const FT avg_error = sum_error / FT(static_cast<double>(num_proxies));

    std::vector<Proxy_error> px_error;
    for (std::size_t i = 0; i < m_proxies.size(); ++i)
      px_error.push_back(Proxy_error(i, m_proxies[i].err));
    // sort partition by error
    std::sort(px_error.begin(), px_error.end());

    // number of proxies to be added to each region
    std::vector<std::size_t> num_to_add(m_proxies.size(), 0);
    if (avg_error == FT(0.0)) {
      // rare case on extremely regular geometry like a cube
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
      std::cerr << "zero error, diffuse w.r.t. number of facets" << std::endl;
#endif
      const FT avg_facet = FT(
        static_cast<double>(num_faces(*m_ptm)) / static_cast<double>(num_proxies));
      std::vector<FT> px_size(m_proxies.size(), FT(0.0));
      BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
        px_size[get(m_fproxy_map, f)] += FT(1.0);
      FT residual(0.0);
      for (std::size_t i = 0; i < m_proxies.size(); ++i) {
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
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm)) {
      const std::size_t px_id = get(m_fproxy_map, f);
      if (m_proxies[px_id].seed == f)
        continue;

      if (num_to_add[px_id] > 0) {
        m_proxies.push_back(fit_new_proxy(f, m_proxies.size()));
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
   * @param if_force true to force the teleportation (no merge test)
   * @return number of proxies teleported.
   */
  std::size_t teleport_proxies(const std::size_t num_proxies,
    const std::size_t num_iterations = 5,
    const bool if_force = false) {
    std::size_t num_teleported = 0;
    while (num_teleported < num_proxies) {
      // find worst proxy
      compute_fitting_error();
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

      // teleport to a facet of the worst region
      // update merged proxies
      std::list<face_descriptor> merged_patch;
      BOOST_FOREACH(face_descriptor f, faces(*m_ptm)) {
        std::size_t px_idx = get(m_fproxy_map, f);
        if (px_idx == px_enlarged || px_idx == px_merged) {
          put(m_fproxy_map, f, px_enlarged);
          merged_patch.push_back(f);
        }
      }
      m_proxies[px_enlarged] = fit_new_proxy(merged_patch.begin(), merged_patch.end(), px_enlarged);
      // replace the merged proxy position to the newly teleported proxy
      m_proxies[px_merged] = fit_new_proxy(tele_to, px_merged);
      put(m_fproxy_map, tele_to, px_merged);

      num_teleported++;
      // coarse re-fitting
      run(num_iterations);

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
      std::cerr << "teleported" << std::endl;
#endif
    }

    return num_teleported;
  }

  /*!
   * @brief Merge two specified adjacent regions.
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

    // merge px1 to px0
    FT err_sum(0.0);
    std::list<face_descriptor> merged_patch;
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm)) {
      std::size_t px_idx = get(m_fproxy_map, f);
      if (px_idx == px1 || px_idx == px0) {
        err_sum += (*m_perror_metric)(f, m_proxies[px_idx].px);
        put(m_fproxy_map, f, px0);
        merged_patch.push_back(f);
      }
    }
    m_proxies[px0] = fit_new_proxy(merged_patch.begin(), merged_patch.end(), px0);

    // erase px1 and maintain proxy index
    m_proxies.erase(m_proxies.begin() + px1);
    for (std::size_t i = 0; i < m_proxies.size(); ++i)
      m_proxies[i].idx = i;
    // keep facet proxy map valid
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm)) {
      if (get(m_fproxy_map, f) > px1)
        put(m_fproxy_map, f, get(m_fproxy_map, f) - 1);
    }

    FT err_merged(0.0);
    BOOST_FOREACH(face_descriptor f, merged_patch)
      err_merged += (*m_perror_metric)(f, m_proxies[px0].px);

    return err_merged - err_sum;
  }

  /*!
   * @brief Simulate merging and local re-fitting of all two adjacent proxies
   * and find the best two regions to merge.
   * @note The <b>best</b> is defined as the minimum merged sum error
   * <b>change</b> (increase of decrease) among all adjacent pairs.
   * @param[out] px_tobe_enlarged the proxy index to be enlarged
   * @param[out] px_tobe_merged the proxy index to be merged,
   * guaranteed to be greater than <em>px_tobe_enlarged</em>.
   * @param if_test true to activate the merge test.
   * The merge test is considered successful if the merged error change
   * is less than the half of the maximum proxy error.
   * @return true if best merge pair found, false otherwise
   */
  bool find_best_merge(std::size_t &px_tobe_enlarged,
    std::size_t &px_tobe_merged,
    const bool if_test) {
    typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
    typedef std::pair<std::size_t, std::size_t> ProxyPair;
    typedef std::set<ProxyPair> MergedPair;

    compute_fitting_error();

    std::vector<std::list<face_descriptor> > px_facets(m_proxies.size());
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
      px_facets[get(m_fproxy_map, f)].push_back(f);

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
      std::list<face_descriptor> merged_patch(px_facets[pxi]);
      BOOST_FOREACH(face_descriptor f, px_facets[pxj])
        merged_patch.push_back(f);
      Proxy_wrapper pxw_tmp = fit_new_proxy(
        merged_patch.begin(), merged_patch.end(), CGAL_VSA_INVALID_TAG);

      FT error_merged(0.0);
      BOOST_FOREACH(face_descriptor f, merged_patch)
        error_merged += (*m_perror_metric)(f, pxw_tmp.px);
      const FT error_change = error_merged - (m_proxies[pxi].err + m_proxies[pxj].err);

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
   * @brief Split one proxy area by default bisection, but N-section is also possible.
   * @pre valid facet proxy map
   * @param px_idx proxy index
   * @param n number of split sections
   * @param nb_relaxations number of relaxation on the confined proxy area
   * @return true if split succeeds, false otherwise
   */
  bool split(const std::size_t px_idx,
    const std::size_t n = 2,
    const std::size_t nb_relaxations = 10) {
    if (px_idx >= m_proxies.size() || n < 2)
      return false;

    // collect confined proxy area
    std::vector<face_descriptor> confined_area;
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
      if (get(m_fproxy_map, f) == px_idx)
        confined_area.push_back(f);
    // not enough facets to split
    if (n > confined_area.size())
      return false;

    // a copy of confined proxies
    std::vector<Proxy_wrapper> confined_proxies;
    confined_proxies.push_back(m_proxies[px_idx]);

    // select seed facets in the confiend area
    std::size_t count = 1;
    BOOST_FOREACH(face_descriptor f, confined_area) {
      if (count >= n)
        break;

      if (get(m_fproxy_map, f) == px_idx && f != m_proxies[px_idx].seed) {
        put(m_fproxy_map, f, m_proxies.size());
        m_proxies.push_back(fit_new_proxy(f, m_proxies.size()));
        ++count;
        // copy
        confined_proxies.push_back(m_proxies.back());
      }
    }

    // relaxation on confined area and proxies
    for (std::size_t i = 0; i < nb_relaxations; ++i) {
      BOOST_FOREACH(face_descriptor f, confined_area)
        put(m_fproxy_map, f, CGAL_VSA_INVALID_TAG);

      partition(confined_proxies.begin(), confined_proxies.end());
      fit(confined_proxies.begin(), confined_proxies.end());
    }

    // copy back
    BOOST_FOREACH(const Proxy_wrapper &pxw, confined_proxies)
      m_proxies[pxw.idx] = pxw;

    return true;
  }

  /*!
   * @brief Extract the approximated surface mesh.
   * @note If the extracted surface mesh contains non-manifold facets, 
   * they are not built into the output polyhedron.
   * @tparam PolyhedronSurface should be `CGAL::Polyhedron_3`
   * @param[out] tm_out output triangle mesh
   * @param chord_error boundary approximation recursively split criterion
   * @param is_relative_to_chord true if the chord_error is relative to the the chord length (relative sense),
   * otherwise it's relative to the average edge length (absolute sense).
   * @param with_dihedral_angle true if add dihedral angle weight to the distance, false otherwise
   * @param if_optimize_anchor_location true if optimize the anchor location, false otherwise
   * @param pca_plane true if use PCA plane fitting, otherwise use the default area averaged plane parameters
   * @return true if the extracted surface mesh is manifold, false otherwise.
   */
  template <typename PolyhedronSurface>
  bool extract_mesh(PolyhedronSurface &tm_out,
    const FT chord_error = FT(5.0),
    const bool is_relative_to_chord = false,
    const bool with_dihedral_angle = false,
    const bool if_optimize_anchor_location = true,
    const bool pca_plane = false) {
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
    find_edges(chord_error, is_relative_to_chord, with_dihedral_angle);
    add_anchors();

    // descrete Delaunay triangulation
    pseudo_cdt();

    if (if_optimize_anchor_location)
      optimize_anchor_location();

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
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
      facet_proxy_map[f] = get(m_fproxy_map, f);
  }

  /*!
   * @brief Get the facet region of the specified proxy.
   * @tparam OutputIterator output iterator with `boost::graph_traits<TriangleMesh>::%face_descriptor` as value type
   * @param px_idx proxy index
   * @param out_itr output iterator
   */
  template <typename OutputIterator>
  void get_proxy_region(const std::size_t px_idx, OutputIterator out_itr) const {
    if (px_idx >= m_proxies.size())
      return;

    BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
      if (get(m_fproxy_map, f) == px_idx)
        *out_itr++ = f;
  }

  /*!
   * @brief Get the proxies.
   * @tparam OutputIterator output iterator with Proxy as value type
   * @param out_itr output iterator
   */
  template <typename OutputIterator>
  void get_proxies(OutputIterator out_itr) const {
    BOOST_FOREACH(const Proxy_wrapper &pxw, m_proxies)
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
    BOOST_FOREACH(const Proxy_wrapper &pxw, m_proxies)
      *out_itr++ = pxw;
  }
#endif

  /*!
   * @brief Get the proxies size.
   * @return number of proxies
   */
  std::size_t get_proxies_size() const { return m_proxies.size(); }

  /*!
   * @brief Get the anchor points, which have the area-averaged position of the projected anchor vertex points on the incident proxies.
   * @tparam OutputIterator output iterator with Point_3 as value type
   * @param out_itr output iterator
   */
  template <typename OutputIterator>
  void get_anchor_points(OutputIterator out_itr) const {
    BOOST_FOREACH(const Anchor &a, m_anchors)
      *out_itr++ = a.pos;
  }

  /*!
   * @brief Get the anchor vertices.
   * @tparam OutputIterator output iterator with vertex_descriptor as value type
   * @param out_itr output iterator
   */
  template <typename OutputIterator>
  void get_anchor_vertices(OutputIterator out_itr) const {
    BOOST_FOREACH(const Anchor &a, m_anchors)
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
    BOOST_FOREACH(const std::vector<std::size_t> &t, m_tris)
      *out_itr++ = t;
  }

  /*!
   * @brief Get the indexed boundary polygon approximation.
   * @tparam OutputIterator output iterator with std::vector<std::size_t> as value type
   */
  template <typename OutputIterator>
  void get_indexed_boundary_polygons(OutputIterator out_itr) const {
    BOOST_FOREACH(const Boundary_cycle &bcycle, m_bcycles) {
      std::vector<std::size_t> plg;
      halfedge_descriptor he = bcycle.he_head;
      do {
        Boundary_chord chord;
        walk_to_next_anchor(he, chord);
        plg.push_back(get(m_vanchor_map, target(he, *m_ptm)));
      } while (he != bcycle.he_head);
      *out_itr++ = plg;
    }
  }

// private member functions
private:
  /*!
   * @brief Random initialize proxies to target number of proxies.
   * @note To ensure the randomness, call `std::srand()` beforehand.
   * @param max_nb_proxies maximum number of proxies, 
   * should be in range (nb_connected_components, num_faces(*m_ptm))
   * @param num_iterations number of re-fitting iterations 
   * @return number of proxies initialized
   */
  std::size_t init_random(const std::size_t max_nb_proxies,
    const std::size_t num_iterations) {
    // random shuffled facets except for the bootstrapped connected component seed facets
    std::vector<face_descriptor> shuffled_facets;
    random_shuffle_non_seed_facets(shuffled_facets);

    // reach to the number of proxies
    for (std::size_t i = 0; i < shuffled_facets.size()
      && m_proxies.size() < max_nb_proxies; ++i)
      m_proxies.push_back(fit_new_proxy(shuffled_facets[i], m_proxies.size()));
    run(num_iterations);

    return m_proxies.size();
  }

  /*!
   * @brief Incremental initialize proxies to target number of proxies.
   * @param max_nb_proxies maximum number of proxies, 
   * should be in range (nb_connected_components, num_faces(*m_ptm))
   * @param num_iterations number of re-fitting iterations 
   * before each incremental proxy insertion
   * @return number of proxies initialized
   */
  std::size_t init_incremental(const std::size_t max_nb_proxies,
    const std::size_t num_iterations) {
    if (m_proxies.size() < max_nb_proxies)
      add_proxies_furthest(max_nb_proxies - m_proxies.size(), num_iterations);

    return m_proxies.size();
  }

  /*!
   * @brief Hierarchical initialize proxies to target number of proxies.
   * @param max_nb_proxies maximum number of proxies, 
   * should be in range (nb_connected_components, num_faces(*m_ptm))
   * @param num_iterations number of re-fitting iterations
   * before each hierarchical proxy insertion
   * @return number of proxies initialized
   */
  std::size_t init_hierarchical(const std::size_t max_nb_proxies,
    const std::size_t num_iterations) {
    while (m_proxies.size() < max_nb_proxies) {
      // try to double current number of proxies each time
      std::size_t target_px = m_proxies.size();
      if (target_px * 2 > max_nb_proxies)
        target_px = max_nb_proxies;
      else
        target_px *= 2;
      // add proxies by error diffusion
      add_proxies_error_diffusion(target_px - m_proxies.size());
      run(num_iterations);
    }

    return m_proxies.size();
  }

  /*!
   * @brief Randomly initialize proxies
   * with both maximum number of proxies and minimum error drop stop criteria,
   * The first criterion met stops the seeding.
   * @note To ensure the randomness, call `std::srand()` beforehand.
   * @param max_nb_proxies maximum number of proxies, should be in range (nb_connected_components, num_faces(tm) / 3)
   * @param min_error_drop minimum error drop, should be in range (0.0, 1.0)
   * @param num_iterations number of re-fitting iterations 
   * @return number of proxies initialized
   */
  std::size_t init_random_error(const std::size_t max_nb_proxies,
    const FT min_error_drop,
    const std::size_t num_iterations) {
    // random shuffled facets except for the bootstrapped connected component seed facets
    std::vector<face_descriptor> shuffled_facets;
    random_shuffle_non_seed_facets(shuffled_facets);

    // keep a copy of the connected components seeds
    std::vector<face_descriptor> cc_seed_facets;
    BOOST_FOREACH(const Proxy_wrapper &pxw, m_proxies)
      cc_seed_facets.push_back(pxw.seed);

    const FT initial_err = compute_fitting_error();
    FT error_drop = min_error_drop * FT(2.0);
    while (m_proxies.size() < max_nb_proxies && error_drop > min_error_drop) {
      // try to double current number of proxies each time
      std::size_t target_px = m_proxies.size();
      if (target_px * 2 > max_nb_proxies)
        target_px = max_nb_proxies;
      else
        target_px *= 2;

      // reset proxies to the bootstrapped connected components
      m_proxies.clear();
      BOOST_FOREACH(face_descriptor f, cc_seed_facets)
        m_proxies.push_back(fit_new_proxy(f, m_proxies.size()));

      for (std::size_t i = 0; i < shuffled_facets.size()
        && m_proxies.size() < target_px; ++i)
        m_proxies.push_back(fit_new_proxy(shuffled_facets[i], m_proxies.size()));
      run(num_iterations);

      const FT err = compute_fitting_error();
      error_drop = err / initial_err;
    }

    return m_proxies.size();
  }

  /*!
   * @brief Incrementally initialize proxies
   * with both maximum number of proxies and minimum error drop stop criteria,
   * The first criterion met stops the seeding.
   * @param max_nb_proxies maximum number of proxies, should be in range (nb_connected_components, num_faces(tm) / 3)
   * @param min_error_drop minimum error drop, should be in range (0.0, 1.0)
   * @param num_iterations number of re-fitting iterations 
   * @return number of proxies initialized
   */
  std::size_t init_incremental_error(const std::size_t max_nb_proxies,
    const FT min_error_drop,
    const std::size_t num_iterations) {
    const FT initial_err = compute_fitting_error();
    FT error_drop = min_error_drop * FT(2.0);
    while (m_proxies.size() < max_nb_proxies && error_drop > min_error_drop) {
      add_proxy_furthest();
      run(num_iterations);
      const FT err = compute_fitting_error();
      error_drop = err / initial_err;
    }

    return m_proxies.size();
  }

  /*!
   * @brief Hierarchically initialize proxies
   * with both maximum number of proxies and minimum error drop stop criteria,
   * The first criterion met stops the seeding.
   * @param max_nb_proxies maximum number of proxies, should be in range (nb_connected_components, num_faces(tm) / 3)
   * @param min_error_drop minimum error drop, should be in range (0.0, 1.0)
   * @param num_iterations number of re-fitting iterations 
   * @return number of proxies initialized
   */
  std::size_t init_hierarchical_error(const std::size_t max_nb_proxies,
    const FT min_error_drop,
    const std::size_t num_iterations) {
    const FT initial_err = compute_fitting_error();
    FT error_drop = min_error_drop * FT(2.0);
    while (m_proxies.size() < max_nb_proxies && error_drop > min_error_drop) {
      // try to double current number of proxies each time
      std::size_t target_px = m_proxies.size();
      if (target_px * 2 > max_nb_proxies)
        target_px = max_nb_proxies;
      else
        target_px *= 2;
      add_proxies_error_diffusion(target_px - m_proxies.size());
      run(num_iterations);
      const FT err = compute_fitting_error();
      error_drop = err / initial_err;
    }

    return m_proxies.size();
  }

  /*!
   * @brief Partition the area tagged with CGAL_VSA_INVALID_TAG with proxies, global facet proxy map is updated.
   * Propagates the proxy seed facets and floods the tagged area to minimize the fitting error.
   * @tparam ProxyWrapperIterator forward iterator with Proxy_wrapper as value type
   * @param beg iterator point to the first element
   * @param end iterator point to the one past the last element
   */
  template<typename ProxyWrapperIterator>
  void partition(const ProxyWrapperIterator beg, const ProxyWrapperIterator end) {
    std::priority_queue<Facet_to_integrate> facet_pqueue;
    for (ProxyWrapperIterator pxw_itr = beg; pxw_itr != end; ++pxw_itr) {
      face_descriptor f = pxw_itr->seed;
      put(m_fproxy_map, f, pxw_itr->idx);

      BOOST_FOREACH(face_descriptor fadj, faces_around_face(halfedge(f, *m_ptm), *m_ptm)) {
        if (fadj != boost::graph_traits<TriangleMesh>::null_face()
            && get(m_fproxy_map, fadj) == CGAL_VSA_INVALID_TAG) {
          facet_pqueue.push(Facet_to_integrate(
            fadj, pxw_itr->idx, (*m_perror_metric)(fadj, pxw_itr->px)));
        }
      }
    }

    while (!facet_pqueue.empty()) {
      const Facet_to_integrate c = facet_pqueue.top();
      facet_pqueue.pop();
      if (get(m_fproxy_map, c.f) == CGAL_VSA_INVALID_TAG) {
        put(m_fproxy_map, c.f, c.px);
        BOOST_FOREACH(face_descriptor fadj, faces_around_face(halfedge(c.f, *m_ptm), *m_ptm)) {
          if (fadj != boost::graph_traits<TriangleMesh>::null_face()
            && get(m_fproxy_map, fadj) == CGAL_VSA_INVALID_TAG) {
            facet_pqueue.push(Facet_to_integrate(
              fadj, c.px, (*m_perror_metric)(fadj, m_proxies[c.px].px)));
          }
        }
      }
    }
  }

  /*!
   * @brief Refitting and update input range of proxies.
   * @tparam ProxyWrapperIterator forward iterator with Proxy_wrapper as value type
   * @param beg iterator point to the first element
   * @param end iterator point to the one past the last element
   */
  template<typename ProxyWrapperIterator>
  void fit(const ProxyWrapperIterator beg, const ProxyWrapperIterator end) {
    std::vector<std::list<face_descriptor> > px_facets(m_proxies.size());
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
      px_facets[get(m_fproxy_map, f)].push_back(f);

    // update proxy parameters and seed
    for (ProxyWrapperIterator pxw_itr = beg; pxw_itr != end; ++pxw_itr) {
      const std::size_t px_idx = pxw_itr->idx;
      *pxw_itr = fit_new_proxy(px_facets[px_idx].begin(), px_facets[px_idx].end(), px_idx);
    }
  }

  /*!
   * @brief Add a proxy seed at the facet with the maximum fitting error.
   * @pre current facet proxy map is valid, proxy error is computed
   * @note No re-fitting is performed. After the operation, the facet proxy map remains valid.
   * @return true add successfully, false otherwise
   */
  bool add_proxy_furthest() {
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
    std::cerr << "add furthest " << m_proxies.size() << std::endl;
#endif
    compute_fitting_error();
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

      FT err = (*m_perror_metric)(f, m_proxies[px_idx].px);
      if (first || max_error < err) {
        first = false;
        max_error = err;
        fworst = f;
      }
    }

    if (first)
      return false;

    put(m_fproxy_map, fworst, m_proxies.size());
    m_proxies.push_back(fit_new_proxy(fworst, m_proxies.size()));

    return true;
  }

  /*!
   * @brief Fitting a new (wrapped) proxy.
   * 1. Compute proxy parameters from a list of facets.
   * 2. Find proxy seed.
   * 3. Sum the proxy error.
   * @tparam FacetIterator face_descriptor container iterator
   * @param beg container begin
   * @param end container end
   * @param px_idx proxy index
   * @return fitted proxy wrapped with internal data
   */
  template<typename FacetIterator>
  Proxy_wrapper fit_new_proxy(const FacetIterator beg,
    const FacetIterator end,
    const std::size_t px_idx) {
    CGAL_assertion(beg != end);

    // use Proxy_fitting functor to fit proxy parameters
    const Proxy px = (*m_pproxy_fitting)(beg, end);

    // find proxy seed and sum error
    face_descriptor seed = *beg;
    FT err_min = (*m_perror_metric)(*beg, px);
    FT sum_error(0.0);
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      const FT err = (*m_perror_metric)(*fitr, px);
      sum_error += err;
      if (err < err_min) {
        err_min = err;
        seed = *fitr;
      }
    }

    return Proxy_wrapper(px, px_idx, seed, sum_error);
  }

  /*!
   * @brief Fitting a new (wrapped) proxy from a single facet.
   * 1. Compute proxy parameters from one facet.
   * 2. Find proxy seed.
   * 3. Sum the proxy error.
   * @param face_descriptor facet
   * @param px_idx proxy index
   * @return fitted proxy wrapped with internal data
   */
  Proxy_wrapper fit_new_proxy(const face_descriptor f, const std::size_t px_idx) {
    // fit proxy parameters
    std::vector<face_descriptor> fvec(1, f);
    const Proxy px = (*m_pproxy_fitting)(fvec.begin(), fvec.end());
    const FT err = (*m_perror_metric)(f, px);

    return Proxy_wrapper(px, px_idx, f, err);
  }

  /*!
   * @brief Random shuffle the non-seed facets into an empty vector,
   * to prepare for consecutive random seed facets selection (no interleaved re-fitting).
   * @param[out] facets shuffled facets vector
   */
  void random_shuffle_non_seed_facets(std::vector<face_descriptor> &facets) {
    std::set<face_descriptor> seed_facets_set;
    BOOST_FOREACH(const Proxy_wrapper &pxw, m_proxies)
      seed_facets_set.insert(pxw.seed);

    if (num_faces(*m_ptm) <= m_proxies.size())
      return;

    const std::size_t nbf = num_faces(*m_ptm) - m_proxies.size();
    facets.reserve(nbf);
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm)) {
      if (seed_facets_set.find(f) != seed_facets_set.end())
        continue;
      facets.push_back(f);
    }

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
   * @brief Initialize proxies from each connected components of the surface.
   * @note This function clears proxy vector and set facet proxy map to initial state,
   * intended only for bootstrapping initialization.
   * Coarse approximation iteration is not performed, because it's inaccurate anyway
   * and may cause serious degenerate cases(e.g. a standard cube mode).
   */
  void bootstrap_from_connected_components() {
    // set all face invalid to mark as unvisited / untagged
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
      put(m_fproxy_map, f, CGAL_VSA_INVALID_TAG);

    // prepare for connected components visiting
    std::vector<face_descriptor> cc_seed_facets;
    bool if_all_visited = false;
    std::size_t cc_idx = 0;
    face_descriptor seed_facet = *(faces(*m_ptm).first);
    while (!if_all_visited) {
      // use current seed facet to traverse the conneceted componnets
      cc_seed_facets.push_back(seed_facet);
      std::stack<face_descriptor> fstack;
      fstack.push(seed_facet);
      put(m_fproxy_map, seed_facet, cc_idx);
      while (!fstack.empty()) {
        face_descriptor active_facet = fstack.top();
        fstack.pop();
        BOOST_FOREACH(face_descriptor fadj,
          faces_around_face(halfedge(active_facet, *m_ptm), *m_ptm)) {
          if (fadj != boost::graph_traits<TriangleMesh>::null_face()
            && get(m_fproxy_map, fadj) == CGAL_VSA_INVALID_TAG) {
            fstack.push(fadj);
            put(m_fproxy_map, fadj, cc_idx);
          }
        }
      }
      // check if all visited
      if_all_visited = true;
      BOOST_FOREACH(face_descriptor f, faces(*m_ptm)) {
        if (get(m_fproxy_map, f) == CGAL_VSA_INVALID_TAG) {
          if_all_visited = false;
          ++cc_idx;
          seed_facet = f;
          break;
        }
      }
    }

    m_proxies.clear();
    BOOST_FOREACH(face_descriptor f, cc_seed_facets)
      m_proxies.push_back(fit_new_proxy(f, m_proxies.size()));
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
    std::cerr << "#cc " << m_proxies.size() << std::endl;
#endif
  }

  /*!
   * @brief Compute proxy planes.
   * The proxy may not contain the plane related properties, so we need these internal planes,
   * used in the chord subdivision and anchor location.
   * @param if_pca_plane true to use the PCA plane fitting
   */
  void compute_proxy_planes(const bool if_pca_plane) {
    // fit proxy planes, areas, normals
    std::vector<std::list<face_descriptor> > px_facets(m_proxies.size());
    BOOST_FOREACH(face_descriptor f, faces(*m_ptm))
      px_facets[get(m_fproxy_map, f)].push_back(f);

    BOOST_FOREACH(const std::list<face_descriptor> &px_patch, px_facets) {
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
   * @param chord_error boundary chord approximation recursive split creterion
   * @param is_relative_to_chord true if the chord_error is relative to the the chord length (relative sense),
   * otherwise it's relative to the average edge length (absolute sense).
   * @param with_dihedral_angle true if add dihedral angle weight to the distance, false otherwise
   */
  void find_edges(const FT chord_error,
    const bool is_relative_to_chord,
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
          chord_error, is_relative_to_chord, with_dihedral_angle);

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
        std::vector<std::size_t> t;
        t.push_back(i);
        t.push_back(j);
        t.push_back(k);
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
   * @param chord_error the chord recursive split error threshold
   * @param is_relative_to_chord true if the chord_error is relative to the the chord length (relative sense),
   * otherwise it's relative to the average edge length (absolute sense).
   * @param with_dihedral_angle true if add dihedral angle weight to the distance, false otherwise
   * @return the number of anchors of the chord apart from the first one
   */
  std::size_t subdivide_chord(
    const Boundary_chord_iterator &chord_begin,
    const Boundary_chord_iterator &chord_end,
    const FT chord_error,
    const bool is_relative_to_chord,
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
      if (is_relative_to_chord)
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

      if (criterion > chord_error)
        if_subdivide = true;
    }

    if (if_subdivide) {
      // subdivide at the most remote vertex
      attach_anchor(*chord_max);

      const std::size_t num_left = subdivide_chord(chord_begin, chord_max + 1,
        chord_error, is_relative_to_chord, with_dihedral_angle);
      const std::size_t num_right = subdivide_chord(chord_max + 1, chord_end,
        chord_error, is_relative_to_chord, with_dihedral_angle);

      return num_left + num_right;
    }

    return 1;
  }

  /*!
   * @brief Return true if the target vertex of a halfedge is attached with an anchor, and false otherwise.
   * @param he a halfedge descriptor
   */
  bool is_anchor_attached(const halfedge_descriptor &he) const {
    return is_anchor_attached(target(he, *m_ptm), m_vanchor_map);
  }

  /*!
   * @brief Check if a vertex is attached with an anchor.
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
   * @brief Attachs an anchor to the vertex.
   * @param vtx vertex
   */
  void attach_anchor(const vertex_descriptor &vtx) {
    put(m_vanchor_map, vtx, m_anchors.size());
    // default anchor location is the vertex point
    m_anchors.push_back(Anchor(vtx, m_vpoint_map[vtx]));
  }

  /*!
   * @brief Attachs an anchor to the target vertex of the halfedge.
   * @param he halfedge
   */
  void attach_anchor(const halfedge_descriptor &he) {
    attach_anchor(target(he, *m_ptm));
  }

  /*!
   * @brief Optimize the anchor location by averaging the projection points of
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
   * @brief Calculate the averaged edge length of a triangle mesh.
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
   * @brief Use an incremental builder to build and test if the indexed triangle surface is manifold
   * @tparam PolyhedronSurface should be `CGAL::Polyhedron_3`
   * @param[out] poly input polyhedorn mesh
   * @return true if build manifold surface successfully
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

} // end namespace VSA
} // end namespace CGAL

#undef CGAL_VSA_INVALID_TAG

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_APPROXIMATION_H
