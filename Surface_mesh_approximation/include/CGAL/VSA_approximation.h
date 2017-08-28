#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_APPROXIMATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_APPROXIMATION_H

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <CGAL/VSA_metrics.h>
#include <CGAL/Default.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/foreach.hpp>

#include <vector>
#include <cmath>
#include <map>
#include <queue>
#include <iterator>

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEGUB
#include <iostream>
#endif

namespace CGAL
{
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
  typename ErrorMetric_ = CGAL::Default,
  typename ProxyFitting_ = CGAL::Default,
  typename GeomTraits_ = CGAL::Default>
class VSA_approximation {
// public typedefs
public:
  // Default typdefs
  /// GeomTraits typdef
  typedef typename CGAL::Default::Get<
    GeomTraits_,
    typename Kernel_traits<
      typename boost::property_traits<VertexPointMap>::value_type
    >::Kernel >::type GeomTraits;
  /// ErrorMetric typdef
  typedef typename CGAL::Default::Get<ErrorMetric_,
    CGAL::L21Metric<TriangleMesh, VertexPointMap, GeomTraits> >::type ErrorMetric;
  /// ProxyFitting typdef
  typedef typename CGAL::Default::Get<ProxyFitting_,
    CGAL::L21ProxyFitting<TriangleMesh, VertexPointMap, GeomTraits> >::type ProxyFitting;
  /// Proxy typdef
  typedef typename ErrorMetric::Proxy Proxy;

  /// Enueration typdef
  enum Initialization {
    /// Random initialization
    RandomInit,
    /// Incremental initialization
    IncrementalInit,
    /// Hierarchical initialization
    HierarchicalInit
  };

// private typedefs and data member
private:
  // GeomTraits typedefs
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Plane_3 Plane_3;
  typedef typename GeomTraits::Construct_vector_3 Construct_vector_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename GeomTraits::Compute_scalar_product_3 Compute_scalar_product_3;

  // graph_traits typedefs
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  // internal typedefs
  typedef boost::associative_property_map<std::map<vertex_descriptor, int> > VertexAnchorMap;

  typedef std::vector<halfedge_descriptor> ChordVector;
  typedef typename ChordVector::iterator ChordVectorIterator;

  // The proxy wrapper for approximation.
  struct ProxyWrapper {
    ProxyWrapper(const Proxy &_p, const face_descriptor &_s)
      : px(_p), seed(_s), err(0) {}

    Proxy px; // parameterized proxy
    face_descriptor seed; // proxy seed
    FT err; // proxy fitting error
  };

  // The proxy fitting plane for meshing.
  struct ProxyPlane {
    ProxyPlane(const Plane_3 &_p, const Vector_3 &_n, const FT &_a)
      : plane(_p), normal(_n), area(_a) {}

    Plane_3 plane;
    Vector_3 normal;
    FT area;
  };

  // The facet candidate to be queued.
  struct FacetToIntegrate {
    FacetToIntegrate(const face_descriptor &_f, const std::size_t &_px, const FT &_err)
      : f(_f), px(_px), err(_err) {}

    bool operator<(const FacetToIntegrate &rhs) const {
      return err > rhs.err;
    }

    face_descriptor f; // facet
    std::size_t px; // proxy index
    FT err; // fitting error
  };

  // Proxy error with its index.
  struct ProxyError {
    ProxyError(const std::size_t &_px, const FT &_err)
      : px(_px), err(_err) {}

    // in ascending order
    bool operator<(const ProxyError &rhs) const {
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
  class TrianglePolyhedronBuilder : public CGAL::Modifier_base<HDS> {
  public:
    const std::vector<Point_3> &vtxs;
    const std::vector<int> &tris;
    bool is_manifold;
    TrianglePolyhedronBuilder(const std::vector<Point_3> &_vtxs,
      const std::vector<int> &_tris)
      : vtxs(_vtxs), tris(_tris), is_manifold(true) {}

    void operator()(HDS &hds) {
      CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);
      typedef typename HDS::Vertex Vertex;
      typedef typename Vertex::Point Point;
      builder.begin_surface(vtxs.size(), tris.size() / 3);
      BOOST_FOREACH(const Point_3 &v, vtxs)
        builder.add_vertex(Point(v));
      for (std::vector<int>::const_iterator itr = tris.begin(); itr != tris.end(); itr += 3) {
        if (builder.test_facet(itr, itr + 3)) {
          builder.begin_facet();
          builder.add_vertex_to_facet(*itr);
          builder.add_vertex_to_facet(*(itr + 1));
          builder.add_vertex_to_facet(*(itr + 2));
          builder.end_facet();
        }
        else {
          builder.end_surface();
          is_manifold = false;
          return;
        }
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
  const ErrorMetric *fit_error;
  // The proxy fitting functor.
  const ProxyFitting *proxy_fitting;

  Construct_vector_3 vector_functor;
  Construct_scaled_vector_3 scale_functor;
  Construct_sum_of_vectors_3 sum_functor;
  Compute_scalar_product_3 scalar_product_functor;

  // The facet proxy index map.
  std::map<face_descriptor, std::size_t> internal_fidx_map;
  boost::associative_property_map<std::map<face_descriptor, std::size_t> > fproxy_map;
  // The attached anchor index of a vertex.
  std::map<vertex_descriptor, int> internal_vidx_map;
  VertexAnchorMap vanchor_map;

  // Proxies.
  std::vector<ProxyWrapper> proxies;
  // Proxy planes
  std::vector<ProxyPlane> px_planes;

  // All anchors.
  std::vector<Anchor> anchors;
  // All borders cycles.
  std::vector<Border> borders;
  // The indexed triangle approximation.
  std::vector<int> tris;

//member functions
public:
  /*!
   * %Default constructor.
   */
  VSA_approximation() :
    m_pmesh(NULL),
    fit_error(NULL),
    proxy_fitting(NULL),
    fproxy_map(internal_fidx_map),
    vanchor_map(internal_vidx_map) {
    GeomTraits traits;
    vector_functor = traits.construct_vector_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scalar_product_functor = traits.compute_scalar_product_3_object();
  }

  /*!
   * Initialize and prepare for the approximation.
   * @param _mesh `CGAL TriangleMesh` on which approximation operate.
   * @param _point_pmap vertex point map of the mesh
   */
  VSA_approximation(const TriangleMesh &_mesh,
    const VertexPointMap &_point_pmap) :
    m_pmesh(&_mesh),
    point_pmap(_point_pmap),
    fit_error(NULL),
    proxy_fitting(NULL),
    fproxy_map(internal_fidx_map),
    vanchor_map(internal_vidx_map) {
    GeomTraits traits;
    vector_functor = traits.construct_vector_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scalar_product_functor = traits.compute_scalar_product_3_object();
  }

  /*!
   * Set the mesh for approximation and rebuild the internal data structure.
   * @pre @a _mesh.is_pure_triangle()
   * @param _mesh `CGAL TriangleMesh` on which approximation operate.
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
  void set_metric(const ErrorMetric &_error_metric,
    const ProxyFitting &_proxy_fitting) {
    fit_error = &_error_metric;
    proxy_fitting = &_proxy_fitting;
  }

  /*!
   * Rebuild the internal data structure.
   */
  void rebuild() {
    // rebuild inter data structure
    proxies.clear();
    internal_fidx_map.clear();
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      internal_fidx_map[f] = 0;

    internal_vidx_map.clear();
    BOOST_FOREACH(vertex_descriptor v, vertices(*m_pmesh))
      internal_vidx_map.insert(std::pair<vertex_descriptor, int>(v, 0));
  }

  /*!
   * @brief Seeding by targeted number of proxies.
   * @param method seeding method
   * @param num_seed target number of proxies seed
   * @param num_iterations number of iterations of coarse re-fitting
   * @return number of proxies initialized
   */
  std::size_t seeding_by_number(
    const Initialization method,
    const std::size_t num_seed,
    const std::size_t num_iterations = 5) {
    switch (method) {
      case RandomInit:
        return seed_random(num_seed);
      case IncrementalInit:
        return seed_incremental(num_seed, num_iterations);
      case HierarchicalInit:
        return seed_hierarchical(num_seed, num_iterations);
      default:
        return 0;
    }
  }

  /*!
   * @brief Seeding by targeted error drop.
   * @param method seeding method
   * @param target_drop targeted error drop to initial state, usually in range (0, 1)
   * @param num_iterations number of iterations of coarse re-fitting
   * @return number of proxies initialized
   */
  std::size_t seeding_error(
    const Initialization method,
    const FT target_drop,
    const std::size_t num_iterations) {
    switch (method) {
      case RandomInit:
        return seed_error_random(target_drop, num_iterations);
      case IncrementalInit:
        return seed_error_incremental(target_drop, num_iterations);
      case HierarchicalInit:
        return seed_error_hierarchical(target_drop, num_iterations);
      default:
        return 0;
    }
  }

  /*!
   * @brief This function run the algorithm by one step,
   * including the partitioning and fitting process.
   * @return the total fitting error of current partition to the proxies.
   */
  FT run_one_step() {
    partition();
    fit();

    return compute_fitting_error();
  }

  /*!
   * @brief This function run the algorithm until the no significant energy drop.
   * @param drop_threshold the percentage of energy drop to between two runs, usually in range [0, 1).
   * @param max_iterations the maximum number of iterations allowed
   * @return true if the algorithm converge, false otherwise.
   */
  bool run_until_convergence(const FT drop_threshold = FT(0.05),
    const std::size_t max_iterations = 100) {
    FT drop_pct(0);
    std::size_t iteration_count = 0;
    FT pre_err = compute_fitting_error();
    do {
      // average 5 steps to have smoother drop curve
      FT avg_sum_err(0);
      for (std::size_t i = 0; i < 5; ++i)
        avg_sum_err += run_one_step();
      avg_sum_err /= FT(5);
      iteration_count += 5;

      drop_pct = (pre_err - avg_sum_err) / pre_err;
      if (drop_pct < FT(0))
        drop_pct = -drop_pct;
      if (drop_pct < drop_threshold)
        return true;

      pre_err = avg_sum_err;
    } while (iteration_count < max_iterations);

    return false;
  }

  /*!
   * @brief Partition the geometry with current proxies.
   * Propagates the proxy seed facets and floods the whole mesh to minimize the fitting error.
   */
  void partition() {
#define CGAL_NOT_TAGGED_ID std::numeric_limits<std::size_t>::max()
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      fproxy_map[f] = CGAL_NOT_TAGGED_ID;

    std::priority_queue<FacetToIntegrate> facet_pqueue;
    for (std::size_t i = 0; i < proxies.size(); ++i) {
      face_descriptor f = proxies[i].seed;
      fproxy_map[f] = i;

      BOOST_FOREACH(face_descriptor fadj, faces_around_face(halfedge(f, *m_pmesh), *m_pmesh)) {
        if (fadj != boost::graph_traits<TriangleMesh>::null_face()
            && fproxy_map[fadj] == CGAL_NOT_TAGGED_ID) {
          facet_pqueue.push(FacetToIntegrate(
            fadj, i, (*fit_error)(fadj, proxies[i].px)));
        }
      }
    }

    while (!facet_pqueue.empty()) {
      const FacetToIntegrate c = facet_pqueue.top();
      facet_pqueue.pop();
      if (fproxy_map[c.f] == CGAL_NOT_TAGGED_ID) {
        fproxy_map[c.f] = c.px;
        BOOST_FOREACH(face_descriptor fadj, faces_around_face(halfedge(c.f, *m_pmesh), *m_pmesh)) {
          if (fadj != boost::graph_traits<TriangleMesh>::null_face()
            && fproxy_map[fadj] == CGAL_NOT_TAGGED_ID) {
            facet_pqueue.push(FacetToIntegrate(
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
   * @brief Adding proxies. The proxies are not updated via fitting process.
   * @param adding_method select one of the adding method: hierarchical or incremental(furthest).
   * @param num_proxies number of proxies
   * @param inner_iteration coarse re-fitting before each insertion if incremental is selected
   * @return number of proxies successfully added.
   */
  std::size_t add_proxies(const Initialization &adding_method,
    const std::size_t &num_proxies = 1,
    const std::size_t inner_iteration = 5) {
    switch (adding_method) {
      case HierarchicalInit:
        return insert_proxy_hierarchical(num_proxies);
      case IncrementalInit:
        return insert_proxy_furthest(num_proxies, inner_iteration);
      default:
        return 0;
    }
  }

  /*!
   * @brief Teleport the local minima to the worst region, this combines the merging and adding processes.
   * The partitioning are updated.
   * Here if we specify more than one proxy this means we teleport in a naive iterative fashion.
   * @param num_proxies number of proxies request to teleport
   * @param if_test true if do the merge test before the teleportation (attempt to escape from local minima).
   * @return number of proxies teleported.
   */
  std::size_t teleport_proxies(const std::size_t &num_proxies, const bool &if_test = true) {
    std::size_t num_teleported = 0;
    while (num_teleported < num_proxies) {
      // find worst proxy
      std::vector<FT> px_error(proxies.size(), FT(0));
      compute_fitting_error(px_error);
      std::size_t px_worst = 0;
      FT max_error = px_error.front();
      for (std::size_t i = 0; i < proxies.size(); ++i) {
        if (max_error < px_error[i]) {
          max_error = px_error[i];
          px_worst = i;
        }
      }
      bool found = false;
      face_descriptor tele_to;
      BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
        if (found)
          break;
        if (fproxy_map[f] == px_worst) {
          if (f != proxies[px_worst].seed) {
            tele_to = f;
            found = true;
          }
        }
      }
      // no where to teleport
      if (!found)
        return num_teleported;

      // find the best merge pair
      std::size_t px_enlarged = 0, px_merged = 0;
      if (!find_best_merge(px_enlarged, px_merged, if_test))
        return num_teleported;
      if (px_worst == px_enlarged || px_worst == px_merged)
        return num_teleported;

      // teleport to a facet to the worst region
      fproxy_map[tele_to] = proxies.size();
      proxies.push_back(fit_new_proxy(tele_to));

      merge(px_enlarged, px_merged);
      num_teleported++;

      // coarse re-fitting
      for (std::size_t i = 0; i < 5; ++i) {
        partition();
        fit();
      }

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEGUB
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
      return FT(0);

    // ensure px0 < px1
    if (px0 > px1)
      std::swap(px0, px1);

    // merge px1 to px0
    FT err_sum(0);
    std::list<face_descriptor> merged_patch;
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
      std::size_t px_idx = fproxy_map[f];
      if (px_idx == px1) {
        err_sum += (*fit_error)(f, proxies[px_idx].px);
        fproxy_map[f] = px0;
        merged_patch.push_back(f);
      }
      else if (px_idx == px0) {
        err_sum += (*fit_error)(f, proxies[px_idx].px);
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

    FT err_merged(0);
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
    // Proxy merged_px;
    FT min_merged_error = FT(0);
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

      ProxyWrapper px = fit_new_proxy(merged_patch.begin(), merged_patch.end());
      FT sum_error(0);
      BOOST_FOREACH(face_descriptor f, merged_patch)
        sum_error += (*fit_error)(f, px.px);
      merged_set.insert(ProxyPair(pxi, pxj));

      if (first_merge || sum_error < min_merged_error) {
        first_merge = false;
        min_merged_error = sum_error;
        // merged_px = px;
        px_enlarged = pxi;
        px_merged = pxj;
      }
    }

    std::vector<FT> px_error(proxies.size(), FT(0));
    compute_fitting_error(px_error);
    FT max_error = px_error.front();
    for (std::size_t i = 0; i < proxies.size(); ++i) {
      if (max_error < px_error[i])
        max_error = px_error[i];
    }

    // test if merge worth it
    if (if_test) {
      const FT merge_thre = max_error / FT(2);
      const FT increase = min_merged_error - (px_error[px_enlarged] + px_error[px_merged]);
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
      return FT(0);

    std::size_t count = 1;
    FT err(0);
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
      if (count >= n)
        break;

      if (fproxy_map[f] == px && f != proxies[px].seed) {
        err += (*fit_error)(f, proxies[px].px);
        fproxy_map[f] = proxies.size();
        proxies.push_back(fit_new_proxy(f));
        ++count;
      }
    }

    return err;
  }

  /*!
   * @brief Meshing, choose the default area weighted or the PCA plane fitting.
   * @tparam PolyhedronSurface should be `CGAL::Polyhedron_3`
   * @param[out] tm_out output triangle mesh
   * @param split_criterion boundary approximation recursively split criterion
   * @param pca_plane if use PCA plane fitting method
   * @return true if output triangle mesh is manifold,false otherwise.
   */
  template <typename PolyhedronSurface>
  bool meshing(PolyhedronSurface &tm_out, const FT split_criterion = FT(0.2), bool pca_plane = false) {
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
    find_edges(split_criterion);
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
  void get_proxy_map(FacetProxyMap &facet_proxy_map) {
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      facet_proxy_map[f] = fproxy_map[f];
  }

  /*!
   * @brief Get the proxies.
   * @tparam OutputIterator output iterator with Proxy as value type
   * @param out_itr output iterator
   */
  template <typename OutputIterator>
  void get_proxies(OutputIterator out_itr) {
    BOOST_FOREACH(const ProxyWrapper &pxw, proxies) {
      *out_itr = pxw.px;
      ++out_itr;
    }
  }

  /*!
   * @brief Get the proxies size.
   * @return number of proxies
   */
  std::size_t get_proxies_size() { return proxies.size(); }

  /*!
   * @brief Get the anchor points, which have the area-averaged position of the projected anchor vertex points on the incident proxies.
   * @tparam OutputIterator output iterator with Point_3 as value type
   * @param out_itr output iterator
   */
  template <typename OutputIterator>
  void get_anchor_points(OutputIterator out_itr) {
    BOOST_FOREACH(const Anchor &a, anchors) {
      *out_itr = a.pos;
      ++out_itr;
    }
  }

  /*!
   * @brief Get the anchor vertices.
   * @tparam OutputIterator output iterator with vertex_descriptor as value type
   * @param out_itr output iterator
   */
  template <typename OutputIterator>
  void get_anchor_vertices(OutputIterator out_itr) {
    BOOST_FOREACH(const Anchor &a, anchors) {
      *out_itr = a.vtx;
      ++out_itr;
    }
  }

  /*!
   * @brief Get the indexed triangles, one triplet of integers per triangles, and that the integers refer to the anchor point indexes.
   * @tparam OutputIterator output iterator with std::size_t as value type
   * @param out_itr output iterator
   */
  template <typename OutputIterator>
  void get_indexed_triangles(OutputIterator out_itr) {
    BOOST_FOREACH(const int &i, tris) {
      *out_itr = i;
      ++out_itr;
    }
  }

  /*!
   * @brief Get the indexed boundary polygon approximation.
   * @return vector of indexed polygons.
   */
  std::vector<std::vector<std::size_t> > get_indexed_boundary_polygons() {
    std::vector<std::vector<std::size_t> > bdrs;
    for (typename std::vector<Border>::iterator bitr = borders.begin();
      bitr != borders.end(); ++bitr) {
      std::vector<std::size_t> bdr;
      const halfedge_descriptor he_mark = bitr->he_head;
      halfedge_descriptor he = he_mark;
      do {
        ChordVector chord;
        walk_to_next_anchor(he, chord);
        bdr.push_back(vanchor_map[target(he, *m_pmesh)]);
      } while(he != he_mark);
      bdrs.push_back(bdr);
    }
    return bdrs;
  }

// private member functions
private:
  /*!
   * @brief Random initialize proxies.
   * @param num_seed number of proxies seed
   * @return number of proxies initialized
   */
  std::size_t seed_random(const std::size_t num_seed) {
    proxies.clear();
    if (num_faces(*m_pmesh) < num_seed)
      return 0;

    const std::size_t interval = num_faces(*m_pmesh) / num_seed;
    std::size_t index = 0;
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
      if ((index++) % interval == 0) {
        proxies.push_back(fit_new_proxy(f));
      }
      if (proxies.size() >= num_seed)
        break;
    }
    return proxies.size();
  }

  /*!
   * @brief Incremental initialize proxies.
   * @param num_seed number of proxies seed
   * @param num_iterations number of iterations of coarse re-fitting
   * before each incremental proxy insertion
   * @return number of proxies initialized
   */
  std::size_t seed_incremental(const std::size_t num_seed,
    const std::size_t num_iterations = 5) {
    proxies.clear();
    if (num_faces(*m_pmesh) < num_seed)
      return 0;
    
    // initialize a proxy and the proxy map to prepare for the insertion
    proxies.push_back(fit_new_proxy(*(faces(*m_pmesh).first)));
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      fproxy_map[f] = 0;

    insert_proxy_furthest(num_seed - 1, num_iterations);
    return proxies.size();
  }

  /*!
   * @brief Hierarchical initialize proxies.
   * @param num_seed number of proxies seed
   * @param num_iterations number of iterations of coarse re-fitting
   * before each hierarchical proxy insertion
   * @return number of proxies initialized
   */
  std::size_t seed_hierarchical(const std::size_t num_seed,
    const std::size_t num_iterations = 5) {
    proxies.clear();
    if (num_faces(*m_pmesh) < num_seed)
      return 0;
    
    // initialize 2 proxy
    typename boost::graph_traits<TriangleMesh>::face_iterator
      fitr = faces(*m_pmesh).first;
    proxies.push_back(fit_new_proxy(*fitr));
    proxies.push_back(fit_new_proxy(*(++fitr)));

    while (proxies.size() < num_seed) {
      for (std::size_t i = 0; i < num_iterations; ++i) {
        partition();
        fit();
      }

      // add proxies by error diffusion
      const std::size_t num_proxies = proxies.size();
      const std::size_t num_proxies_to_be_added =
        (num_proxies * 2 < num_seed) ? num_proxies : (num_seed - num_proxies);
      insert_proxy_hierarchical(num_proxies_to_be_added);
    }
    return proxies.size();
  }
  
  /*!
   * @brief Initialize by targeted error drop.
   * @param target_drop targeted error drop to initial state, usually in range (0, 1)
   * @param num_iterations number of iterations of coarse re-fitting
   * @return number of proxies initialized
   */
  std::size_t seed_error_random(const FT target_drop,
    const std::size_t num_iterations) {
    // maximum allowed number of proxies
    const std::size_t max_proxies = num_faces(*m_pmesh) / 3;
    if (max_proxies < 1)
      return 0;

    // initialize a proxy and the proxy map to prepare for the insertion
    proxies.clear();
    proxies.push_back(fit_new_proxy(*(faces(*m_pmesh).first)));
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      fproxy_map[f] = 0;
    const FT initial_err = compute_fitting_error();

    FT sum_err(0);
    FT drop(0);
    std::size_t target_px = 2;
    do {
      proxies.clear();
      seed_random(target_px);
      for (std::size_t i = 0; i < num_iterations; ++i) {
        partition();
        fit();
      }
      sum_err = compute_fitting_error();
      target_px *= 2;
      drop = sum_err / initial_err;
    } while(drop > target_drop && proxies.size() < max_proxies);

    return proxies.size();
  }

  /*!
   * @brief Initialize by targeted error drop.
   * @param target_drop targeted error drop to initial state, usually in range (0, 1)
   * @param num_iterations number of iterations of coarse re-fitting
   * @return number of proxies initialized
   */
  std::size_t seed_error_incremental(const FT target_drop,
    const std::size_t num_iterations) {
    // maximum allowed number of proxies
    const std::size_t max_proxies = num_faces(*m_pmesh) / 3;
    if (max_proxies < 1)
      return 0;

    // initialize a proxy and the proxy map to prepare for the insertion
    proxies.clear();
    proxies.push_back(fit_new_proxy(*(faces(*m_pmesh).first)));
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      fproxy_map[f] = 0;
    const FT initial_err = compute_fitting_error();

    FT sum_err(0);
    FT drop(0);
    do {
      insert_proxy_furthest();
      for (std::size_t i = 0; i < num_iterations; ++i) {
        partition();
        fit();
      }
      sum_err = compute_fitting_error();
      drop = sum_err / initial_err;
    } while (drop > target_drop && proxies.size() < max_proxies);

    return proxies.size();
  }

  /*!
   * @brief Initialize by targeted error drop.
   * @param target_drop targeted error drop to initial state, usually in range (0, 1)
   * @param num_iterations number of iterations of coarse re-fitting
   * @return number of proxies initialized
   */
  std::size_t seed_error_hierarchical(const FT target_drop,
    const std::size_t num_iterations) {
    // maximum allowed number of proxies
    const std::size_t max_proxies = num_faces(*m_pmesh) / 3;
    if (max_proxies < 1)
      return 0;

    // initialize a proxy and the proxy map to prepare for the insertion
    proxies.clear();
    proxies.push_back(fit_new_proxy(*(faces(*m_pmesh).first)));
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      fproxy_map[f] = 0;
    const FT initial_err = compute_fitting_error();

    FT sum_err(0);
    FT drop(0);
    std::size_t target_px = 1;
    do {
      insert_proxy_hierarchical(target_px);
      for (std::size_t i = 0; i < num_iterations; ++i) {
        partition();
        fit();
      }
      sum_err = compute_fitting_error();
      target_px *= 2;
      drop = sum_err / initial_err;
    } while(drop > target_drop && proxies.size() < max_proxies);

    return proxies.size();
  }

  /*!
   * @brief Inserts a proxy at the furthest facet of the region with the maximum fitting error.
   * No re-fitting is performed.
   * @return true if insertion success, false otherwise
   */
  bool insert_proxy_furthest() {
    std::vector<FT> px_error(proxies.size(), FT(0.0));
    std::vector<FT> max_facet_error(proxies.size(), FT(0.0));
    std::vector<face_descriptor> max_facet(proxies.size());

    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
      std::size_t px_idx = fproxy_map[f];
      FT err = (*fit_error)(f, proxies[px_idx].px);
      px_error[px_idx] += err;

      if (err > max_facet_error[px_idx]) {
        max_facet_error[px_idx] = err;
        max_facet[px_idx] = f;
      }
    }

    FT max_px_error = px_error.front();
    std::size_t max_px_idx = 0;
    for (std::size_t i = 0; i < proxies.size(); ++i) {
      if (px_error[i] > max_px_error) {
        max_px_error = px_error[i];
        max_px_idx = i;
      }
    }
    if (max_facet[max_px_idx] == proxies[max_px_idx].seed)
      return false;

    proxies.push_back(fit_new_proxy(max_facet[max_px_idx]));
    return true;
  }

  /*!
   * @brief Inserts more than one proxies to the regions with the maximum fitting error.
   * Except for the first one, a coarse re-fitting is performed before each proxy is inserted.
   * @param num_proxies number of proxies to be inserted
   * @param inner_iteration the number of iterations of coarse re-fitting
   * @return number of proxies inserted
   */
  std::size_t insert_proxy_furthest(const std::size_t num_proxies,
    const std::size_t inner_iteration = 5) {
    // when insert only one proxy, it has the same effect of insert_proxy_furthest()
    if (num_proxies == 0 || !insert_proxy_furthest())
      return 0;

    std::size_t num_inserted = 1;
    for (; num_inserted < num_proxies; ++num_inserted) {
      for (std::size_t i = 0; i < inner_iteration; ++i) {
        partition();
        fit();
      }
      if (!insert_proxy_furthest())
        return num_inserted;
    }
    return num_inserted;
  }

  /*!
   * @brief Add proxies by diffusing fitting error into current partitions.
   * Each partition is added with the number of proxies in proportional to its fitting error.
   * Note that the number of inserted proxies doesn't necessarily equal the requested number.
   * @param num_proxies_to_be_added added number of proxies
   * @return inserted number of proxies
   */
  std::size_t insert_proxy_hierarchical(const std::size_t num_proxies_to_be_added) {

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEGUB
    std::cerr << "#px " << proxies.size() << std::endl;
#endif

    std::vector<FT> err(proxies.size(), FT(0));
    const FT sum_error = compute_fitting_error(err);
    const FT avg_error = sum_error / FT(static_cast<double>(num_proxies_to_be_added));

    std::vector<ProxyError> px_error;
    for (std::size_t i = 0; i < proxies.size(); ++i)
      px_error.push_back(ProxyError(i, err[i]));
    // sort partition by error
    std::sort(px_error.begin(), px_error.end());

    // number of proxies to be added to each region
    std::vector<std::size_t> num_to_add(proxies.size(), 0);
    // residual from previous proxy in range (-0.5, 0.5] * avg_error
    FT residual(0);
    BOOST_FOREACH(const ProxyError &pxe, px_error) {
      // add error residual from previous proxy
      // to_add maybe negative but greater than -0.5
      FT to_add = (residual + pxe.err) / avg_error;
      // floor_to_add maybe negative but no less than -1
      FT floor_to_add = FT(std::floor(CGAL::to_double(to_add)));
      const std::size_t q_to_add = static_cast<std::size_t>(CGAL::to_double(
        ((to_add - floor_to_add) > FT(0.5)) ? (floor_to_add + FT(1)) : floor_to_add));
      residual = (to_add - FT(static_cast<double>(q_to_add))) * avg_error;
      num_to_add[pxe.px] = q_to_add;
    }

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEGUB
    for (std::size_t i = 0; i < px_error.size(); ++i)
      std::cerr << "#px " << px_error[i].px
        << ", #error " << px_error[i].err
        << ", #num_to_add " << num_to_add[px_error[i].px] << std::endl;
#endif

    std::size_t num_inserted = 0;
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
      const std::size_t px_id = fproxy_map[f];
      if (proxies[px_id].seed == f)
        continue;

      if (num_to_add[px_id] > 0) {
        proxies.push_back(fit_new_proxy(f));
        --num_to_add[px_id];
        ++num_inserted;
      }
    }

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEGUB
    std::cerr << "#requested/inserted "
      << num_proxies_to_be_added << '/' << num_inserted << std::endl;
#endif

    return num_inserted;
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
  ProxyWrapper fit_new_proxy(const FacetIterator &beg, const FacetIterator &end) {
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

    return ProxyWrapper(px, seed);
  }

  /*!
   * @brief Fitting a new proxy from a single facet.
   * 1. Fit proxy parameters from one facet.
   * 2. Set seed.
   * @param face_descriptor facet
   */
  ProxyWrapper fit_new_proxy(const face_descriptor &f) {
    std::vector<face_descriptor> fvec(1, f);
    // fit proxy parameters
    Proxy px = (*proxy_fitting)(fvec.begin(), fvec.end());

    return ProxyWrapper(px, f);
  }

  /*!
   * @brief Computes fitting error of a current partition and proxies.
   * @return total fitting error
   */
  FT compute_fitting_error() {
    FT sum_error(0);
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
      sum_error += (*fit_error)(f, proxies[fproxy_map[f]].px);
    return sum_error;
  }

  /*!
   * @brief Computes fitting error of a current partition and proxies.
   * @param px_error vector of error of each proxy
   * @return total fitting error
   */
  FT compute_fitting_error(std::vector<FT> &px_error) {
    FT sum_error(0);
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
      const std::size_t px_idx = fproxy_map[f];
      FT err = (*fit_error)(f, proxies[px_idx].px);
      px_error[px_idx] += err;
      sum_error += err;
    }
    return sum_error;
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
      FT area(0);
      BOOST_FOREACH(face_descriptor f, px_patch) {
        halfedge_descriptor he = halfedge(f, *m_pmesh);
        const Point_3 &p0 = point_pmap[source(he, *m_pmesh)];
        const Point_3 &p1 = point_pmap[target(he, *m_pmesh)];
        const Point_3 &p2 = point_pmap[target(next(he, *m_pmesh), *m_pmesh)];
        FT farea(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
        Vector_3 fnorm = CGAL::unit_normal(p0, p1, p2);

        norm = sum_functor(norm, scale_functor(fnorm, farea));
        area += farea;
      }
      norm = scale_functor(norm, FT(1.0 / std::sqrt(CGAL::to_double(norm.squared_length()))));

      px_planes.push_back(ProxyPlane(fit_plane, norm, area));
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
   * @brief Finds and approximates the edges connecting the anchors.
   * @param split_criterion edge approximation recursive split creterion
   */
  void find_edges(const FT split_criterion) {
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

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEGUB
      std::cerr << "#border " << borders.size() << std::endl;
#endif

      const halfedge_descriptor he_mark = he_start;
      do {
        ChordVector chord;
        walk_to_next_anchor(he_start, chord);
        borders.back().num_anchors += subdivide_chord(chord.begin(), chord.end(), split_criterion);

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEGUB
        std::cerr << "#chord_anchor " << borders.back().num_anchors << std::endl;
#endif

        for (ChordVectorIterator citr = chord.begin(); citr != chord.end(); ++citr)
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
      ChordVector chord;
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
      } while(he != he_mark);

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
      for (ChordVectorIterator citr = chord.begin(); citr != chord.end(); ++citr) {
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

    typedef std::map<vertex_descriptor, sg_vertex_descriptor> VertexMap;
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
          add_edge(superv, v, FT(0), gmain);
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
        ChordVector chord;
        walk_to_next_anchor(he, chord);

        std::vector<FT> vdist;
        vdist.push_back(FT(0));
        BOOST_FOREACH(halfedge_descriptor h, chord) {
          FT elen = global_eweight_map[edge(
            to_sgv_map[source(h, *m_pmesh)],
            to_sgv_map[target(h, *m_pmesh)],
            gmain).first];
          vdist.push_back(vdist.back() + elen);
        }

        FT half_chord_len = vdist.back() / FT(2);
        const int anchorleft = vanchor_map[source(chord.front(), *m_pmesh)];
        const int anchorright = vanchor_map[target(chord.back(), *m_pmesh)];
        typename std::vector<FT>::iterator ditr = vdist.begin() + 1;
        for (typename ChordVector::iterator hitr = chord.begin();
          hitr != chord.end() - 1; ++hitr, ++ditr) {
          if (*ditr < half_chord_len)
            global_vtag_map[to_sgv_map[target(*hitr, *m_pmesh)]] = anchorleft;
          else
            global_vtag_map[to_sgv_map[target(*hitr, *m_pmesh)]] = anchorright;
        }
      } while(he != he_mark);
    }

    // collect triangles
    BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
      halfedge_descriptor he = halfedge(f, *m_pmesh);
      int i = global_vtag_map[to_sgv_map[source(he, *m_pmesh)]];
      int j = global_vtag_map[to_sgv_map[target(he, *m_pmesh)]];
      int k = global_vtag_map[to_sgv_map[target(next(he, *m_pmesh), *m_pmesh)]];
      if (i != j && i != k && j != k) {
        tris.push_back(i);
        tris.push_back(j);
        tris.push_back(k);
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
  void walk_to_next_anchor(halfedge_descriptor &he_start, ChordVector &chord) {
    do {
      walk_to_next_border_halfedge(he_start);
      chord.push_back(he_start);
    } while (!is_anchor_attached(he_start));
  }

  /*!
   * @brief Walks to next border halfedge.
   * @param[in/out] he_start region border halfedge
   */
  void walk_to_next_border_halfedge(halfedge_descriptor &he_start) {
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
   * @param thre the recursive split threshold
   * @return the number of anchors of the chord apart from the first one
   */
  std::size_t subdivide_chord(
    const ChordVectorIterator &chord_begin,
    const ChordVectorIterator &chord_end,
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
    ChordVectorIterator chord_max;
    const Point_3 &pt_begin = point_pmap[source(he_first, *m_pmesh)];
    const Point_3 &pt_end = point_pmap[target(he_last, *m_pmesh)];
    if (anchor_first == anchor_last) {
      // circular chord
      CGAL_assertion(chord_size > 2);

      FT dist_max(0.0);
      for (ChordVectorIterator citr = chord_begin; citr != chord_end; ++citr) {
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

      for (ChordVectorIterator citr = chord_begin; citr != chord_end; ++citr) {
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
   * @brief Check if the target vertex of a halfedge is attached with an anchor.
   * @param he halfedge
   */
  bool is_anchor_attached(const halfedge_descriptor &he) {
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
    const VertexAnchorIndexMap &vertex_anchor_map) {
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
    FT avgx(0), avgy(0), avgz(0), sum_area(0);
    const Point_3 vtx_pt = point_pmap[v];
    for (std::set<std::size_t>::iterator pxitr = px_set.begin();
      pxitr != px_set.end(); ++pxitr) {
      std::size_t px_idx = *pxitr;
      Point_3 proj = px_planes[px_idx].plane.projection(vtx_pt);
      FT area = px_planes[px_idx].area;
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
    TrianglePolyhedronBuilder<HDS> tpbuilder(vtx, tris);
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
    FT sum_area(0);
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
    cent = scale_functor(cent, FT(1) / sum_area);

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

    typedef typename GeomTraits::Triangle_3 Triangle_3;
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

} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_APPROXIMATION_H
