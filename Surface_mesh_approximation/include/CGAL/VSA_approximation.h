#ifndef CGAL_VSA_APPROXIMATION
#define CGAL_VSA_APPROXIMATION

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/foreach.hpp>

#include <vector>
#include <cmath>
#include <map>
#include <queue>
#include <iostream>
#include <iterator>

#define CGAL_NOT_TAGGED_ID std::numeric_limits<std::size_t>::max()

namespace CGAL
{
/*!
 * @brief Main class for Variational Shape Approximation algorithm.
 * @tparam TriangleMesh a CGAL TriangleMesh
 * @tparam ErrorMetric error metric type
 * @tparam ProxyFitting proxy fitting type
 * @tparam GeomTraits geometric traits type
 */
template <typename TriangleMesh,
  typename Proxy,
  typename ErrorMetric,
  typename ProxyFitting,
  typename GeomTraits = typename TriangleMesh::Traits>
class VSA_approximation {
  // type definitions
private:
  typedef typename GeomTraits::FT FT;
  // typedef typename ErrorMetric::Proxy Proxy;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  // The facet candidate to be queued.
  struct FacetToIntegrate {
    face_descriptor f;
    std::size_t i;
    FT fit_error;
    bool operator<(const FacetToIntegrate &rhs) const {
      return fit_error > rhs.fit_error;
    }
  };

  // proxy error with proxy index
  struct ProxyError {
    ProxyError(const std::size_t &id, const FT &er)
      : px_idx(id), fit_error(er) {}
    // in ascending order
    bool operator<(const ProxyError &rhs) const {
      return fit_error < rhs.fit_error;
    }
    std::size_t px_idx;
    FT fit_error;
  };

public:
  enum Initialization {
    RandomInit,
    IncrementalInit,
    HierarchicalInit
  };

  // member variables
private:
  // TODO, update mesh
  const TriangleMesh &mesh;

  // Proxy.
  std::vector<Proxy> proxies;

  // facet proxy index map
  std::map<face_descriptor, std::size_t> internal_fidx_map;
  boost::associative_property_map<std::map<face_descriptor, std::size_t> > seg_pmap;

  // The error metric.
  const ErrorMetric &fit_error;

  // The proxy fitting functor.
  const ProxyFitting &proxy_fitting;
  //member functions
public:
  /*!
   * Initialize and prepare for the approximation.
   * @param _fit_error error calculation functor.
   * @param _proxy_fitting proxy fitting functor.
   */
  VSA_approximation(const ErrorMetric &_fit_error,
    const ProxyFitting &_proxy_fitting) :
    fit_error(_fit_error),
    proxy_fitting(_proxy_fitting),
    seg_pmap(internal_fidx_map) {}

  /*!
   * Set the mesh for approximation and rebuild the internal data structure.
   * @pre @a _mesh.is_pure_triangle()
   * @param _mesh `CGAL TriangleMesh` on which approximation operate.
   */
  void set_mesh(const TriangleMesh &_mesh) {
    mesh = _mesh;
    rebuild();
  }

  /*!
   * Rebuild the internal data structure.
   */
  void rebuild() {
    proxies.clear();
    internal_fidx_map.clear();
    BOOST_FOREACH(face_descriptor f, faces(mesh))
      internal_fidx_map[f] = 0;
  }

  /*!
   * @brief Initialize by number of proxies.
   * @param num_proxy number of proxies
   * @param seeding_method select one of the seeding method: random, hierarchical, incremental
   * @return #proxies initialized
   */
  std::size_t init_proxies(const std::size_t num_proxy, const Initialization &seeding_method) {
    proxies.clear();
    if (num_faces(mesh) < num_proxy)
      return 0;

    switch (seeding_method) {
      case IncrementalInit:
        return seed_incremental(num_proxy);
      case HierarchicalInit:
        return seed_hierarchical(num_proxy);
      default:
        return seed_random(num_proxy);
    }
  }

  /*!
   * @brief Initialize by targeted error.
   * @param target_error targeted error
   * @param seeding_method select one of the seeding method: random, hierarchical, incremental
   */
  void init_error(const FT &target_error, const Initialization &seeding_method) {
    // TODO
  }

  /*!
   * @brief This function run the algorithm by one step,
   * including the partitioning and fitting process.
   * @return true if partitioning has changed, false otherwise.
   */
  bool run_one_step() {
    partition();
    fit();

    // TODO
    return false;
  }

  /*!
   * @brief This function run the algorithm until the stop criterion is met.
   * @return true if the algorithm converge, false otherwise.
   */
  bool run_until_convergence(stop_criterion) {
    while (stop_criterion) {
      run_one_step();
    }

    // TODO
    return false;
  }

  /*!
   * @brief Partition the geometry with current proxies.
   * Propagates the proxy seed facets and floods the whole mesh to minimize the fitting error.
   */
  void partition() {
    BOOST_FOREACH(face_descriptor f, faces(mesh))
      seg_pmap[f] = CGAL_NOT_TAGGED_ID;

    std::priority_queue<FacetToIntegrate> facet_pqueue;
    for (std::size_t i = 0; i < proxies.size(); ++i) {
      face_descriptor f = proxies[i].seed;
      seg_pmap[f] = i;

      BOOST_FOREACH(face_descriptor fadj, faces_around_face(halfedge(f, mesh), mesh)) {
        if (fadj != boost::graph_traits<TriangleMesh>::null_face()
            && seg_pmap[fadj] == CGAL_NOT_TAGGED_ID) {
          FacetToIntegrate cand;
          cand.f = fadj;
          cand.fit_error = fit_error(fadj, proxies[i]);
          cand.i = i;
          facet_pqueue.push(cand);
        }
      }
    }

    while (!facet_pqueue.empty()) {
      const FacetToIntegrate c = facet_pqueue.top();
      facet_pqueue.pop();
      if (seg_pmap[c.f] == CGAL_NOT_TAGGED_ID) {
        seg_pmap[c.f] = c.i;
        BOOST_FOREACH(face_descriptor fadj, faces_around_face(halfedge(c.f, mesh), mesh)) {
          if (fadj != boost::graph_traits<TriangleMesh>::null_face()
            && seg_pmap[fadj] == CGAL_NOT_TAGGED_ID) {
            FacetToIntegrate cand;
            cand.f = fadj;
            cand.fit_error = fit_error(fadj, proxies[c.i]);
            cand.i = c.i;
            facet_pqueue.push(cand);
          }
        }
      }
    }
  }

  /*!
   * @brief Refitting of current partitioning, update proxy parameters.
   * Calculates and updates the best fitting proxies of current partition.
   * @return the total fitting error.
   */
  FT fit() {
    std::vector<std::list<face_descriptor> > px_facets(proxies.size());
    BOOST_FOREACH(face_descriptor f, faces(mesh))
      px_facets[seg_pmap[f]].push_back(f);

    // update proxy parameters and seed
    for (std::size_t i = 0; i < proxies.size(); ++i)
      proxies[i] = fit_new_proxy(px_facets[i].begin(), px_facets[i].end());
  }

  /*!
   * @brief Adding proxies. The proxies are not updated via fitting process.
   * @param num_proxies number of proxies
   * @param adding_method select one of the adding method: hierarchical or incremental(furthest).
   * @return #proxies successfully added.
   */
  std::size_t add_proxies(const Initialization &adding_method, const std::size_t &num_proxies = 1) {
    switch (adding_method) {
      case HierarchicalInit:
        return insert_proxy_hierarchical(num_proxies);
      case IncrementalInit:
        return insert_proxy_furthest(num_proxies);
      default:
        return 0;
    }
  }

  /*!
   * @brief Teleport the local minima to other places, this combines the merging and adding processes
   * The proxies are not updated. Here if we specify more than one proxy this means we
   * need to select on one side the set of N best pairs of adjacent regions to merge,
   * and match them with N regions in most need of being split. The matching can be done in a naive iterative fashion for now.
   * @param force true if force the teleportation (attempt to escape from local minima), else teleport only if it is beneficial to the error.
   * @return #proxies teleported.
   */
  std::size_t teleport_proxies(const std::size_t &num_proxies, const bool &force) {
    typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
    typedef std::pair<std::size_t, std::size_t> ProxyPair;
    typedef std::set<ProxyPair> MergedPair;

    // find worst proxy
    std::vector<FT> px_error(proxies.size(), FT(0));
    fitting_error(seg_pmap, px_error);
    std::size_t px_worst = 0;
    FT max_error = px_error.front();
    for (std::size_t i = 0; i < proxies.size(); ++i) {
      if (max_error < px_error[i]) {
        max_error = px_error[i];
        px_worst = i;
      }
    }
    std::list<face_descriptor> worst_patch;
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      if (seg_pmap[px_worst] == px_worst)
        worst_patch.push_back(f);
    }
    // not enough facets to teleport
    if (worst_patch.size() < num_proxies + 1)
      return 0;

    std::size_t num_merged = 0;
    while (num_merged < num_proxies) {
      // find the best merge pair
      std::size_t px_enlarged = 0, px_merged = 0;
      if (!find_best_merge(px_enlarged, px_merged, !force))
        return num_merged;
      if (px_worst == px_enlarged || px_worst == px_merged)
        return num_merged;

      // teleport to a facet in the worst region
      BOOST_FOREACH(face_descriptor f, worst_patch) {
        if (f != proxies[px_worst].seed) {
          proxies.push_back(fit_new_proxy(f));
          worst_patch.remove(f);
          break;
        }
      }

      if (merge(px_enlarged, px_merged))
        num_merged++;
      else
        return num_merged;

      std::cerr << "teleported" << std::endl;
    }

    return num_merged;
  }

  /*!
   * @brief Merge two specified adjacent regions, the re-fitting is performed.
   * @return change of error
   */
  FT merge(const std::size_t &px_enlarged, const std::size_t &px_merged) {
    // merge
    BOOST_FOREACH(face_descriptor f, px_facets[px_merged])
      seg_pmap[f] = px_enlarged;
    proxies[px_enlarged] = merged_px;
    // update facet proxy map
    proxies.erase(proxies.begin() + px_merged);
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      if (seg_pmap[f] > px_merged)
        --seg_pmap[f];
    }
  }

  /*!
   * @brief Merge adjacent regions, the re-fitting is performed.
   * @param range_of_proxies range of proxies, must be adjacent
   * @return change of error
   */
  FT merge(range_of_proxies) {}
  // Document what happens when the model has only one proxy,.

  /*!
   * @brief Split one proxy by default bisection, but N-section is also possible.
   * @param p proxy
   * @return change of error
   */
  FT split(Proxy &p, int n = 2) {}

  /*!
   * @brief Split range of proxies by default bisection, but N-section is also possible.
   * @param range_of_proxies range of proxies
   * @return change of error
   */
  FT split(range_of_proxies, int n = 2) {}

  /*!
   * @brief Meshing, choose the default area weighted or the PCA plane fitting.
   * @param[out] tm_out output triangle mesh
   * @param split_criterion boundary approximation recursively split criterion
   * @param pca_plane if use PCA plane fitting method
   * @return true if output triangle mesh is manifold,false otherwise.
   */
  bool meshing(TriangleMesh &tm_out, const FT split_criterion = 1, bool pca_plane = false) {

  }

  /*!
   * @brief Get the facet-proxy index map.
   * @tparam FacetProxyMap `WritablePropertyMap` with
   * `boost::graph_traits<TriangleMesh>::face_descriptor` as key and `std::size_t` as value type
   */
  template <typename FacetProxyMap>
  void get_proxy_map(FacetProxyMap &facet_proxy_map) {

  }

  /*!
   * @brief Get the proxies.
   * @return range of proxies.
   */
  std::pair<PxIterator, PxIterator> get_proxies();

// add function to return the proxy errors

  /*!
   * @brief Get the anchor points, which have the area-averaged position of the projected anchor vertex points on the incident proxies.
   * @return vector of anchor position points
   */
  const std::vector<Point_3> &get_anchor_points();

  /*!
   * @brief Get the anchor vertices.
   * @return vector of anchor vertices
   */
  const std::vector<vertex_descriptor> &get_anchor_vertices();

  /*!
   * @brief Get the indexed triangles, one triplet of integers per triangles, and that the integers refer to the anchor point indexes.
   * @return vector of indexed triangles.
   */
  const std::vector<std::size_t> &get_indexed_triangles();

private:
  /*!
   * @brief Random initialize proxies.
   * @param initial_px number of proxies
   * @return #proxies initialized
   */
  std::size_t seed_random(const std::size_t initial_px) {
    const std::size_t interval = num_faces(mesh) / initial_px;
    std::size_t index = 0;
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      if ((index++) % interval == 0) {
        proxies.push_back(fit_new_proxy(f));
      }
      if (proxies.size() >= initial_px)
        break;
    }
    return proxies.size();
    // std::cerr << initial_px << ' ' << proxies.size() << std::endl;
  }

  /*!
   * @brief Incremental initialize proxies.
   * @param initial_px number of proxies
   * @return #proxies initialized
   */
  std::size_t seed_incremental(const std::size_t initial_px) {
    // initialize a proxy and the proxy map to prepare for the insertion
    proxies.push_back(fit_new_proxy(*(faces(mesh).first)));
    BOOST_FOREACH(face_descriptor f, faces(mesh))
      seg_pmap[f] = 0;

    insert_proxy_furthest(initial_px - 1);
    return proxies.size();
  }

  /*!
   * @brief Hierarchical initialize proxies.
   * @param initial_px number of proxies
   * @return #proxies initialized
   */
  std::size_t seed_hierarchical(const std::size_t initial_px) {
    // initialize 2 proxy
    typename boost::graph_traits<TriangleMesh>::face_iterator
      fitr = faces(mesh).first;
    proxies.push_back(fit_new_proxy(*fitr));
    proxies.push_back(fit_new_proxy(*(++fitr)));

    const std::size_t num_steps = 5;
    while (proxies.size() < initial_px) {
      for (std::size_t i = 0; i < num_steps; ++i) {
        partition();
        fit();
      }

      // add proxies by error diffusion
      const std::size_t num_proxies = proxies.size();
      const std::size_t num_proxies_to_be_added =
        (num_proxies * 2 < initial_px) ? num_proxies : (initial_px - num_proxies);
      insert_proxy_hierarchical(num_proxies_to_be_added);
    }
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

    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      std::size_t px_idx = seg_pmap[f];
      FT err = fit_error(f, proxies[px_idx]);
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
   * @return #proxies inserted
   */
  std::size_t insert_proxy_furthest(const std::size_t num_proxies) {
    // when insert only one proxy, it has the same effect of insert_proxy_furthest()
    if (num_proxies == 0 || !insert_proxy_furthest())
      return 0;

    std::size_t num_inserted = 1;
    const std::size_t num_steps = 5;
    for (; num_inserted < num_proxies; ++num_inserted) {
      for (std::size_t i = 0; i < num_steps; ++i) {
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
    std::cout << "#px " << proxies.size() << std::endl;
    std::vector<FT> err(proxies.size(), FT(0));
    const FT sum_error = fitting_error(seg_pmap, err);
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
      FT to_add = (residual + pxe.fit_error) / avg_error;
      // floor_to_add maybe negative but no less than -1
      FT floor_to_add = FT(std::floor(CGAL::to_double(to_add)));
      const std::size_t q_to_add = static_cast<std::size_t>(CGAL::to_double(
        ((to_add - floor_to_add) > FT(0.5)) ? (floor_to_add + FT(1)) : floor_to_add));
      residual = (to_add - FT(static_cast<double>(q_to_add))) * avg_error;
      num_to_add[pxe.px_idx] = q_to_add;
    }
    for (std::size_t i = 0; i < px_error.size(); ++i)
      std::cout << "#px_id " << px_error[i].px_idx
        << ", #fit_error " << px_error[i].fit_error
        << ", #num_to_add " << num_to_add[px_error[i].px_idx] << std::endl;

    std::size_t num_inserted = 0;
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      const std::size_t px_id = seg_pmap[f];
      if (proxies[px_id].seed == f)
        continue;

      if (num_to_add[px_id] > 0) {
        proxies.push_back(fit_new_proxy(f));
        --num_to_add[px_id];
        ++num_inserted;
      }
    }
    std::cout << "#requested/inserted "
      << num_proxies_to_be_added << '/' << num_inserted << std::endl;

    return num_inserted;
  }

  /**
   * Find the best two regions to merge.
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
    BOOST_FOREACH(face_descriptor f, faces(mesh))
      px_facets[seg_pmap[f]].push_back(f);

    // find best merge
    MergedPair merged_set;
    // Proxy merged_px;
    FT min_merged_error = FT(0);
    bool first_merge = true;
    BOOST_FOREACH(edge_descriptor e, edges(mesh)) {
      if (CGAL::is_border(e, mesh))
        continue;
      std::size_t pxi = seg_pmap[face(halfedge(e, mesh), mesh)];
      std::size_t pxj = seg_pmap[face(opposite(halfedge(e, mesh), mesh), mesh)];
      if (pxi == pxj)
        continue;
      if (pxi > pxj)
        std::swap(pxi, pxj);
      if (merged_set.find(ProxyPair(pxi, pxj)) != merged_set.end())
        continue;

      std::list<face_descriptor> merged_patch(px_facets[pxi]);
      BOOST_FOREACH(face_descriptor f, px_facets[pxj])
        merged_patch.push_back(f);

      Proxy px = fit_new_proxy(merged_patch.begin(), merged_patch.end());
      FT sum_error(0);
      BOOST_FOREACH(face_descriptor f, merged_patch)
        sum_error += fit_error(f, px);
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
    fitting_error(seg_pmap, px_error);
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
   * @brief Fitting a new proxy.
   * 1. Fit proxy parameters from a list of facets.
   * 2. Set seed.
   * @tparam FacetIterator face_descriptor container iterator
   * @param beg container begin
   * @param end container end
   */
  template<typename FacetIterator>
  Proxy fit_new_proxy(const FacetIterator &beg, const FacetIterator &end) {
    CGAL_assertion(beg != end);

    // use proxy_fitting functor to fit proxy parameters
    Proxy px = proxy_fitting(beg, end);

    // find proxy seed
    px.seed = *beg;
    FT err_min = fit_error(*beg, px);
    std::pair<FacetIterator, FacetIterator> facets(beg, end);
    BOOST_FOREACH(face_descriptor f, facets) {
      FT err = fit_error(f, px);
      if (err < err_min) {
        err_min = err;
        px.seed = f;
      }
    }

    return px;
  }

  /*!
   * @brief Fitting a new proxy from a single facet.
   * 1. Fit proxy parameters from one facet.
   * 2. Set seed.
   * @param face_descriptor facet
   */
  Proxy fit_new_proxy(const face_descriptor &f) {
    std::vector<face_descriptor> fvec(1, f);
    // fit proxy parameters
    Proxy px = proxy_fitting(fvec.begin(), fvec.end());
    // find proxy seed
    px.seed = f;

    return px;
  }

  /*!
   * @brief Computes fitting error of a given partition @a seg_pmap.
   * @param seg_map facet partition index
   * @return total fitting error
   */
  FT fitting_error(const FacetSegmentMap &seg_pmap) {
    FT sum_error(0);
    BOOST_FOREACH(face_descriptor f, faces(mesh))
      sum_error += fit_error(f, proxies[seg_pmap[f]]);
    return sum_error;
  }

  /*!
   * @brief Computes fitting error of a given partition @a seg_pmap.
   * @param seg_map facet partition index
   * @param px_error vector of error of each proxy
   * @return total fitting error
   */
  FT fitting_error(const FacetSegmentMap &seg_pmap, std::vector<FT> &px_error) {
    FT sum_error(0);
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      const std::size_t px_idx = seg_pmap[f];
      FT err = fit_error(f, proxies[px_idx]);
      px_error[px_idx] += err;
      sum_error += err;
    }
    return sum_error;
  }
};

} // end namespace CGAL

#undef CGAL_NOT_TAGGED_ID

#endif // CGAL_VSA_APPROXIMATION
