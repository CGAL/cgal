#ifndef CGAL_VSA_H
#define CGAL_VSA_H

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
#include <set>
#include <queue>
#include <iostream>
#include <iterator>

#define CGAL_NOT_TAGGED_ID std::numeric_limits<std::size_t>::max()

namespace CGAL
{
/// @cond CGAL_DOCUMENT_INTERNAL
namespace internal
{
/**
 * @brief Main class for Variational Shape Approximation algorithm.
 *
 * @tparam TriangleMesh a CGAL TriangleMesh
 * @tparam GeomTraits a model of ApproximationGeomTraits
 */
template <typename TriangleMesh,
  typename ApproximationTrait,
  typename FacetAreaMap,
  typename VertexPointPmap>
  class VSA
{
  // type definitions
private:
  typedef typename ApproximationTrait::GeomTraits GeomTraits;
  typedef typename ApproximationTrait::Proxy Proxy;
  typedef typename ApproximationTrait::ErrorMetric ErrorMetric;
  typedef typename ApproximationTrait::ProxyFitting ProxyFitting;
  typedef typename ApproximationTrait::PlaneFitting PlaneFitting;

  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point;
  typedef typename GeomTraits::Vector_3 Vector;
  typedef typename GeomTraits::Plane_3 Plane;
  typedef typename GeomTraits::Construct_vector_3 Construct_vector_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename GeomTraits::Compute_scalar_product_3 Compute_scalar_product_3;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_iterator halfedge_iterator;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_iterator face_iterator;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_iterator vertex_iterator;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_iterator edge_iterator;

  typedef boost::associative_property_map<std::map<vertex_descriptor, int> > VertexStatusPMap;
  typedef boost::associative_property_map<std::map<halfedge_descriptor, int> > HalfedgeStatusPMap;

  typedef std::vector<halfedge_descriptor> ChordVector;
  typedef typename ChordVector::iterator ChordVectorIterator;

public:
  enum Initialization {
    RandomInit,
    IncrementalInit,
    HierarchicalInit
  };
  
  // The average positioned anchor attached to a vertex.
  struct Anchor {
    vertex_descriptor vtx; // The associated vertex.
    Point pos; // The position of the anchor.
  };

  // The border cycle of a region.
  // One region may have multiple border cycles.
  struct Border {
    Border(const halfedge_descriptor &h)
      : he_head(h), num_anchors(0) {}

    halfedge_descriptor he_head; // The heading halfedge of the border cylce.
    std::size_t num_anchors; // The number of anchors on the border.
  };

  // member variables
private:
  const TriangleMesh &mesh;
  const VertexPointPmap &vertex_point_pmap;
  Construct_vector_3 vector_functor;
  Construct_scaled_vector_3 scale_functor;
  Construct_sum_of_vectors_3 sum_functor;
  Compute_scalar_product_3 scalar_product_functor;

  // Proxy.
  std::vector<Proxy> proxies;
  
  // The attached anchor index of a vertex.
  std::map<vertex_descriptor, int> vertex_status_map;
  VertexStatusPMap vertex_status_pmap;

  // Mesh facet area map, for anchor position average.
  const FacetAreaMap area_pmap;

  // All anchors.
  std::size_t anchor_index;
  std::vector<Anchor> anchors;

  // All borders cycles.
  std::vector<Border> borders;

  // The error metric.
  ErrorMetric fit_error;

  // The proxy fitting functor.
  ProxyFitting proxy_fitting;

  // The proxy plane fitting functor.
  PlaneFitting plane_fitting;

  //member functions
public:
  /**
   * Initialize and prepare for the approximation.
   * @pre @a TriangleMesh.is_pure_triangle()
   * @param _mesh `CGAL TriangleMesh` on which approximation operate.
   * @param _vertex_point_map vertex point map of the input mesh.
   * @param _traits geometric trait object.
   */
  VSA(const TriangleMesh &_mesh,
    const ApproximationTrait &_appx_trait,
    const VertexPointPmap &_vertex_point_map,
    const FacetAreaMap &_facet_area_map)
    : mesh(_mesh),
    vertex_point_pmap(_vertex_point_map),
    area_pmap(_facet_area_map),
    vertex_status_pmap(vertex_status_map),
    fit_error(_appx_trait.construct_fit_error_functor()),
    proxy_fitting(_appx_trait.construct_proxy_fitting_functor()),
    plane_fitting(_appx_trait.construct_plane_fitting_functor()) {

    GeomTraits traits;
    vector_functor = traits.construct_vector_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scalar_product_functor = traits.compute_scalar_product_3_object();

    // initialize all vertex anchor status
    enum Vertex_status { NO_ANCHOR = -1 };
    BOOST_FOREACH(vertex_descriptor v, vertices(mesh))
      vertex_status_map.insert(std::pair<vertex_descriptor, int>(v, static_cast<int>(NO_ANCHOR)));
  }

  /**
   * Partitions the mesh into the designated number of regions, and stores them in @a seg_map.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param number_of_segments number of designated proxies
   * @param number_of_iterations number of iterations when fitting the proxy
   * @param[out] seg_map facet partition index
   */
  template<typename FacetSegmentMap>
  void partition(const std::size_t number_of_segments, const std::size_t number_of_iterations, FacetSegmentMap &seg_pmap) {
    random_seed(number_of_segments);
    //random_seed(number_of_segments / 2);
    for (std::size_t i = 0; i < number_of_iterations; ++i) {
      flooding(seg_pmap);
      fitting(seg_pmap);
    }
  }

  /**
   * Partitions the mesh incrementally into the designated number of regions, and stores them in @a seg_map.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param number_of_segments number of designated proxies
   * @param number_of_iterations number of iterations when fitting the proxy
   * @param[out] seg_map facet partition index
   */
  template<typename FacetSegmentMap>
  void partition_incre(const std::size_t number_of_segments, const std::size_t number_of_iterations, FacetSegmentMap &seg_pmap) {
    // random_seed(number_of_segments);
    random_seed(number_of_segments / 2);

    for (std::size_t i = 0; i < number_of_iterations; ++i) {
      flooding(seg_pmap);
      fitting(seg_pmap);
    }

    while (proxies.size() < number_of_segments) {
      insert_proxy(seg_pmap);
      for (std::size_t i = 0; i < number_of_iterations; ++i) {
        flooding(seg_pmap);
        fitting(seg_pmap);
      }
    }
  }

  /**
   * Partitions the mesh into the designated number of regions, and stores them in @a seg_map.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param number_of_segments number of designated proxies
   * @param number_of_iterations number of iterations when fitting the proxy
   * @param[out] seg_map facet partition index
   */
  template<typename FacetSegmentMap>
  void partition_hierarchical(const std::size_t number_of_segments, const std::size_t number_of_iterations, FacetSegmentMap &seg_pmap) {
    hierarchical_seed(number_of_segments, seg_pmap);
    for (std::size_t i = 0; i < number_of_iterations; ++i) {
      flooding(seg_pmap);
      fitting(seg_pmap);
    }
  }

  /**
   * Extracts the approximated triangle mesh from a partition @a seg_pmap, and stores the triangles in @a tris.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param seg_map facet partition index
   * @param[out] tris indexed triangles
   */
  template<typename FacetSegmentMap>
  void extract_mesh(const FacetSegmentMap &seg_pmap, std::vector<int> &tris) {
    anchor_index = 0;
    find_anchors(seg_pmap);
    find_edges(seg_pmap);
    add_anchors(seg_pmap);

    pseudo_CDT(seg_pmap, tris);

    std::vector<Plane> px_planes(proxies.size());
    std::vector<std::list<face_descriptor> > px_facets(proxies.size());
    BOOST_FOREACH(face_descriptor f, faces(mesh))
      px_facets[seg_pmap[f]].push_back(f);
    for (std::size_t i = 0; i < proxies.size(); ++i)
      px_planes[i] = plane_fitting(px_facets[i].begin(), px_facets[i].end());
    compute_anchor_position(seg_pmap, px_planes);
    std::vector<Point> vtx;
    BOOST_FOREACH(const Anchor &a, anchors)
      vtx.push_back(a.pos);
    if (is_manifold_surface(tris, vtx))
      std::cout << "Manifold surface." << std::endl;
    else
      std::cout << "Non-manifold surface." << std::endl;
  }

  /**
   * Use a incremental builder to test if the indexed triangle surface is manifold
   * @param tris indexed triangles
   * @param vtx vertex positions
   * @return true if build successfully
   */
  bool is_manifold_surface(const std::vector<int> &tris, const std::vector<Point> &vtx) {
    typedef CGAL::Polyhedron_3<GeomTraits> PolyhedronSurface;
    typedef typename PolyhedronSurface::HalfedgeDS HDS;
    
    HDS hds;
    CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);
    builder.begin_surface(vtx.size(), tris.size() / 3);
    BOOST_FOREACH(const Point &v, vtx)
      builder.add_vertex(v);
    for (std::vector<int>::const_iterator itr = tris.begin(); itr != tris.end(); itr += 3) {
      if (builder.test_facet(itr, itr + 3)) {
        builder.begin_facet();
        builder.add_vertex_to_facet(*itr);
        builder.add_vertex_to_facet(*(itr + 1));
        builder.add_vertex_to_facet(*(itr + 2));
        builder.end_facet();
      }
      else {
        // std::cerr << "test_facet failed" << std::endl;
        builder.end_surface();
        return false;
      }
    }
    builder.end_surface();

    return true;
  }

  /**
   * Collect the anchors.
   * @return vector of anchors
   */
  std::vector<Anchor> collect_anchors() {
    return anchors;
  }

  /**
   * Collect the partition @a seg_pmap approximated borders.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param seg_map facet partition index
   * @return anchor indexes of each border
   */
  template<typename FacetSegmentMap>
  std::vector<std::vector<std::size_t> >
  collect_borders(const FacetSegmentMap &seg_pmap) {
    std::vector<std::vector<std::size_t> > bdrs;
    for (typename std::vector<Border>::iterator bitr = borders.begin();
      bitr != borders.end(); ++bitr) {
      std::vector<std::size_t> bdr;
      const halfedge_descriptor he_mark = bitr->he_head;
      halfedge_descriptor he = he_mark;
      do {
        ChordVector chord;
        walk_to_next_anchor(he, chord, seg_pmap);
        bdr.push_back(vertex_status_pmap[target(he, mesh)]);
      } while(he != he_mark);
      bdrs.push_back(bdr);
    }
    return bdrs;
  }

private:
  /**
   * Random initialize proxies.
   * @param initial_px number of proxies
   */
  void random_seed(const std::size_t initial_px) {
    proxies.clear();
    const std::size_t interval = num_faces(mesh) / initial_px;
    std::size_t index = 0;
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      if ((index++) % interval == 0) {
        // Use proxy_fitting functor to create a proxy
        std::vector<face_descriptor> fvec(1, f);
        proxies.push_back(proxy_fitting(fvec.begin(), fvec.end()));
      }
      if (proxies.size() >= initial_px)
        break;
    }
    std::cerr << initial_px << ' ' << proxies.size() << std::endl;
  }

  /**
   * Hierarchical initialize proxies.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param initial_px number of proxies
   * @param[out] seg_map facet partition index
   */
  template<typename FacetSegmentMap>
  void hierarchical_seed(const std::size_t initial_px, FacetSegmentMap &seg_pmap) {
    if (initial_px == 0)
      return;
    if (num_faces(mesh) < 2)
      return;

    proxies.clear();
    // generate 2 seeds
    face_iterator f0 = faces(mesh).first, f1 = ++f0;
    std::vector<face_descriptor> fvec(1, *f0);
    proxies.push_back(proxy_fitting(fvec.begin(), fvec.end()));
    std::vector<face_descriptor> fvec2(1, *f1);
    proxies.push_back(proxy_fitting(fvec2.begin(), fvec2.end()));

    const std::size_t num_steps = 5;
    while (proxies.size() < initial_px) {
      for (std::size_t i = 0; i < num_steps; ++i) {
        flooding(seg_pmap);
        fitting(seg_pmap);
      }

      // add proxies by error diffusion
      const std::size_t num_proxies = proxies.size();
      const std::size_t num_proxies_to_be_added =
        (num_proxies * 2 < initial_px) ? num_proxies : (initial_px - num_proxies);
      insert_proxy_error_diffusion(num_proxies_to_be_added, seg_pmap);
    }
  }

  /**
   * Propagates the proxy seed facets and floods the whole mesh to minimize the fitting error.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param[out] seg_map facet partition index
   */
  template<typename FacetSegmentMap>
  void flooding(FacetSegmentMap &seg_pmap) {
    // The facet candidate to be queued.
    struct FacetToIntegrate {
      face_descriptor f;
      std::size_t i;
      FT fit_error;
      bool operator<(const FacetToIntegrate &rhs) const {
        return fit_error > rhs.fit_error;
      }
    };

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

  /**
   * Calculates and updates the best fitting proxies of a given partition @a seg_pmap.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param seg_map facet partition index
   */
  template <typename FacetSegmentMap>
  void fitting(const FacetSegmentMap &seg_pmap) {
    std::vector<std::list<face_descriptor> > px_facets(proxies.size());
    BOOST_FOREACH(face_descriptor f, faces(mesh))
      px_facets[seg_pmap[f]].push_back(f);

    for (std::size_t i = 0; i < proxies.size(); ++i)
      proxies[i] = proxy_fitting(px_facets[i].begin(), px_facets[i].end());
  }

  /**
   * Inserts a proxy of a given partition @a seg_pmap at the furthest facet of the region with the maximum fitting error.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param seg_map facet partition index
   */
  template <typename FacetSegmentMap>
  void insert_proxy(const FacetSegmentMap &seg_pmap) {
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

    // create new proxy
    std::vector<face_descriptor> fvec(1, max_facet[max_px_idx]);
    proxies.push_back(proxy_fitting(fvec.begin(), fvec.end()));
  }

  /**
   * Add proxies by diffusing fitting error into current partitions.
   * Each partition is added with the number of proxies in proportional to its fitting error.
   * Note that the number of inserted proxies doesn't necessarily equal the requested number.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param num_proxies_to_be_added added number of proxies
   * @param seg_map facet partition index
   * @return inserted number of proxies
   */
  template<typename FacetSegmentMap>
  std::size_t insert_proxy_error_diffusion(const std::size_t num_proxies_to_be_added, const FacetSegmentMap &seg_pmap) {
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
        std::vector<face_descriptor> fvec(1, f);
        proxies.push_back(proxy_fitting(fvec.begin(), fvec.end()));
        --num_to_add[px_id];
        ++num_inserted;
      }
    }
    std::cout << "#requested/inserted "
      << num_proxies_to_be_added << '/' << num_inserted << std::endl;

    return num_inserted;
  }

  /**
   * Finds the anchors of a given partition @a seg_pmap.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param seg_map facet partition index
   */
  template<typename FacetSegmentMap>
  void find_anchors(const FacetSegmentMap &seg_pmap) {
    anchors.clear();

    BOOST_FOREACH(vertex_descriptor vtx, vertices(mesh)) {
      std::size_t border_count = 0;

      BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(vtx, mesh)) {
        if (CGAL::is_border_edge(h, mesh))
          ++border_count;
        else if (seg_pmap[face(h, mesh)] != seg_pmap[face(opposite(h, mesh), mesh)])
          ++border_count;
      }
      if (border_count >= 3)
        attach_anchor(vtx);
    }
  }

  /**
   * Finds and approximates the edges connecting the anchors of a given partition @a seg_pmap.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param seg_map facet partition index
   */
  template<typename FacetSegmentMap>
  void find_edges(const FacetSegmentMap &seg_pmap) {
    // collect candidate halfedges in a set
    std::set<halfedge_descriptor> he_candidates;
    BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh)) {
      if (!CGAL::is_border(h, mesh)
        && (CGAL::is_border(opposite(h, mesh), mesh)
          || seg_pmap[face(h, mesh)] != seg_pmap[face(opposite(h, mesh), mesh)]))
        he_candidates.insert(h);
    }

    // pick up one candidate halfedge each time and traverse the connected border
    borders.clear();
    while (!he_candidates.empty()) {
      halfedge_descriptor he_start = *he_candidates.begin();
      walk_to_first_anchor(he_start, seg_pmap);
      // no anchor in this connected border, make a new anchor
      if (!is_anchor_attached(he_start))
        attach_anchor(he_start);

      // a new connected border
      borders.push_back(Border(he_start));
      std::cerr << "#border " << borders.size() << std::endl;
      const halfedge_descriptor he_mark = he_start;
      do {
        ChordVector chord;
        walk_to_next_anchor(he_start, chord, seg_pmap);
        borders.back().num_anchors += subdivide_chord(chord.begin(), chord.end(), seg_pmap);
        std::cerr << "#chord_anchor " << borders.back().num_anchors << std::endl;

        for (ChordVectorIterator citr = chord.begin(); citr != chord.end(); ++citr)
          he_candidates.erase(*citr);
      } while (he_start != he_mark);
    }
  }

  /**
   * Adds anchors to the border cycles with only 2 anchors of a given partition @a seg_pmap.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param seg_map facet partition index
   */
  template<typename FacetSegmentMap>
  void add_anchors(const FacetSegmentMap &seg_pmap) {
    typedef typename std::vector<Border>::iterator BorderIterator;
    for (BorderIterator bitr = borders.begin(); bitr != borders.end(); ++bitr) {
      if (bitr->num_anchors > 2)
        continue;

      // 2 initial anchors at least
      CGAL_assertion(bitr->num_anchors == 2);
      // borders with only 2 initial anchors
      const halfedge_descriptor he_mark = bitr->he_head;
      Point pt_begin = vertex_point_pmap[target(he_mark, mesh)];
      Point pt_end = pt_begin;

      halfedge_descriptor he = he_mark;
      ChordVector chord;
      std::size_t count = 0;
      do {
        walk_to_next_border_halfedge(he, seg_pmap);
        if (!is_anchor_attached(he))
          chord.push_back(he);
        else {
          if (count == 0)
            pt_end = vertex_point_pmap[target(he, mesh)];
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
      Vector chord_vec = vector_functor(pt_begin, pt_end);
      chord_vec = scale_functor(chord_vec,
        FT(1.0 / std::sqrt(CGAL::to_double(chord_vec.squared_length()))));
      for (ChordVectorIterator citr = chord.begin(); citr != chord.end(); ++citr) {
        Vector vec = vector_functor(pt_begin, vertex_point_pmap[target(*citr, mesh)]);
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

  /**
   * Runs the pseudo Constrained Delaunay Triangulation at each region of a given partition @a seg_pmap, and stores the extracted indexed triangles in @a tris.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param seg_map facet partition index
   * @param tris extracted tirangles, index of anchors
   */
  template<typename FacetSegmentMap>
  void pseudo_CDT(const FacetSegmentMap &seg_pmap,
    std::vector<int> &tris) {
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
    typedef typename SubGraph::edge_descriptor sg_edge_descriptor;
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
    BOOST_FOREACH(vertex_descriptor v, vertices(mesh)) {
      sg_vertex_descriptor sgv = add_vertex(gmain);
      global_vanchor_map[sgv] = vertex_status_pmap[v];
      global_vtag_map[sgv] = vertex_status_pmap[v];
      vmap.insert(std::pair<vertex_descriptor, sg_vertex_descriptor>(v, sgv));
    }
    BOOST_FOREACH(edge_descriptor e, edges(mesh)) {
      vertex_descriptor vs = source(e, mesh);
      vertex_descriptor vt = target(e, mesh);
      FT len(std::sqrt(CGAL::to_double(
        CGAL::squared_distance(vertex_point_pmap[vs], vertex_point_pmap[vt]))));
      add_edge(to_sgv_map[vs], to_sgv_map[vt], len, gmain);
    }

    std::vector<VertexVector> vertex_patches(proxies.size());
    BOOST_FOREACH(vertex_descriptor v, vertices(mesh)) {
      std::set<std::size_t> px_set;
      BOOST_FOREACH(face_descriptor f, faces_around_target(halfedge(v, mesh), mesh)) {
        if (f != boost::graph_traits<TriangleMesh>::null_face())
          px_set.insert(seg_pmap[f]);
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
        walk_to_next_anchor(he, chord, seg_pmap);

        std::vector<FT> vdist;
        vdist.push_back(FT(0));
        BOOST_FOREACH(halfedge_descriptor h, chord) {
          FT elen = global_eweight_map[edge(
            to_sgv_map[source(h, mesh)],
            to_sgv_map[target(h, mesh)],
            gmain).first];
          vdist.push_back(vdist.back() + elen);
        }

        FT half_chord_len = vdist.back() / FT(2);
        const int anchorleft = vertex_status_pmap[source(chord.front(), mesh)];
        const int anchorright = vertex_status_pmap[target(chord.back(), mesh)];
        typename std::vector<FT>::iterator ditr = vdist.begin() + 1;
        for (typename ChordVector::iterator hitr = chord.begin();
          hitr != chord.end() - 1; ++hitr, ++ditr) {
          if (*ditr < half_chord_len)
            global_vtag_map[to_sgv_map[target(*hitr, mesh)]] = anchorleft;
          else
            global_vtag_map[to_sgv_map[target(*hitr, mesh)]] = anchorright;
        }
      } while(he != he_mark);
    }

    // collect triangles
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      halfedge_descriptor he = halfedge(f, mesh);
      int i = global_vtag_map[to_sgv_map[source(he, mesh)]];
      int j = global_vtag_map[to_sgv_map[target(he, mesh)]];
      int k = global_vtag_map[to_sgv_map[target(next(he, mesh), mesh)]];
      if (i != j && i != k && j != k) {
        tris.push_back(i);
        tris.push_back(j);
        tris.push_back(k);
      }
    }
  }

  /**
   * Computes fitting error of a given partition @a seg_pmap.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param seg_map facet partition index
   * @return total fitting error
   */
  template<typename FacetSegmentMap>
  FT fitting_error(const FacetSegmentMap &seg_pmap) {
    FT sum_error(0);
    BOOST_FOREACH(face_descriptor f, faces(mesh))
      sum_error += fit_error(f, proxies[seg_pmap[f]]);
    return sum_error;
  }

  /**
   * Computes fitting error of a given partition @a seg_pmap.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param seg_map facet partition index
   * @param px_error vector of error of each proxy
   * @return total fitting error
   */
  template<typename FacetSegmentMap>
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

  /**
   * Computes and updates the proxy centers of a given partition @a seg_pmap.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param seg_map facet partition index
   * @param proxies_center proxy centers
   */
  template<typename FacetSegmentMap>
  void compute_proxy_center(const FacetSegmentMap &seg_pmap, std::vector<Point> &proxies_center) {
    proxies_center = std::vector<Point>(proxies.size());

    std::vector<Vector> centers(proxies.size(), CGAL::NULL_VECTOR);
    std::vector<FT> areas(proxies.size(), FT(0.0));
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      const std::size_t px_idx = seg_pmap[f];
      const halfedge_descriptor he = halfedge(f, mesh);
      Point pt = CGAL::centroid(
        vertex_point_pmap[source(he, mesh)],
        vertex_point_pmap[target(he, mesh)],
        vertex_point_pmap[target(next(he, mesh), mesh)]);
      Vector vec = vector_functor(CGAL::ORIGIN, pt);
      FT area = area_pmap[f];
      areas[px_idx] += area;
      centers[px_idx] = sum_functor(centers[px_idx], scale_functor(vec, area));
    }

    for (std::size_t i = 0; i < proxies.size(); ++i) {
      Vector vec = scale_functor(centers[i], FT(1.0) / areas[i]);
      proxies_center[i] = Point(vec.x(), vec.y(), vec.z());
    }
  }

  /**
   * Computes and updates the proxy areas of a given partition @a seg_pmap.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param seg_map facet partition index
   * @param proxies_area proxy areas
   */
  template<typename FacetSegmentMap>
  void compute_proxy_area(const FacetSegmentMap &seg_pmap, std::vector<FT> &proxies_area) {
    std::vector<FT> areas(proxies.size(), FT(0));
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      areas[seg_pmap[f]] += area_pmap[f];
    }
    proxies_area.swap(areas);
  }

  /**
   * Walks along the region border to the first halfedge pointing to a vertex associated with an anchor of a given partition @a seg_pmap.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param[in/out] he_start region border halfedge
   * @param seg_map facet partition index
   */
  template<typename FacetSegmentMap>
  void walk_to_first_anchor(halfedge_descriptor &he_start, const FacetSegmentMap &seg_pmap) {
    const halfedge_descriptor start_mark = he_start;
    while (!is_anchor_attached(he_start)) {
      // no anchor attached to the halfedge target
      walk_to_next_border_halfedge(he_start, seg_pmap);
      if (he_start == start_mark) // back to where started, a circular border
        return;
    }
  }

  /**
   * Walks along the region border to the next anchor and records the path as @a chord of a given partition @a seg_pmap.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param[in/out] he_start starting region border halfedge pointing to a vertex associated with an anchor
   * @param[out] chord recorded path chord
   * @param seg_map facet partition index
   */
  template<typename FacetSegmentMap>
  void walk_to_next_anchor(
    halfedge_descriptor &he_start,
    ChordVector &chord,
    const FacetSegmentMap &seg_pmap) {
    do {
      walk_to_next_border_halfedge(he_start, seg_pmap);
      chord.push_back(he_start);
    } while (!is_anchor_attached(he_start));
  }

  /**
   * Walks to next border halfedge of a given partition @a seg_pmap.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param[in/out] he_start region border halfedge
   * @param seg_map facet partition index
   */
  template<typename FacetSegmentMap>
  void walk_to_next_border_halfedge(halfedge_descriptor &he_start, const FacetSegmentMap &seg_pmap) {
    const std::size_t px_idx = seg_pmap[face(he_start, mesh)];
    BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(he_start, mesh)) {
      if (CGAL::is_border(h, mesh) || seg_pmap[face(h, mesh)] != px_idx) {
        he_start = opposite(h, mesh);
        return;
      }
    }
  }

  /**
   * Subdivides a chord recursively in range [@a chord_begin, @a chord_end) of a given partition @a seg_pmap.
   * @tparam FacetSegmentMap `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::face_handle` as key and `std::size_t` as value type
   * @param chord_begin begin iterator of the chord
   * @param chord_end end iterator of the chord
   * @param seg_map facet partition index
   * @return the number of anchors of the chord apart from the first one
   */
  template<typename FacetSegmentMap>
  std::size_t subdivide_chord(
    const ChordVectorIterator &chord_begin,
    const ChordVectorIterator &chord_end,
    const FacetSegmentMap &seg_pmap,
    const FT thre = FT(0.2)) {
    const std::size_t chord_size = std::distance(chord_begin, chord_end);
    // do not subdivide trivial chord
    if (chord_size < 4)
      return 1;

    halfedge_descriptor he_start = *chord_begin;
    std::size_t px_left = seg_pmap[face(he_start, mesh)];
    std::size_t px_right = px_left;
    if (!CGAL::is_border(opposite(he_start, mesh), mesh))
      px_right = seg_pmap[face(opposite(he_start, mesh), mesh)];

    // suppose the proxy normal angle is acute
    FT norm_sin(1.0);
    if (!CGAL::is_border(opposite(he_start, mesh), mesh)) {
      Vector vec = CGAL::cross_product(proxies[px_left].normal, proxies[px_right].normal);
      norm_sin = FT(std::sqrt(CGAL::to_double(scalar_product_functor(vec, vec))));
    }
    FT criterion = thre + FT(1.0);

    ChordVectorIterator he_max;
    const ChordVectorIterator chord_last = chord_end - 1;
    std::size_t anchor_begin = vertex_status_pmap[source(he_start, mesh)];
    std::size_t anchor_end = vertex_status_pmap[target(*chord_last, mesh)];
    const Point &pt_begin = vertex_point_pmap[source(he_start, mesh)];
    const Point &pt_end = vertex_point_pmap[target(*chord_last, mesh)];
    if (anchor_begin == anchor_end) {
      // circular chord
      CGAL_assertion(chord_size > 2);
      // if (chord_size < 3)
      //   return;

      FT dist_max(0.0);
      for (ChordVectorIterator citr = chord_begin; citr != chord_last; ++citr) {
        FT dist = CGAL::squared_distance(pt_begin, vertex_point_pmap[target(*citr, mesh)]);
        dist = FT(std::sqrt(CGAL::to_double(dist)));
        if (dist > dist_max) {
          he_max = citr;
          dist_max = dist;
        }
      }
    }
    else {
      FT dist_max(0.0);
      Vector chord_vec = vector_functor(pt_begin, pt_end);
      FT chord_len(std::sqrt(CGAL::to_double(chord_vec.squared_length())));
      chord_vec = scale_functor(chord_vec, FT(1.0) / chord_len);

      for (ChordVectorIterator citr = chord_begin; citr != chord_last; ++citr) {
        Vector vec = vector_functor(pt_begin, vertex_point_pmap[target(*citr, mesh)]);
        vec = CGAL::cross_product(chord_vec, vec);
        FT dist(std::sqrt(CGAL::to_double(vec.squared_length())));
        if (dist > dist_max) {
          he_max = citr;
          dist_max = dist;
        }
      }

      criterion = dist_max * norm_sin / chord_len;
    }

    if (criterion > thre) {
      // subdivide at the most remote vertex
      attach_anchor(*he_max);

      std::size_t num0 = subdivide_chord(chord_begin, he_max + 1, seg_pmap);
      std::size_t num1 = subdivide_chord(he_max + 1, chord_end, seg_pmap);

      return num0 + num1;
    }

    return 1;
  }

  /**
   * Check if the target vertex of a halfedge is attached with an anchor.
   * @param he halfedge
   */
  bool is_anchor_attached(const halfedge_descriptor &he) {
    return is_anchor_attached(target(he, mesh), vertex_status_pmap);
  }

  /**
   * Check if a vertex is attached with an anchor.
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

  /**
   * Attachs an anchor to the vertex.
   * @param vtx vertex
   */
  void attach_anchor(const vertex_descriptor &vtx) {
    vertex_status_pmap[vtx] = static_cast<int>(anchor_index++);
  }

  /**
   * Attachs an anchor to the target vertex of the halfedge.
   * @param he halfedge
   */
  void attach_anchor(const halfedge_descriptor &he) {
    vertex_descriptor vtx = target(he, mesh);
    attach_anchor(vtx);
  }

  /**
   * Calculate the anchor positions.
   * @param 
   */
  template<typename FacetSegmentMap>
  void compute_anchor_position(
    const FacetSegmentMap &seg_pmap,
    const std::vector<Plane> &px_planes) {
    std::vector<FT> proxies_area; // The proxy area.
    compute_proxy_area(seg_pmap, proxies_area);

    anchors = std::vector<Anchor>(anchor_index);
    BOOST_FOREACH(vertex_descriptor v, vertices(mesh)) {
      if (is_anchor_attached(v, vertex_status_pmap)) {
        // construct an anchor from vertex and the incident proxies
        std::set<std::size_t> px_set;
        BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(v, mesh)) {
          if (!CGAL::is_border(h, mesh))
            px_set.insert(seg_pmap[face(h, mesh)]);
        }

        // construct an anchor from vertex and the incident proxies
        FT avgx(0), avgy(0), avgz(0), sum_area(0);
        const Point vtx_pt = vertex_point_pmap[v];
        for (std::set<std::size_t>::iterator pxitr = px_set.begin();
          pxitr != px_set.end(); ++pxitr) {
          std::size_t px_idx = *pxitr;
          Point proj = px_planes[px_idx].projection(vtx_pt);
          FT area = proxies_area[px_idx];
          avgx += proj.x() * area;
          avgy += proj.y() * area;
          avgz += proj.z() * area;
          sum_area += area;
        }
        Point pos = Point(avgx / sum_area, avgy / sum_area, avgz / sum_area);
        std::size_t aidx = vertex_status_pmap[v];
        anchors[aidx].vtx = v;
        anchors[aidx].pos = pos;

        // FT sum_area(0);
        // Vector pos = CGAL::NULL_VECTOR;
        // const Point vtx_pt = vertex_point_pmap[v];
        // BOOST_FOREACH(const std::size_t &px_idx, px_set) {
        //   Vector proj = vector_functor(CGAL::ORIGIN, proxies[px_idx].fit_plane.projection(vtx_pt));
        //   FT area = proxies_area[px_idx];
        //   pos = sum_functor(pos, scale_functor(proj, area));
        //   sum_area += area;
        // }
        // pos = scale_functor(pos, FT(1) / sum_area);
        // std::size_t aidx = vertex_status_pmap[v];
        // anchors[aidx].vtx = v;
        // anchors[aidx].pos = Point(CGAL::ORIGIN) + pos;
      }
    }
  }
}; // end class VSA
} // end namespace internal
} // end namespace CGAL

#undef CGAL_NOT_TAGGED_ID

#endif // CGAL_VSA_H
