#ifndef CGAL_VSA_SEGMENTATION_H
#define CGAL_VSA_SEGMENTATION_H

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/foreach.hpp>

#include <vector>
#include <cmath>
#include <map>
#include <set>
#include <iostream>
#include <iterator>

#define CGAL_NOT_TAGGED_ID std::numeric_limits<std::size_t>::max()

namespace CGAL
{
/// @cond CGAL_DOCUMENT_INTERNAL
namespace internal
{
/**
 * @brief Main entry point for VSA mesh segmentation algorithm.
 *
 * blah blah...
 *
 * @tparam Polyhedron a CGAL polyhedron
 * @tparam GeomTraits a model of SegmentationGeomTraits
 */
template <typename Polyhedron,
  typename GeomTraits,
  //typename FacetSegmentMap,
  typename VertexPointPmap>
  class VSA_segmentation
{
private:
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point;
  typedef typename GeomTraits::Vector_3 Vector;
  typedef typename GeomTraits::Plane_3 Plane;
  typedef typename GeomTraits::Construct_vector_3 Construct_vector_3;
  typedef typename GeomTraits::Construct_normal_3 Construct_normal_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename GeomTraits::Compute_squared_area_3 Compute_squared_area_3;
  typedef typename GeomTraits::Compute_scalar_product_3 Compute_scalar_product_3;

  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::halfedge_iterator halfedge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::face_iterator face_iterator;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_iterator vertex_iterator;
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::edge_iterator edge_iterator;

  //typedef boost::associative_property_map<std::map<face_descriptor, std::size_t> > FacetSegmentMap;
  typedef boost::associative_property_map<std::map<face_descriptor, Vector> > FacetNormalMap;
  typedef boost::associative_property_map<std::map<face_descriptor, FT> > FacetAreaMap;
  typedef boost::associative_property_map<std::map<vertex_descriptor, int> > VertexStatusPMap;
  typedef boost::associative_property_map<std::map<halfedge_descriptor, int> > HalfedgeStatusPMap;

  typedef std::vector<halfedge_descriptor> ChordVector;
  typedef typename ChordVector::iterator ChordVectorIterator;

  enum Vertex_status {
    NO_ANCHOR = -1 // vertex v has no anchor attached
  };

  enum Halfedge_status {
    OFF_BORDER, // halfedge h is off proxy border
    CANDIDATE, // halfedge h is on proxy border, waiting to be visited
    ON_BORDER // halfedge h is on proxy border
  };

public:
  struct PlaneProxy {
    Vector normal;
    face_descriptor seed;
  };

  struct FacetToIntegrate {
    face_descriptor f;
    std::size_t i;
    FT fit_error;

    bool operator<(const FacetToIntegrate &rhs) const {
      return this->fit_error < rhs.fit_error;
    }
  };

  /*template <typename ShapeProxy,
    typename FacetNormalMap = boost::associative_property_map<std::map<face_descriptor, Vector> >,
    typename FacetAreaMap = boost::associative_property_map<std::map<face_descriptor, FT> >>*/
  /*template <typename FacetNormalMap,
    typename FacetAreaMap>*/
  struct L21Metric {
    //L21Metric() {}

    L21Metric(GeomTraits traits,
      const FacetNormalMap &normal_pmap,
      const FacetAreaMap &area_pmap)
      : scalar_product_functor(traits.compute_scalar_product_3_object()),
      sum_functor(traits.construct_sum_of_vectors_3_object()),
      scale_functor(traits.construct_scaled_vector_3_object()),
      normal_pmap(normal_pmap),
      area_pmap(area_pmap)
      {}

    FT operator()(const face_descriptor &f, const PlaneProxy &px) {
      Vector v = sum_functor(normal_pmap[f], scale_functor(px.normal, FT(-1)));
      return area_pmap[f] * scalar_product_functor(v, v);
    }

    const FacetNormalMap &normal_pmap;
    const FacetAreaMap &area_pmap;
    Construct_scaled_vector_3 scale_functor;
    Compute_scalar_product_3 scalar_product_functor;
    Construct_sum_of_vectors_3 sum_functor;
  };

  // Anchor
  struct Anchor {
    Anchor(const vertex_descriptor &vtx_, const Point &pos_)
    : vtx(vtx_), pos(pos_) {}

    vertex_descriptor vtx; // The associated vertex
    Point pos; // The position of the anchor
  };

  // Border
  struct Border {
    Border(const halfedge_descriptor &h)
      : he_head(h), num_anchors(0) {}

    halfedge_descriptor he_head; // The heading halfedge of the border
    std::size_t num_anchors; // The number of anchors of the border
  };

// member variables
private:
  const Polyhedron &mesh;
  const VertexPointPmap &vertex_point_pmap;
  GeomTraits traits;
  Construct_vector_3 vector_functor;
  Construct_normal_3 normal_functor;
  Construct_scaled_vector_3 scale_functor;
  Construct_sum_of_vectors_3 sum_functor;
  Compute_scalar_product_3 scalar_product_functor;
  Compute_squared_area_3 area_functor;

  std::vector<PlaneProxy> proxies;
  // proxy auxiliary information
  std::vector<Point> proxies_center; // proxy center
  std::vector<FT> proxies_area; // proxy area
  
  // the lifetime of the container must encompass the use of the adaptor
  std::map<face_descriptor, Vector> facet_normals;
  FacetNormalMap normal_pmap;
  std::map<face_descriptor, FT> facet_areas;
  FacetAreaMap area_pmap;

  std::map<vertex_descriptor, int> vertex_status_map;
  VertexStatusPMap vertex_status_pmap;
  std::map<halfedge_descriptor, int> halfedge_status_map;
  HalfedgeStatusPMap halfedge_status_pmap;

  // All anchors
  std::vector<Anchor> anchors;

  // All borders, some regions may have multiple borders
  std::vector<Border> borders;

  L21Metric fit_error;
  //L21Metric<FacetNormalMap, FacetAreaMap> fit_error;
  //L21Metric<PlaneProxy> fit_error;

  //boost::associative_property_map<std::map<face_descriptor, std::size_t> > seg_pmap;
  //boost::property_map<face_descriptor, std::size_t> &seg_pmap;

//member functions
public:
  /**
   * @pre @a polyhedron.is_pure_triangle()
   * @param mesh `CGAL Polyhedron` on which other functions operate.
   */
  VSA_segmentation(const Polyhedron &_mesh,
    VertexPointPmap _vertex_point_map,
    GeomTraits _traits)
    : mesh(_mesh),
    vertex_point_pmap(_vertex_point_map),
    traits(_traits),
    vector_functor(traits.construct_vector_3_object()),
    normal_functor(traits.construct_normal_3_object()),
    scale_functor(traits.construct_scaled_vector_3_object()),
    sum_functor(traits.construct_sum_of_vectors_3_object()),
    scalar_product_functor(traits.compute_scalar_product_3_object()),
    area_functor(traits.compute_squared_area_3_object()),
    normal_pmap(facet_normals),
    area_pmap(facet_areas),
    vertex_status_pmap(vertex_status_map),
    halfedge_status_pmap(halfedge_status_map),
    fit_error(traits, normal_pmap, area_pmap) {
    // CGAL_precondition(is_pure_triangle(mesh));
    // construct facet normal & area map
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      const halfedge_descriptor he = halfedge(f, mesh);
      const Point p1 = get(vertex_point_pmap, source(he, mesh));
      const Point p2 = get(vertex_point_pmap, target(he, mesh));
      const Point p3 = get(vertex_point_pmap, target(next(he, mesh), mesh));
      //const Point center = centroid_functor(p1, p2, p3);
      Vector normal = normal_functor(p1, p2, p3);
      normal = scale_functor(normal,
        FT(1.0 / std::sqrt(CGAL::to_double(normal.squared_length()))));
      facet_normals.insert(std::pair<face_descriptor, Vector>(f, normal));

      FT area(std::sqrt(CGAL::to_double(area_functor(p1, p2, p3))));
      facet_areas.insert(std::pair<face_descriptor, FT>(f, area));
    }

    // initialize all vertex anchor status
    BOOST_FOREACH(vertex_descriptor v, vertices(mesh)) {
      vertex_status_map.insert(std::pair<vertex_descriptor, int>(v, static_cast<int>(NO_ANCHOR)));
    }

    // tag all halfedge off proxy border
    BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh)) {
      halfedge_status_map.insert(std::pair<halfedge_descriptor, int>(h, static_cast<int>(OFF_BORDER)));
    }
  }

  template<typename FacetSegmentMap>
  void partition(const std::size_t number_of_segments, const std::size_t number_of_iterations, FacetSegmentMap &seg_pmap) {
    random_seed(number_of_segments);
    //random_seed(number_of_segments / 2);
    for (std::size_t i = 0; i < number_of_iterations; ++i) {
      flooding(seg_pmap);
      fitting(seg_pmap);
    }
  }

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

  // extract the approximated mesh from a partition
  template<typename FacetSegmentMap>
  void extract_mesh(FacetSegmentMap &seg_pmap, std::vector<int> &tris) {
    compute_proxy_area(seg_pmap);
    compute_proxy_center(seg_pmap);

    find_anchors(seg_pmap);
    find_edges(seg_pmap);
    add_anchors(seg_pmap);

    pseudo_CDT(seg_pmap, tris);
  }

  std::vector<Anchor> collect_anchors() {
    return anchors;
  }

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
  void random_seed(const std::size_t initial_px) {
    proxies.clear();
    std::size_t number_of_faces = num_faces(mesh);
    std::size_t interval = number_of_faces / initial_px;
    face_iterator fitr, fend;
    std::size_t index = 0;
    for (boost::tie(fitr, fend) = faces(mesh);
      (fitr != fend) && (proxies.size() < initial_px);
      ++fitr, ++index) {
      if (index % interval == 0) {
        // a proxy is created
        // PlaneProxy(face_descriptor)
        PlaneProxy px;
        px.normal = normal_pmap[*fitr];
        px.seed = *fitr;
        proxies.push_back(px);
      }
    }
    std::cerr << initial_px << ' ' << proxies.size() << std::endl;
  }

  template<typename FacetSegmentMap>
  void flooding(FacetSegmentMap &seg_pmap) {
    typedef std::multiset<FacetToIntegrate> CandidateSet;

    BOOST_FOREACH(face_descriptor f, faces(mesh))
      seg_pmap[f] = CGAL_NOT_TAGGED_ID;

    const std::size_t num_proxies = proxies.size();
    CandidateSet facet_candidates;
    for (std::size_t i = 0; i < num_proxies; ++i) {
      face_descriptor f = proxies[i].seed;
      seg_pmap[f] = i;

      BOOST_FOREACH(face_descriptor fadj, faces_around_face(halfedge(f, mesh), mesh)) {
        if (fadj != boost::graph_traits<Polyhedron>::null_face()) {
          FacetToIntegrate cand;
          cand.f = fadj;
          cand.fit_error = fit_error(fadj, proxies[i]);
          cand.i = i;
          facet_candidates.insert(cand);
        }
      }
    }

    while (!facet_candidates.empty()) {
      const FacetToIntegrate &c = *(facet_candidates.begin());
      facet_candidates.erase(facet_candidates.begin());
      if (seg_pmap[c.f] == CGAL_NOT_TAGGED_ID) {
        seg_pmap[c.f] = c.i;
        BOOST_FOREACH(face_descriptor fadj, faces_around_face(halfedge(c.f, mesh), mesh)) {
          if (fadj != boost::graph_traits<Polyhedron>::null_face()
            && seg_pmap[fadj] == CGAL_NOT_TAGGED_ID) {
            FacetToIntegrate cand;
            cand.f = fadj;
            cand.fit_error = fit_error(fadj, proxies[c.i]);
            cand.i = c.i;
            facet_candidates.insert(cand);
          }
        }
      }
    }
  }

  template <typename FacetSegmentMap>
  void fitting(FacetSegmentMap &seg_pmap) {
    // update normal
    std::vector<Vector> px_normals(proxies.size(), CGAL::NULL_VECTOR);
    std::vector<FT> px_areas(proxies.size(), FT(0));
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      std::size_t px_idx = seg_pmap[f];
      px_normals[px_idx] = sum_functor(px_normals[px_idx],
        scale_functor(normal_pmap[f], area_pmap[f]));
      px_areas[px_idx] += area_pmap[f];
    }
    
    for (std::size_t i = 0; i < proxies.size(); ++i) {
      Vector norm = scale_functor(px_normals[i], FT(1.0 / CGAL::to_double(px_areas[i]))); // redundant
      norm = scale_functor(norm, FT(1.0 / std::sqrt(CGAL::to_double(norm.squared_length()))));
      proxies[i].normal = norm;
    }

    // update seed
    std::vector<std::size_t> facet_px_idx;
    facet_px_idx.reserve(num_faces(mesh));
    std::vector<FT> facet_px_err;
    facet_px_err.reserve(num_faces(mesh));
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      std::size_t px_idx = seg_pmap[f];
      facet_px_idx.push_back(px_idx);
      facet_px_err.push_back(fit_error(f, proxies[px_idx]));
    }
    FT max_facet_error = facet_px_err.front();
    for (std::size_t i = 0; i < facet_px_err.size(); ++i) {
      if (max_facet_error < facet_px_err[i])
        max_facet_error = facet_px_err[i];
    }
    std::vector<FT> distance_min(proxies.size(), max_facet_error);
    std::size_t fidx = 0;
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      std::size_t px_idx = facet_px_idx[fidx];
      FT err = facet_px_err[fidx];
      if (err < distance_min[px_idx]) {
        proxies[px_idx].seed = f;
        distance_min[px_idx] = err;
      }
      ++fidx;
    }
  }

  // insert proxy at the facet with the maximum fitting error in the proxy with maximum error
  template <typename FacetSegmentMap>
  void insert_proxy(FacetSegmentMap &seg_pmap) {
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
    PlaneProxy new_px;
    new_px.normal = normal_pmap[max_facet[max_px_idx]];
    new_px.seed = max_facet[max_px_idx];
    proxies.push_back(new_px);
  }

  template<typename FacetSegmentMap>
  void find_anchors(FacetSegmentMap &seg_pmap) {
    anchors.clear();

    BOOST_FOREACH(vertex_descriptor vtx, vertices(mesh)) {
      std::set<std::size_t> px_set;
      std::size_t border_count = 0;

      BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(vtx, mesh)) {
        if (CGAL::is_border_edge(h, mesh)) {
          ++border_count;
          if (!CGAL::is_border(h, mesh))
            px_set.insert(seg_pmap[face(h, mesh)]);
        }
        else if (seg_pmap[face(h, mesh)] != seg_pmap[face(opposite(h, mesh), mesh)]) {
          ++border_count;
          px_set.insert(seg_pmap[face(h, mesh)]);
        }
      }
      if (border_count >= 3)
        attach_anchor(vtx, px_set);
    }
  }

  template<typename FacetSegmentMap>
  void find_edges(FacetSegmentMap &seg_pmap) {
    // tag all proxy border halfedges as candidate
    tag_halfedges_status(seg_pmap);

    // collect candidate halfedges in a set
    std::set<halfedge_descriptor> he_candidates;
    BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh)) {
      if (halfedge_status_pmap[h] == static_cast<int>(CANDIDATE))
        he_candidates.insert(h);
    }

    // pick up one candidate halfedge each time and traverse the connected border
    borders.clear();
    while (!he_candidates.empty()) {
      halfedge_descriptor he_start = *he_candidates.begin();
      walk_to_first_anchor(he_start, seg_pmap);
      if (!is_anchor_attached(he_start)) {
        // no anchor in this connected border, make a new anchor
        std::set<std::size_t> px_set;
        px_set.insert(seg_pmap[face(he_start, mesh)]);
        halfedge_descriptor he_oppo = opposite(he_start, mesh);
        if (!CGAL::is_border(he_oppo, mesh))
          px_set.insert(seg_pmap[face(he_oppo, mesh)]);
        attach_anchor(he_start, px_set);
      }

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

  // add anchors to the borders with only 2 anchors
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

      std::set<std::size_t> px_set;
      halfedge_descriptor he_oppo = opposite(he_max, mesh);
      px_set.insert(seg_pmap[face(he_max, mesh)]);
      if (!CGAL::is_border(he_oppo, mesh))
        px_set.insert(seg_pmap[face(he_oppo, mesh)]);
      attach_anchor(he_max, px_set);

      // increase border anchors by one
      bitr->num_anchors++;
    }
  }

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

    // mapping the Polyhedron mesh into a SubGraph
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
        if (f != boost::graph_traits<Polyhedron>::null_face())
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
      VertexIndex1Map &local_vanchor_map = get(boost::vertex_index1, glocal);
      VertexIndex2Map &local_vtag_map = get(boost::vertex_index2, glocal);
      EdgeWeightMap &local_eweight_map = get(boost::edge_weight, glocal);

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

  template<typename FacetSegmentMap>
  void compute_proxy_center(const FacetSegmentMap &seg_pmap) {
    proxies_center.clear();
    proxies_center.swap(std::vector<Point>(proxies.size()));
    
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

  template<typename FacetSegmentMap>
  void compute_proxy_area(const FacetSegmentMap &seg_pmap) {
    std::vector<FT> areas(proxies.size(), FT(0));
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      areas[seg_pmap[f]] += area_pmap[f];
    }
    proxies_area.swap(areas);
  }

  template<typename FacetSegmentMap>
  void tag_halfedges_status(FacetSegmentMap &seg_pmap) {
    BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh)) {
      halfedge_status_pmap[h] = static_cast<int>(OFF_BORDER);
      if (!CGAL::is_border(h, mesh)
        && (CGAL::is_border(opposite(h, mesh), mesh)
          || seg_pmap[face(h, mesh)] != seg_pmap[face(opposite(h, mesh), mesh)])) {
        halfedge_status_pmap[h] = static_cast<int>(CANDIDATE);
      }
    }
  }

  // walk along the proxy border to the first halfedge pointing to a vertex associated with an anchor
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

  // starting at a border halfedge pointing to a vertex associated with an anchor
  // walk along the proxy border to the next anchor, which may be itself
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

  // walk to next border halfedge, the input halfedge must be a border halfedge
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

  // subdivide a chord in range [chord_begin, chord_end)
  // returns the number of anchors of the chord apart from the first one
  template<typename FacetSegmentMap>
  std::size_t subdivide_chord(
    const ChordVectorIterator &chord_begin,
    const ChordVectorIterator &chord_end,
    const FacetSegmentMap &seg_pmap,
    const FT thre = FT(0.2)) {
    const std::size_t chord_size = std::distance(chord_begin, chord_end);
    if (chord_size < 2)
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
      std::set<std::size_t> px_set;
      px_set.insert(px_left);
      px_set.insert(px_right);
      attach_anchor(*he_max, px_set);

      std::size_t num0 = subdivide_chord(chord_begin, he_max + 1, seg_pmap);
      std::size_t num1 = subdivide_chord(he_max + 1, chord_end, seg_pmap);

      return num0 + num1;
    }

    return 1;
  }

  // check if halfedge target vertex is attached with an anchor
  bool is_anchor_attached(const halfedge_descriptor &he) {
    return is_anchor_attached(target(he, mesh), vertex_status_pmap);
  }

  // check if a vertex is attached with an anchor
  // suppose anchor index starts with 0
  template<typename VertexAnchorIndexMap>
  bool is_anchor_attached(
    const typename boost::property_traits<VertexAnchorIndexMap>::key_type &v,
    const VertexAnchorIndexMap &vertex_anchor_map) {
    return vertex_anchor_map[v] >= 0;
  }

  // attach anchor to vertex
  void attach_anchor(const vertex_descriptor &vtx, const std::set<std::size_t> &px_set) {
    // construct an anchor from vertex and the incident proxies
    FT avgx(0), avgy(0), avgz(0), sum_area(0);
    const Point vtx_pt = vertex_point_pmap[vtx];
    for (std::set<std::size_t>::iterator pxitr = px_set.begin();
      pxitr != px_set.end(); ++pxitr) {
      std::size_t px_idx = *pxitr;
      Plane px_plane(proxies_center[px_idx], proxies[px_idx].normal);
      Point proj = px_plane.projection(vtx_pt);
      FT area = proxies_area[px_idx];
      avgx += proj.x() * area;
      avgy += proj.y() * area;
      avgz += proj.z() * area;
      sum_area += area;
    }
    Point pos = Point(avgx / sum_area, avgy / sum_area, avgz / sum_area);
    vertex_status_pmap[vtx] = static_cast<int>(anchors.size());
    anchors.push_back(Anchor(vtx, pos));
  }

  // attach anchor to the target vertex of the halfedge
  void attach_anchor(const halfedge_descriptor &he, const std::set<std::size_t> &px_set) {
    vertex_descriptor vtx = target(he, mesh);
    attach_anchor(vtx, px_set);
  }
}; // end class VSA_segmentation
} // end namespace internal
} // end namespace CGAL

#undef CGAL_NOT_TAGGED_ID

#endif // CGAL_VSA_SEGMENTATION_H
