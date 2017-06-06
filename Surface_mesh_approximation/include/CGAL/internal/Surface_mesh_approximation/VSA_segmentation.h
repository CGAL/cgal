#ifndef CGAL_VSA_SEGMENTATION_H
#define CGAL_VSA_SEGMENTATION_H

#include <CGAL/boost/graph/helpers.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>

#include <vector>
#include <cmath>
#include <map>
#include <set>
#include <iostream>

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

  //typedef boost::associative_property_map<std::map<face_descriptor, std::size_t> > FacetSegmentMap;
  typedef boost::associative_property_map<std::map<face_descriptor, Vector> > FacetNormalMap;
  typedef boost::associative_property_map<std::map<face_descriptor, FT> > FacetAreaMap;
  typedef boost::associative_property_map<std::map<vertex_descriptor, int> > VertexStatusPMap;
  typedef boost::associative_property_map<std::map<halfedge_descriptor, int> > HalfedgeStatusPMap;

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
    //Point center;
    Vector normal;
    face_descriptor seed;
  };

  struct FacetToIntegrate {
    face_descriptor f;
    std::size_t i;
    FT fit_error;
  };

  struct CompFacet {
    bool operator()(const FacetToIntegrate &f1, const FacetToIntegrate &f2) const {
      return f1.fit_error < f2.fit_error;
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
    // construct an anchor from vertex and the incident proxies
    Anchor(const vertex_descriptor &vtx_, const Point &vtx_pt_, const std::set<std::size_t> &px_set_)
    : vtx(vtx_) {
      FT avgx(0), avgy(0), avgz(0), sum_area(0);
      for (std::set<std::size_t>::iterator pxitr = px_set_.begin();
        pxitr != px_set_.end(); ++pxitr) {
        std::size_t px_idx = *pxitr;
        // TODO: Plane px_plane(proxies[px_idx].center, proxies[px_idx].normal);
        Plane px_plane;
        Point proj = px_plane.projection(vtx_pt_);
        // TODO: FT area = proxies[px_idx];
        FT area = FT(0.0);
        avgx += proj.x() * area;
        avgy += proj.y() * area;
        avgz += proj.z() * area;
        sum_area += area;
      }
      pos = Point(avgx / sum_area, avgy / sum_area, avgz / sum_area);
    }

    vertex_descriptor vtx; // The associated vertex
    Point pos; // The position of the anchor
  };

// member variables
private:
  const Polyhedron &mesh;
  const VertexPointPmap &vertex_point_map;
  GeomTraits traits;
  Construct_normal_3 normal_functor;
  Construct_scaled_vector_3 scale_functor;
  Construct_sum_of_vectors_3 sum_functor;
  Compute_scalar_product_3 scalar_product_functor;
  Compute_squared_area_3 area_functor;

  std::vector<PlaneProxy> proxies;
  // the lifetime of the container must encompass the use of the adaptor
  std::map<face_descriptor, Vector> facet_normals;
  FacetNormalMap normal_pmap;
  std::map<face_descriptor, FT> facet_areas;
  FacetAreaMap area_pmap;

  std::map<vertex_descriptor, int> vertex_status_map;
  VertexStatusPMap vertex_status_pmap;
  std::map<halfedge_descriptor, int> halfedge_status_map;
  HalfedgeStatusPMap halfedge_status_pmap;

  std::vector<Anchor> anchors;

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
    vertex_point_map(_vertex_point_map),
    traits(_traits),
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
    // construct facet normal map
    face_iterator fitr, fend;
    for (boost::tie(fitr, fend) = faces(mesh); fitr != fend; ++fitr) {
      const Point p1 = get(vertex_point_map, target(halfedge(*fitr, mesh), mesh));
      const Point p2 = get(vertex_point_map, target(next(halfedge(*fitr, mesh), mesh), mesh));
      const Point p3 = get(vertex_point_map, target(prev(halfedge(*fitr, mesh), mesh), mesh));
      //const Point center = centroid_functor(p1, p2, p3);
      Vector normal = normal_functor(p1, p2, p3);
      normal = scale_functor(normal,
        FT(1.0 / std::sqrt(to_double(normal.squared_length()))));

      facet_normals.insert(std::pair<face_descriptor, Vector>(*fitr, normal));
    }

    // construct facet area map
    for (boost::tie(fitr, fend) = faces(mesh); fitr != fend; ++fitr) {
      const Point p1 = get(vertex_point_map, target(halfedge(*fitr, mesh), mesh));
      const Point p2 = get(vertex_point_map, target(next(halfedge(*fitr, mesh), mesh), mesh));
      const Point p3 = get(vertex_point_map, target(prev(halfedge(*fitr, mesh), mesh), mesh));
      FT area(std::sqrt(to_double(area_functor(p2, p1, p3))));
      facet_areas.insert(std::pair<face_descriptor, FT>(*fitr, area));
    }

    // tag all vertex without anchor
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
  void extract_mesh(FacetSegmentMap &seg_pmap) {
    find_anchors(seg_pmap);
    tag_halfedges_status(seg_pmap);
    find_edges(seg_pmap);

    pseudo_CDT();
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
    face_iterator fitr, fend;
    for (boost::tie(fitr, fend) = faces(mesh); fitr != fend; ++fitr)
      seg_pmap[*fitr] = CGAL_NOT_TAGGED_ID;

    typedef std::multiset<FacetToIntegrate, CompFacet> CandidateSet;

    const std::size_t num_proxies = proxies.size();
    CandidateSet facet_candidates;
    for (std::size_t i = 0; i < num_proxies; ++i) {
      face_descriptor f = proxies[i].seed;
      seg_pmap[f] = i;

      Halfedge_around_face_circulator<Polyhedron> facet_circulator(halfedge(f, mesh), mesh), done(facet_circulator);
      do {
        if (face(opposite(*facet_circulator, mesh), mesh) != boost::graph_traits<Polyhedron>::null_face()) {
          FacetToIntegrate cand;
          cand.f = face(opposite(*facet_circulator, mesh), mesh);
          cand.fit_error = fit_error(cand.f, proxies[i]);
          cand.i = i;
          facet_candidates.insert(cand);
        }
      } while (++facet_circulator != done);
    }

    for (CandidateSet::iterator citr = facet_candidates.begin(); citr != facet_candidates.end();) {
      if (seg_pmap[citr->f] == CGAL_NOT_TAGGED_ID) {
        seg_pmap[citr->f] = citr->i;
        Halfedge_around_face_circulator<Polyhedron> facet_circulator(halfedge(citr->f, mesh), mesh), done(facet_circulator);
        do {
          face_descriptor oppo_facet = face(opposite(*facet_circulator, mesh), mesh);
          if (oppo_facet != boost::graph_traits<Polyhedron>::null_face()
            && seg_pmap[oppo_facet] == CGAL_NOT_TAGGED_ID) {
            FacetToIntegrate cand;
            cand.f = oppo_facet;
            cand.fit_error = fit_error(oppo_facet, proxies[citr->i]);
            cand.i = citr->i;
            facet_candidates.insert(cand);
          }
        } while (++facet_circulator != done);
      }
      facet_candidates.erase(citr);
      citr = facet_candidates.begin();
    }
  }

  template <typename FacetSegmentMap>
  void fitting(FacetSegmentMap &seg_pmap) {
    // update normal
    std::vector<Vector> px_normals(proxies.size(), CGAL::NULL_VECTOR);
    std::vector<FT> px_areas(proxies.size(), FT());
    face_iterator fitr, fend;
    for (boost::tie(fitr, fend) = faces(mesh); fitr != fend; ++fitr) {
      std::size_t px_idx = seg_pmap[*fitr];
      px_normals[px_idx] = sum_functor(px_normals[px_idx],
        scale_functor(normal_pmap[*fitr], area_pmap[*fitr]));
      px_areas[px_idx] += area_pmap[*fitr];
    }
    for (std::size_t i = 0; i < proxies.size(); ++i) {
      Vector norm = scale_functor(px_normals[i], FT(1.0 / to_double(px_areas[i]))); // redundant
      norm = scale_functor(norm, FT(1.0 / std::sqrt(to_double(norm.squared_length()))));
      proxies[i].normal = norm;
    }

    // update seed
    std::vector<std::size_t> facet_px_idx;
    facet_px_idx.reserve(num_faces(mesh));
    std::vector<FT> facet_px_err;
    facet_px_err.reserve(num_faces(mesh));
    for (boost::tie(fitr, fend) = faces(mesh); fitr != fend; ++fitr) {
      std::size_t px_idx = seg_pmap[*fitr];
      facet_px_idx.push_back(px_idx);
      facet_px_err.push_back(fit_error(*fitr, proxies[px_idx]));
    }
    FT max_facet_error = facet_px_err.front();
    for (std::size_t i = 0; i < facet_px_err.size(); ++i) {
      if (max_facet_error < facet_px_err[i])
        max_facet_error = facet_px_err[i];
    }
    std::vector<FT> distance_min(proxies.size(), max_facet_error);
    std::size_t fidx = 0;
    for (boost::tie(fitr, fend) = faces(mesh); fitr != fend; ++fitr, ++fidx) {
      std::size_t px_idx = facet_px_idx[fidx];
      FT err = facet_px_err[fidx];
      if (err < distance_min[px_idx]) {
        proxies[px_idx].seed = *fitr;
        distance_min[px_idx] = err;
      }
    }
  }

  // insert proxy at the facet with the maximum fitting error in the proxy with maximum error
  template <typename FacetSegmentMap>
  void insert_proxy(FacetSegmentMap &seg_pmap) {
    std::vector<FT> px_error(proxies.size(), FT(0.0));
    std::vector<FT> max_facet_error(proxies.size(), FT(0.0));
    std::vector<face_descriptor> max_facet(proxies.size());
    face_iterator fitr, fend;
    for (boost::tie(fitr, fend) = faces(mesh); fitr != fend; ++fitr) {
      std::size_t px_idx = seg_pmap[*fitr];
      FT err = fit_error(*fitr, proxies[px_idx]);
      px_error[px_idx] += err;

      if (err > max_facet_error[px_idx]) {
        max_facet_error[px_idx] = err;
        max_facet[px_idx] = *fitr;
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

    vertex_iterator vitr, vend;
    for (boost::tie(vitr, vend) = vertices(mesh); vitr != vend; ++vitr) {
      std::set<std::size_t> px_set;
      std::size_t border_count = 0;

      BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(*vitr, mesh)) {
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
      if (border_count >= 3) {
        // make an anchor and attach it to the vertex
        vertex_status_pmap[*vitr] = static_cast<int>(anchors.size());
        anchors.push_back(Anchor(*vitr, vertex_point_map[*vitr], px_set));
      }
    }
  }

  template<typename FacetSegmentMap>
  void tag_halfedges_status(FacetSegmentMap &seg_pmap) {
    BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh)) {
      if (!CGAL::is_border(h, mesh)
        && (CGAL::is_border(opposite(h, mesh), mesh)
          || seg_pmap[face(h, mesh)] != seg_pmap[face(opposite(h, mesh), mesh)])) {
        halfedge_status_pmap[h] = static_cast<int>(ON_BORDER);
      }
    }
  }

  template<typename FacetSegmentMap>
  void find_edges(FacetSegmentMap &seg_pmap) {

  }

  void pseudo_CDT() {}
};
}
}

#undef CGAL_NOT_TAGGED_ID

#endif // CGAL_VSA_SEGMENTATION_H
