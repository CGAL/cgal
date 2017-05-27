#ifndef CGAL_VSA_SEGMENTATION_H
#define CGAL_VSA_SEGMENTATION_H

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
    std::vector<FT> distance_min(proxies.size(), FT(100000.0));
    //std::vector<FT> distance_min(proxies.size(), FT::max());
    for (boost::tie(fitr, fend) = faces(mesh); fitr != fend; ++fitr) {
      std::size_t px_idx = seg_pmap[*fitr];
      FT err = fit_error(*fitr, proxies[px_idx]);
      if (err < distance_min[px_idx]) {
        proxies[px_idx].seed = *fitr;
        distance_min[px_idx] = err;
      }
    }
  }
};
}
}

#undef CGAL_NOT_TAGGED_ID

#endif // CGAL_VSA_SEGMENTATION_H
