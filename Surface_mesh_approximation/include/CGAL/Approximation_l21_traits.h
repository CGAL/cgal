#ifndef CGAL_APPROXIMATION_L21_TRAITS_H
#define CGAL_APPROXIMATION_L21_TRAITS_H

#include <CGAL/license/Surface_mesh_approximation.h>

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/unordered_map.hpp>

namespace CGAL {

/*!
 * \ingroup PkgTSMA
 * @brief Approximation traits for L21 metric.
 *
 * \cgalModels`Approximation_traits`
 *
 * @tparam TriangleMesh a triangle `FaceListGraph`
 * @tparam VertexPointMap a property map with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
    as key type, GeomTraits::Point_3 as value type
 * @tparam with_area_weighing set `true` to activate area weighing
 * @tparam GeomTraits geometric traits
 */
template <typename TriangleMesh,
  typename VertexPointMap
    = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  bool with_area_weighing = true,
  typename GeomTraits
    = typename TriangleMesh::Traits>
class Approximation_l21_traits {
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename GeomTraits::Compute_scalar_product_3 Compute_scalar_product_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef CGAL::internal::face_property_t<Vector_3> Face_normal_tag;
  typedef CGAL::internal::face_property_t<FT> Face_area_tag;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_normal_tag >::type Face_normal_map;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_area_tag >::type Face_area_map;

public:
  // type required by the `Approximation_traits` concept
  typedef typename GeomTraits::Vector_3 Proxy;

  // constructor
  Approximation_l21_traits(const TriangleMesh &tm, const VertexPointMap &vpmap) {
    GeomTraits traits;
    m_scalar_product_functor = traits.compute_scalar_product_3_object();
    m_sum_functor = traits.construct_sum_of_vectors_3_object();
    m_scale_functor = traits.construct_scaled_vector_3_object();
    
    m_fnmap = CGAL::internal::add_property(
      Face_normal_tag("VSA-face_normal"), const_cast<TriangleMesh &>(tm));
    m_famap = CGAL::internal::add_property(
      Face_area_tag("VSA-face_area"), const_cast<TriangleMesh &>(tm));
    // construct internal facet normal & area map
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = vpmap[source(he, tm)];
      const Point_3 &p1 = vpmap[target(he, tm)];
      const Point_3 &p2 = vpmap[target(next(he, tm), tm)];
      put(m_fnmap, f, CGAL::unit_normal(p0, p1, p2));
      put(m_famap, f, std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    }
  }

  // member function required by the `Approximation_traits` concept
  // It is a function that takes a facet and a proxy, returns the L21 error between them.
  FT compute_error(const face_descriptor &f, const Proxy &px) const {
    Vector_3 v = m_sum_functor(get(m_fnmap, f), m_scale_functor(px, FT(-1.0)));
    return get(m_famap, f) * m_scalar_product_functor(v, v);
  }

  // member function required by the `Approximation_traits` concept
  // It returns the proxy fitted from the facets from beg to end.
  template <typename FacetIterator>
  Proxy fit_proxy(const FacetIterator beg, const FacetIterator end) const {
    CGAL_assertion(beg != end);

    // fitting normal
    Vector_3 norm = CGAL::NULL_VECTOR;
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      norm = m_sum_functor(norm,
        m_scale_functor(get(m_fnmap, *fitr), get(m_famap, *fitr)));
    }
    norm = m_scale_functor(norm,
      FT(1.0 / std::sqrt(CGAL::to_double(norm.squared_length()))));

    return norm;
  }

private:
  Face_normal_map m_fnmap;
  Face_area_map m_famap;
  Construct_scaled_vector_3 m_scale_functor;
  Compute_scalar_product_3 m_scalar_product_functor;
  Construct_sum_of_vectors_3 m_sum_functor;
};

// specialization without area weighing
template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
class Approximation_l21_traits<TriangleMesh,
  VertexPointMap,
  false,
  GeomTraits> {
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename GeomTraits::Compute_scalar_product_3 Compute_scalar_product_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type Vertex_point_map;

  typedef CGAL::internal::face_property_t<Vector_3> Face_normal_tag;
  typedef CGAL::internal::face_property_t<FT> Face_area_tag;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_normal_tag >::type Face_normal_map;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Face_area_tag >::type Face_area_map;

public:
  // type required by the `Approximation_traits` concept
  typedef typename GeomTraits::Vector_3 Proxy;

  // constructor
  Approximation_l21_traits(const TriangleMesh &tm, const VertexPointMap &vpmap) {
    GeomTraits traits;
    m_scalar_product_functor = traits.compute_scalar_product_3_object();
    m_sum_functor = traits.construct_sum_of_vectors_3_object();
    m_scale_functor = traits.construct_scaled_vector_3_object();

    m_fnmap = CGAL::internal::add_property(
      Face_normal_tag("VSA-face_normal"), const_cast<TriangleMesh &>(tm));

    // construct internal facet normal & area map
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      const halfedge_descriptor he = halfedge(f, tm);
      const Point_3 &p0 = vpmap[source(he, tm)];
      const Point_3 &p1 = vpmap[target(he, tm)];
      const Point_3 &p2 = vpmap[target(next(he, tm), tm)];
      put(m_fnmap, f, CGAL::unit_normal(p0, p1, p2));
    }
  }

  // member function required by the `Approximation_traits` concept
  // It is a function that takes a facet and a proxy, returns the L21 error between them.
  FT compute_error(const face_descriptor &f, const Proxy &px) const {
    Vector_3 v = m_sum_functor(get(m_fnmap, f), m_scale_functor(px, FT(-1.0)));
    return m_scalar_product_functor(v, v);
  }

  // member function required by the `Approximation_traits` concept
  // It returns the proxy fitted from the facets from beg to end.
  template <typename FacetIterator>
  Proxy fit_proxy(const FacetIterator beg, const FacetIterator end) const {
    CGAL_assertion(beg != end);

    // fitting normal
    Vector_3 norm = CGAL::NULL_VECTOR;
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      norm = m_sum_functor(norm, get(m_fnmap, *fitr));
    }
    norm = m_scale_functor(norm,
      FT(1.0 / std::sqrt(CGAL::to_double(norm.squared_length()))));

    return norm;
  }

private:
  Face_normal_map m_fnmap;
  Construct_scaled_vector_3 m_scale_functor;
  Compute_scalar_product_3 m_scalar_product_functor;
  Construct_sum_of_vectors_3 m_sum_functor;
};

} // namespace CGAL

#endif // CGAL_APPROXIMATION_L21_TRAITS_H
