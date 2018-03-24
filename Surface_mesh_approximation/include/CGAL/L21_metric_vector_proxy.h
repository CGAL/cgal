#ifndef CGAL_L21_METRIC_PLANE_PROXY_H
#define CGAL_L21_METRIC_PLANE_PROXY_H

#include <CGAL/license/Surface_mesh_approximation.h>

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/unordered_map.hpp>

namespace CGAL {
namespace VSA {

/// \ingroup PkgTSMA
/// @brief Approximation L21 metric of vector proxy.
///
/// \cgalModels `ErrorMetricProxy`
///
/// @tparam TriangleMesh a triangle `FaceListGraph`
/// @tparam VertexPointMap a property map with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///    as key type, GeomTraits::Point_3 as value type
/// @tparam GeomTraits geometric traits
template <typename TriangleMesh,
  typename VertexPointMap
    = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  typename GeomTraits
    = typename TriangleMesh::Traits>
class L21_metric_vector_proxy {
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
  /// \name Types
  /// @{
  typedef typename GeomTraits::Vector_3 Proxy;
  /// @}

  /// \name Constructor
  /// @{
  /*!
   * @brief Constructor
   * @param tm triangle mesh
   * @param vpmap vertex point map
   */
  L21_metric_vector_proxy(const TriangleMesh &tm, const VertexPointMap &vpmap) {
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
  /// @}

  /// \name Operations
  /*!
   * @brief Computes the L2,1 error from a facet to a proxy. 
   * @param f face_descriptor of a face
   * @param px proxy
   * @return computed error
   */
  FT compute_error(const face_descriptor &f, const Proxy &px) const {
    Vector_3 v = m_sum_functor(get(m_fnmap, f), m_scale_functor(px, FT(-1.0)));
    return get(m_famap, f) * m_scalar_product_functor(v, v);
  }

  /*!
   * @brief Fits a proxy to a face range.
   * @param beg face range begin
   * @param end face range end
   * @return fitted proxy
   */
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
  /// @}

private:
  Face_normal_map m_fnmap;
  Face_area_map m_famap;
  Construct_scaled_vector_3 m_scale_functor;
  Compute_scalar_product_3 m_scalar_product_functor;
  Construct_sum_of_vectors_3 m_sum_functor;
};

} // namespace VSA
} // namespace CGAL

#endif // CGAL_L21_METRIC_PLANE_PROXY_H
