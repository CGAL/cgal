#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_TRAITS_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_TRAITS_H

#include <list>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <boost/graph/graph_traits.hpp>

namespace CGAL
{
/**
 * @brief Plane proxy class for the Variational Shape Approximation algorithm.
 *
 * It is simply a class containing few proxy parameters.
 * It is used as the default proxy for the `L21Metric` and `L2Metric`
 *
 * @tparam TriangleMesh a triangle `FaceGraph`
 * @tparam GeomTraits geometric traits
 */
template<typename TriangleMesh,
  typename GeomTraits = typename TriangleMesh::Traits>
class PlaneProxy
{
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Plane_3 Plane_3;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

public:
  // The proxy seed.
  face_descriptor seed;
  // The proxy normal used in the `L21Metric`.
  Vector_3 normal;
  // The fitting plane of the proxy used in the `L2Metric`.
  Plane_3 fit_plane;
};

/**
 * @brief L21 metric class for the Variational Shape Approximation algorithm.
 *
 * It is simply a functor that takes a facet and a proxy, returns the L21 error between them.
 *
 * @tparam TriangleMesh a triangle `FaceGraph`
 * @tparam FacetNormalMap a property map containing the facet normals,
           and `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type,
           GeomTraits::Vector_3 as value type
 * @tparam FacetAreaMap a property map containing the facet areas,
           and `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type,
           GeomTraits::FT as value type
 * @tparam GeomTraits geometric traits
 * @tparam PlaneProxy a model of `PlaneProxy`
 */
template<typename TriangleMesh,
  typename FacetNormalMap,
  typename FacetAreaMap,
  typename GeomTraits = typename TriangleMesh::Traits,
  typename PlaneProxy = CGAL::PlaneProxy<TriangleMesh, GeomTraits> >
class L21Metric
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename GeomTraits::Compute_scalar_product_3 Compute_scalar_product_3;
  typedef typename FacetAreaMap::key_type face_descriptor;

public:
  // The type define required by the `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // default constructor
  L21Metric() {
    GeomTraits traits;
    scalar_product_functor = traits.compute_scalar_product_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
  }

  // construct L21 metric functor from a facet normal map and a facet area map.
  L21Metric(const FacetNormalMap &normal_pmap, const FacetAreaMap &area_pmap)
    : normal_pmap(normal_pmap),
    area_pmap(area_pmap) {
    GeomTraits traits;
    scalar_product_functor = traits.compute_scalar_product_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
  }

  // returns L21 error of a facet f to a proxy px.
  FT operator()(const face_descriptor &f, const Proxy &px) const {
    Vector_3 v = sum_functor(normal_pmap[f], scale_functor(px.normal, FT(-1)));
    return area_pmap[f] * scalar_product_functor(v, v);
  }

private:
  const FacetNormalMap normal_pmap;
  const FacetAreaMap area_pmap;
  Construct_scaled_vector_3 scale_functor;
  Compute_scalar_product_3 scalar_product_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

/**
 * @brief L21 proxy fitting class for the Variational Shape Approximation algorithm.
 *
 * It is simply a functor that takes a range of facets, fitting the proxy parameters.
 *
 * @tparam TriangleMesh a triangle `FaceGraph`
 * @tparam FacetNormalMap a property map containing the facet normals,
           and `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type,
           GeomTraits::Vector_3 as value type
 * @tparam FacetAreaMap a property map containing the facet areas,
           and `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type,
           GeomTraits::FT as value type
 * @tparam GeomTraits geometric traits
 * @tparam PlaneProxy a model of `PlaneProxy`
 */
template<typename TriangleMesh,
  typename FacetNormalMap,
  typename FacetAreaMap,
  typename GeomTraits = typename TriangleMesh::Traits,
  typename PlaneProxy = CGAL::PlaneProxy<TriangleMesh, GeomTraits> >
class L21ProxyFitting
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;

public:
  // The type define required by the `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // default constructor
  L21ProxyFitting() {
    GeomTraits traits;
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
  }

  // construct L21 proxy fitting functor from a facet normal map and a facet area map.
  L21ProxyFitting(const FacetNormalMap &normal_pmap, const FacetAreaMap &area_pmap)
    : normal_pmap(normal_pmap), area_pmap(area_pmap) {
    GeomTraits traits;
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
  }

  // returns the proxy fitted from the facets from beg to end.
  template<typename FacetIterator>
  Proxy operator()(const FacetIterator beg, const FacetIterator end) const {
    CGAL_assertion(beg != end);

    // fitting normal
    Vector_3 norm = CGAL::NULL_VECTOR;
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      norm = sum_functor(norm,
        scale_functor(normal_pmap[*fitr], area_pmap[*fitr]));
    }
    norm = scale_functor(norm,
      FT(1.0 / std::sqrt(CGAL::to_double(norm.squared_length()))));

    // construct proxy
    Proxy px;
    px.normal = norm;

    return px;
  }

private:
  const FacetNormalMap normal_pmap;
  const FacetAreaMap area_pmap;
  Construct_scaled_vector_3 scale_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

/**
 * @brief Area weighted plane fitting class for the Variational Shape Approximation algorithm.
 *
 * It is simply a functor class that fits a plane from a range of facets.
 *
 * @tparam TriangleMesh a triangle `FaceGraph`
 * @tparam VertexPointMap a property map containing the vertex points,
           and `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type,
           GeomTraits::Point_3 as value type
 * @tparam GeomTraits geometric traits
 */
template<typename TriangleMesh,
  typename VertexPointMap
    = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  typename GeomTraits = typename TriangleMesh::Traits>
class PlaneFitting
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Plane_3 Plane_3;
  typedef typename GeomTraits::Construct_vector_3 Construct_vector_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

public:
  // construct plane fitting functor from a triangle mesh and its vertex point map.
  PlaneFitting(const TriangleMesh &_mesh, const VertexPointMap &_point_pmap)
    : mesh(_mesh), point_pmap(_point_pmap) {
    GeomTraits traits;
    vector_functor = traits.construct_vector_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
  }

  // construct plane fitting functor from a triangle mesh.
  PlaneFitting(const TriangleMesh &_mesh)
    : mesh(_mesh),
    point_pmap(get(boost::vertex_point, const_cast<TriangleMesh &>(_mesh))) {
    GeomTraits traits;
    vector_functor = traits.construct_vector_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
  }

  // fitting a plane from facets in range beg to end.
  template<typename FacetIterator>
  Plane_3 operator()(const FacetIterator &beg, const FacetIterator &end) const {
    CGAL_assertion(beg != end);
    // area average normal and centroid
    Vector_3 norm = CGAL::NULL_VECTOR;
    Vector_3 cent = CGAL::NULL_VECTOR;
    FT sum_area(0);
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      const halfedge_descriptor he = halfedge(*fitr, mesh);
      const Point_3 p0 = point_pmap[source(he, mesh)];
      const Point_3 p1 = point_pmap[target(he, mesh)];
      const Point_3 p2 = point_pmap[target(next(he, mesh), mesh)];

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

private:
  const TriangleMesh &mesh;
  const VertexPointMap point_pmap;
  Construct_vector_3 vector_functor;
  Construct_scaled_vector_3 scale_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

/**
 * @brief L2 metric class for the Variational Shape Approximation algorithm.
 *
 * It is simply a functor that takes a facet and a proxy, returns the L2 error between them.
 *
 * @tparam TriangleMesh a triangle `FaceGraph`
 * @tparam FacetAreaMap a property map containing the facet areas,
           and `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type,
           GeomTraits::FT as value type
 * @tparam VertexPointMap a property map containing the vertex points,
           and `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type,
           GeomTraits::Point_3 as value type
 * @tparam GeomTraits geometric traits
 * @tparam PlaneProxy a model of `PlaneProxy`
 */
template<typename TriangleMesh,
  typename FacetAreaMap,
  typename VertexPointMap
    = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  typename GeomTraits = typename TriangleMesh::Traits,
  typename PlaneProxy = CGAL::PlaneProxy<TriangleMesh, GeomTraits> >
class L2Metric
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

public:
  // The type define required by the `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // default constructor
  L2Metric() : mesh(nullptr) {}

  // construct L2 metric functor from a triangle mesh, a facet area map and the vertex point map.
  L2Metric(const TriangleMesh &_mesh,
    const FacetAreaMap &_area_pmap,
    const VertexPointMap &_point_pmap)
    : mesh(&_mesh), area_pmap(_area_pmap), point_pmap(_point_pmap) {}


  // construct L2 metric functor from a triangle mesh and a facet area map.
  L2Metric(const TriangleMesh &_mesh,
    const FacetAreaMap &_area_pmap)
    : mesh(&_mesh), area_pmap(_area_pmap),
    point_pmap(get(boost::vertex_point, const_cast<TriangleMesh &>(_mesh))) {}

  // returns the l2 fitting error of a facet f to proxy px.
  FT operator()(const face_descriptor &f, const PlaneProxy &px) const {
    halfedge_descriptor he = halfedge(f, *mesh);
    const Point_3 &p0 = point_pmap[source(he, *mesh)];
    const Point_3 &p1 = point_pmap[target(he, *mesh)];
    const Point_3 &p2 = point_pmap[target(next(he, *mesh), *mesh)];
    FT sq_d0 = CGAL::squared_distance(p0, px.fit_plane);
    FT sq_d1 = CGAL::squared_distance(p1, px.fit_plane);
    FT sq_d2 = CGAL::squared_distance(p2, px.fit_plane);
    FT d0(std::sqrt(CGAL::to_double(sq_d0)));
    FT d1(std::sqrt(CGAL::to_double(sq_d1)));
    FT d2(std::sqrt(CGAL::to_double(sq_d2)));

    return (sq_d0 + sq_d1 + sq_d2 + d0 * d1 + d1 * d2 + d2 * d0) * area_pmap[f] / FT(6);
  }

private:
  const FacetAreaMap area_pmap;
  const VertexPointMap point_pmap;
  const TriangleMesh *mesh;
};

/**
 * @brief L2 proxy fitting class for the Variational Shape Approximation algorithm.
 *
 * It is simply a functor that takes a range of facets, fitting the L2 proxy parameters.
 * It uses the PCA algorithm to fit the proxy parameters.
 *
 * @tparam TriangleMesh a triangle `FaceGraph`
 * @tparam VertexPointMap a property map containing the vertex points,
           and `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type,
           GeomTraits::Point_3 as value type
 * @tparam GeomTraits geometric traits
 * @tparam PlaneProxy a model of `PlaneProxy`
 */
template<typename TriangleMesh,
  typename VertexPointMap
    = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  typename GeomTraits = typename TriangleMesh::Traits,
  typename PlaneProxy = CGAL::PlaneProxy<TriangleMesh, GeomTraits> >
class L2ProxyFitting
{
private:
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Triangle_3 Triangle_3;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

public:
  // The type define required by the `ErrorMetric` concept
  typedef PlaneProxy Proxy;

  // default constructor
  L2ProxyFitting() : mesh(nullptr) {}

  // construct L2 proxy fitting functor from a triangle mesh and the vertex point map.
  L2ProxyFitting(const TriangleMesh &_mesh, const VertexPointMap &_point_pmap)
    : mesh(&_mesh), point_pmap(_point_pmap) {}

  // construct L2 proxy fitting functor from a triangle mesh.
  L2ProxyFitting(const TriangleMesh &_mesh)
    : mesh(&_mesh),
    point_pmap(get(boost::vertex_point, const_cast<TriangleMesh &>(_mesh))) {}

  // returns the proxy fitted from a range of facets.
  template<typename FacetIterator>
  Proxy operator()(const FacetIterator beg, const FacetIterator end) const {
    CGAL_assertion(beg != end);

    std::list<Triangle_3> tris;
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      halfedge_descriptor he = halfedge(*fitr, *mesh);
      const Point_3 &p0 = point_pmap[source(he, *mesh)];
      const Point_3 &p1 = point_pmap[target(he, *mesh)];
      const Point_3 &p2 = point_pmap[target(next(he, *mesh), *mesh)];
      tris.push_back(Triangle_3(p0, p1, p2));
    }

    // construct and fit proxy plane
    Proxy px;
    CGAL::linear_least_squares_fitting_3(
      tris.begin(),
      tris.end(),
      px.fit_plane,
      CGAL::Dimension_tag<2>());

    return px;
  }

private:
  const TriangleMesh *mesh;
  const VertexPointMap point_pmap;
};

/**
 * @brief PCA plane fitting class.
 *
 * It is simply a functor class that uses the PCA algorithm to fit a plane from a range of facets.
 *
 * @tparam TriangleMesh a triangle `FaceGraph`
 * @tparam VertexPointMap a property map containing the vertex points,
           and `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type,
           GeomTraits::Point_3 as value type
 * @tparam GeomTraits geometric traits
 */
template<typename TriangleMesh,
  typename VertexPointMap
    = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  typename GeomTraits = typename TriangleMesh::Traits>
class PCAPlaneFitting
{
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Plane_3 Plane_3;
  typedef typename GeomTraits::Triangle_3 Triangle_3;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

public:
  // construct PCA plane fitting functor from a triangle mesh and the vertex point map.
  PCAPlaneFitting(const TriangleMesh &_mesh,  const VertexPointMap &_point_pmap)
    : mesh(_mesh), point_pmap(_point_pmap) {}

  // construct PCA plane fitting functor from a triangle mesh.
  PCAPlaneFitting(const TriangleMesh &_mesh)
    : mesh(_mesh),
    point_pmap(get(boost::vertex_point, const_cast<TriangleMesh &>(_mesh))) {}

  // returns a plane fitted from facets in range beg to end by the PCA algorithm.
  template<typename FacetIterator>
  Plane_3 operator()(const FacetIterator beg, const FacetIterator end) const {
    CGAL_assertion(beg != end);

    std::list<Triangle_3> tris;
    for (FacetIterator fitr = beg; fitr != end; ++fitr) {
      halfedge_descriptor he = halfedge(*fitr, mesh);
      const Point_3 &p0 = point_pmap[source(he, mesh)];
      const Point_3 &p1 = point_pmap[target(he, mesh)];
      const Point_3 &p2 = point_pmap[target(next(he, mesh), mesh)];
      tris.push_back(Triangle_3(p0, p1, p2));
    }

    // construct and fit proxy plane
    Plane_3 fit_plane;
    CGAL::linear_least_squares_fitting_3(
      tris.begin(),
      tris.end(),
      fit_plane,
      CGAL::Dimension_tag<2>());
    
    return fit_plane;
  }

private:
  const TriangleMesh &mesh;
  const VertexPointMap point_pmap;
};
} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_TRAITS_H
