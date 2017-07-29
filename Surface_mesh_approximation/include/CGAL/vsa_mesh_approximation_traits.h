#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_TRAITS_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_TRAITS_H

#include <list>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <boost/graph/graph_traits.hpp>

namespace CGAL
{
template<typename TriangleMesh,
  typename GeomTraits = typename TriangleMesh::Traits>
struct PlaneProxy
{
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Plane_3 Plane_3;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  face_descriptor seed;
  Vector_3 normal;
  Plane_3 fit_plane;
};

template<typename TriangleMesh,
  typename FacetNormalMap,
  typename FacetAreaMap,
  typename GeomTraits = typename TriangleMesh::Traits,
  typename PlaneProxy = CGAL::PlaneProxy<TriangleMesh, GeomTraits> >
struct L21Metric
{
  L21Metric(const FacetNormalMap &normal_pmap, const FacetAreaMap &area_pmap)
    : normal_pmap(normal_pmap),
    area_pmap(area_pmap) {
    GeomTraits traits;
    scalar_product_functor = traits.compute_scalar_product_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
  }

  typedef PlaneProxy Proxy;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename GeomTraits::Compute_scalar_product_3 Compute_scalar_product_3;
  typedef typename FacetAreaMap::key_type face_descriptor;

  FT operator()(const face_descriptor &f, const Proxy &px) const {
    Vector_3 v = sum_functor(normal_pmap[f], scale_functor(px.normal, FT(-1)));
    return area_pmap[f] * scalar_product_functor(v, v);
  }

  const FacetNormalMap normal_pmap;
  const FacetAreaMap area_pmap;
  Construct_scaled_vector_3 scale_functor;
  Compute_scalar_product_3 scalar_product_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

template<typename TriangleMesh,
  typename FacetNormalMap,
  typename FacetAreaMap,
  typename GeomTraits = typename TriangleMesh::Traits,
  typename PlaneProxy = CGAL::PlaneProxy<TriangleMesh, GeomTraits> >
struct L21ProxyFitting
{
  L21ProxyFitting(const FacetNormalMap &normal_pmap, const FacetAreaMap &area_pmap)
    : normal_pmap(normal_pmap), area_pmap(area_pmap) {
    GeomTraits traits;
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
  }

  typedef PlaneProxy Proxy;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;

  // Fit and construct a proxy
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

  const FacetNormalMap normal_pmap;
  const FacetAreaMap area_pmap;
  Construct_scaled_vector_3 scale_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

template<typename TriangleMesh,
  typename VertexPointMap
    = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  typename GeomTraits = typename TriangleMesh::Traits>
struct PlaneFitting
{
  PlaneFitting(const TriangleMesh &_mesh, const VertexPointMap &_point_pmap)
    : mesh(_mesh), point_pmap(_point_pmap) {
    GeomTraits traits;
    vector_functor = traits.construct_vector_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
  }

  PlaneFitting(const TriangleMesh &_mesh)
    : mesh(_mesh),
    point_pmap(get(boost::vertex_point, const_cast<TriangleMesh &>(_mesh))) {
    GeomTraits traits;
    vector_functor = traits.construct_vector_3_object();
    sum_functor = traits.construct_sum_of_vectors_3_object();
    scale_functor = traits.construct_scaled_vector_3_object();
  }

  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Plane_3 Plane_3;
  typedef typename GeomTraits::Construct_vector_3 Construct_vector_3;
  typedef typename GeomTraits::Construct_scaled_vector_3 Construct_scaled_vector_3;
  typedef typename GeomTraits::Construct_sum_of_vectors_3 Construct_sum_of_vectors_3;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

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

  const TriangleMesh &mesh;
  const VertexPointMap point_pmap;
  Construct_vector_3 vector_functor;
  Construct_scaled_vector_3 scale_functor;
  Construct_sum_of_vectors_3 sum_functor;
};

template<typename TriangleMesh,
  typename FacetAreaMap,
  typename VertexPointMap
    = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  typename GeomTraits = typename TriangleMesh::Traits,
  typename PlaneProxy = CGAL::PlaneProxy<TriangleMesh, GeomTraits> >
struct L2Metric
{
  L2Metric(const TriangleMesh &_mesh,
    const FacetAreaMap &_area_pmap,
    const VertexPointMap &_point_pmap)
    : mesh(_mesh), area_pmap(_area_pmap), point_pmap(_point_pmap) {}

  L2Metric(const TriangleMesh &_mesh,
    const FacetAreaMap &_area_pmap)
    : mesh(_mesh), area_pmap(_area_pmap),
    point_pmap(get(boost::vertex_point, const_cast<TriangleMesh &>(_mesh))) {}

  typedef PlaneProxy Proxy;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  FT operator()(const face_descriptor &f, const PlaneProxy &px) const {
    halfedge_descriptor he = halfedge(f, mesh);
    const Point_3 &p0 = point_pmap[source(he, mesh)];
    const Point_3 &p1 = point_pmap[target(he, mesh)];
    const Point_3 &p2 = point_pmap[target(next(he, mesh), mesh)];
    FT sq_d0 = CGAL::squared_distance(p0, px.fit_plane);
    FT sq_d1 = CGAL::squared_distance(p1, px.fit_plane);
    FT sq_d2 = CGAL::squared_distance(p2, px.fit_plane);
    FT d0(std::sqrt(CGAL::to_double(sq_d0)));
    FT d1(std::sqrt(CGAL::to_double(sq_d1)));
    FT d2(std::sqrt(CGAL::to_double(sq_d2)));

    return (sq_d0 + sq_d1 + sq_d2 + d0 * d1 + d1 * d2 + d2 * d0) * area_pmap[f] / FT(6);
  }

  const FacetAreaMap area_pmap;
  const VertexPointMap point_pmap;
  const TriangleMesh &mesh;
};

template<typename TriangleMesh,
  typename VertexPointMap
    = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  typename GeomTraits = typename TriangleMesh::Traits,
  typename PlaneProxy = CGAL::PlaneProxy<TriangleMesh, GeomTraits> >
struct L2ProxyFitting
{
  typedef PlaneProxy Proxy;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Triangle_3 Triangle_3;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  L2ProxyFitting(const TriangleMesh &_mesh, const VertexPointMap &_point_pmap)
    : mesh(_mesh), point_pmap(_point_pmap) {}

  L2ProxyFitting(const TriangleMesh &_mesh)
    : mesh(_mesh),
    point_pmap(get(boost::vertex_point, const_cast<TriangleMesh &>(_mesh))) {}

  template<typename FacetIterator>
  Proxy operator()(const FacetIterator beg, const FacetIterator end) const {
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
    Proxy px;
    CGAL::linear_least_squares_fitting_3(
      tris.begin(),
      tris.end(),
      px.fit_plane,
      CGAL::Dimension_tag<2>());

    return px;
  }

  const TriangleMesh &mesh;
  const VertexPointMap point_pmap;
};

template<typename TriangleMesh,
  typename VertexPointMap
    = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  typename GeomTraits = typename TriangleMesh::Traits>
struct PCAPlaneFitting
{
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Plane_3 Plane_3;
  typedef typename GeomTraits::Triangle_3 Triangle_3;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  PCAPlaneFitting(const TriangleMesh &_mesh,  const VertexPointMap &_point_pmap)
    : mesh(_mesh), point_pmap(_point_pmap) {}

  PCAPlaneFitting(const TriangleMesh &_mesh)
    : mesh(_mesh),
    point_pmap(get(boost::vertex_point, const_cast<TriangleMesh &>(_mesh))) {}

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

  const TriangleMesh &mesh;
  const VertexPointMap point_pmap;
};
} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_TRAITS_H
