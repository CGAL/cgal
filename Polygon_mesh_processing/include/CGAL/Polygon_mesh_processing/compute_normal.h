// Copyright (c) 2015-2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois
//                 Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H
#define CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H

#include <CGAL/license/Polygon_mesh_processing/Compute_normal.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Origin.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>
#include <limits>
#include <utility>
#include <vector>

#ifdef CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP
# ifndef CGAL_PMP_COMPUTE_NORMAL_DEBUG
#  define CGAL_PMP_COMPUTE_NORMAL_DEBUG
# endif
#endif

#ifdef CGAL_PMP_COMPUTE_NORMAL_DEBUG
#include <fstream>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#endif

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template <class GT>
void normalize(typename GT::Vector_3& v, const GT& traits)
{
  typedef typename GT::FT                                             FT;

  //If the vector is small enough, approx_sqrt might return 0, and then we get nan values.
  //To avoid that, we check the resulted norm. If it is 0, we don't normalize.
  const FT norm = CGAL::approximate_sqrt(traits.compute_squared_length_3_object()(v));
  if(norm != FT(0))
  {
    v = traits.construct_divided_vector_3_object()(v, norm);
  }
}

template<typename Point, typename GT>
typename GT::Vector_3
triangle_normal(const Point& p0, const Point& p1, const Point& p2, const GT& traits)
{
  typedef typename GT::FT                                             FT;

  typename GT::Vector_3 n = traits.construct_cross_product_vector_3_object()(
                              traits.construct_vector_3_object()(p1, p2),
                              traits.construct_vector_3_object()(p1, p0));

  //cross-product(AB, AC)'s norm is the area of the parallelogram
  //formed by these 2 vectors.
  //the triangle's area is half of it
  return traits.construct_scaled_vector_3_object()(n, FT(1)/FT(2));
}

template<typename Point, typename PM, typename VertexPointMap, typename Vector, typename GT>
void sum_normals(const PM& pmesh,
                 typename boost::graph_traits<PM>::face_descriptor f,
                 VertexPointMap vpmap,
                 Vector& sum,
                 const GT& traits)
{
  typedef typename boost::graph_traits<PM>::vertex_descriptor           vertex_descriptor;
  typedef typename boost::graph_traits<PM>::halfedge_descriptor         halfedge_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::reference    Point_ref;

  halfedge_descriptor he = halfedge(f, pmesh);
  vertex_descriptor v = source(he, pmesh);
  const Point_ref pv = get(vpmap, v);

  while(v != target(next(he, pmesh), pmesh))
  {
    const Point_ref pvn = get(vpmap, target(he, pmesh));
    const Point_ref pvnn = get(vpmap, target(next(he, pmesh), pmesh));

    const Vector n = internal::triangle_normal(pv, pvn, pvnn, traits);

#ifdef CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP
    std::cout << "Normal of " << f << " pts: " << pv << " ; " << pvn << " ; " << pvnn << std::endl;
    std::cout << " --> " << n << std::endl;
#endif

    sum = traits.construct_sum_of_vectors_3_object()(sum, n);

    he = next(he, pmesh);
  }
}

} // namespace internal

/**
* \ingroup PMP_normal_grp
* computes the outward unit vector normal to face `f`.
* @tparam PolygonMesh a model of `FaceGraph`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param f the face whose normal is computed
* @param pmesh the polygon mesh containing `f`
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* @return the computed normal. The return type is a 3D vector type. It is
* either deduced from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
* or from the geometric traits class deduced from the point property map
* of `pmesh`.
*
* \warning This function involves a square root computation.
* If the field type (`FT`) of the traits does not support the `sqrt()` operation,
* the square root computation will be performed approximately.
*/
template <typename PolygonMesh, typename NamedParameters>
#ifdef DOXYGEN_RUNNING
Vector_3
#else
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Vector_3
#endif
compute_face_normal(typename boost::graph_traits<PolygonMesh>::face_descriptor f,
                    const PolygonMesh& pmesh,
                    const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type               GT;
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type     VPMap;
  VPMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                 get_const_property_map(vertex_point, pmesh));

  typedef typename GT::Point_3                                                     Point;
  typedef typename GT::Vector_3                                                    Vector_3;

  Vector_3 normal = traits.construct_vector_3_object()(CGAL::NULL_VECTOR);
  internal::sum_normals<Point>(pmesh, f, vpmap, normal, traits);

  if(!traits.equal_3_object()(normal, CGAL::NULL_VECTOR))
    internal::normalize(normal, traits);

  return normal;
}

template <typename PolygonMesh>
typename GetGeomTraits<PolygonMesh>::type::Vector_3
compute_face_normal(typename boost::graph_traits<PolygonMesh>::face_descriptor f,
                    const PolygonMesh& pmesh)
{
  return compute_face_normal(f, pmesh, CGAL::parameters::all_default());
}

/**
* \ingroup PMP_normal_grp
* computes the outward unit vector normal for all faces of the polygon mesh.
* @tparam PolygonMesh a model of `FaceGraph`
* @tparam Face_normal_map a model of `WritablePropertyMap` with
    `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
    `Kernel::Vector_3` as value type.
*
* @param pmesh the polygon mesh
* @param face_normals the property map in which the normals are written
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \warning This function involves a square root computation.
* If the field type (`FT`) of the traits does not support the `sqrt()` operation,
* the square root computation will be performed approximately.
*/
template <typename PolygonMesh, typename Face_normal_map, typename NamedParameters>
void compute_face_normals(const PolygonMesh& pmesh,
                          Face_normal_map face_normals,
                          const NamedParameters& np)
{
  typedef typename GetGeomTraits<PolygonMesh,NamedParameters>::type Kernel;

  for(typename boost::graph_traits<PolygonMesh>::face_descriptor f : faces(pmesh))
  {
    typename Kernel::Vector_3 vec = compute_face_normal(f, pmesh, np);
    put(face_normals, f, vec);
#ifdef CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP
    std::cout << "normal at face " << f << " is " << get(face_normals, f) << std::endl;
#endif
  }
}

template <typename PolygonMesh, typename Face_normal_map>
void compute_face_normals(const PolygonMesh& pmesh, Face_normal_map face_normals)
{
  compute_face_normals(pmesh, face_normals, CGAL::parameters::all_default());
}

namespace internal {

enum Vertex_normal_type {
  NO_WEIGHT = 0,
  SIN_WEIGHT,
  MOST_VISIBLE
};

template <typename GT>
bool almost_equal(const typename GT::Vector_3& v1, const typename GT::Vector_3& v2,
                  const GT& traits)
{
  typedef typename GT::FT                                                     FT;

  // We are doing a lot of likely inexact constructions to compute min circles on the sphere and
  // these degenerate cases often happen (e.g. almost flat surfaces), so we don't want to to the usual
  // "switch to the exact kernel" in degenerate cases robustness trick. Rather, use a fixed tolerance
  // to assimilate vectors and scalar products.
  //
  // Tolerance is (arbitrarily) chosen as theta := 0.01°, thus we want sp(v1,v2) >= cos(theta)
  const FT cos_theta = 0.99999998476912910;
  return traits.compute_scalar_product_3_object()(v1, v2) >= cos_theta;
}

template <typename PolygonMesh, typename FaceNormalVector, typename K>
bool does_enclose_other_normals(const std::size_t i, const std::size_t j, const std::size_t k,
                                const typename K::Vector_3& nb,
                                const typename K::FT sp_bi,
                                const std::vector<typename boost::graph_traits<PolygonMesh>::face_descriptor>& incident_faces,
                                const FaceNormalVector& face_normals,
                                const K& traits)
{
  typedef typename K::FT                                                   FT;
  typedef typename boost::property_traits<FaceNormalVector>::reference     Vector_ref;

  typename K::Compute_scalar_product_3 sp = traits.compute_scalar_product_3_object();

  const FT nbn = CGAL::approximate_sqrt(traits.compute_squared_length_3_object()(nb));

  // check that this min circle defined by the diameter contains the other points
  const std::size_t nif = incident_faces.size();
  for(std::size_t l=0; l<nif; ++l)
  {
    if(l == i || l == j || l == k)
      continue;

    const Vector_ref nl = get(face_normals, incident_faces[l]);
    if(nl == CGAL::NULL_VECTOR)
      continue;

    // this is a bound on how much the scalar product between (v1,v2) and (v1, v3) can change
    // when the angle changes theta_bound := 0.01°
    // The weird number is thus := max_(theta_i, theta_j)[abs(std::cos(theta_i), std::cos(theta_j))]
    // with theta_j - theta_i = theta_bound
    const FT sp_diff_bound = nbn * 0.00017453292431333;
    const FT sp_bl = sp(nb, nl);

    // norm of nl is 1 by construction
    if(CGAL::abs(sp_bi - sp_bl) <= sp_diff_bound)
      continue;

    if(sp_bl < sp_bi)
      return false;
  }

  return true;
}

template <typename GT>
typename GT::Vector_3 compute_normals_bisector(const typename GT::Vector_3& ni,
                                               const typename GT::Vector_3& nj,
                                               const GT& traits)
{
  if(traits.equal_3_object()(ni, nj))
    return ni;

  return traits.construct_sum_of_vectors_3_object()(ni, nj); // not normalized
}

template <typename GT>
typename GT::Vector_3 compute_normals_bisector(const typename GT::Vector_3& ni,
                                               const typename GT::Vector_3& nj,
                                               const typename GT::Vector_3& nk,
                                               const GT& traits)
{
  typedef typename GT::FT                                                  FT;
  typedef typename GT::Point_3                                             Point_3;
  typedef typename GT::Vector_3                                            Vector_3;

  typename GT::Construct_scaled_vector_3 cslv_3 = traits.construct_scaled_vector_3_object();
  typename GT::Construct_sum_of_vectors_3 csv_3 = traits.construct_sum_of_vectors_3_object();
  typename GT::Construct_vector_3 cv_3 = traits.construct_vector_3_object();

  Vector_3 nb = cv_3(CGAL::NULL_VECTOR);

  if(almost_equal(ni, nj, traits) || nk == CGAL::NULL_VECTOR)
  {
    if(almost_equal(nj, nk, traits))
      nb = ni;
    else // ni == nj, but nij != nk
      nb = compute_normals_bisector(nj, nk, traits);
  }
  else if(almost_equal(ni, nk, traits) || nj == CGAL::NULL_VECTOR) // ni != nj
  {
    nb = compute_normals_bisector(nj, nk, traits);
  }
  else if(almost_equal(nj, nk, traits) || ni == CGAL::NULL_VECTOR) // ni != nj, ni != nk
  {
    nb = compute_normals_bisector(ni, nk, traits);
  }
  else
  {
    CGAL_assertion(ni != nj);
    CGAL_assertion(ni != nk);
    CGAL_assertion(nj != nk);
    CGAL_assertion(!traits.collinear_3_object()(CGAL::ORIGIN + ni, CGAL::ORIGIN + nj, CGAL::ORIGIN + nk));

#ifdef CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP
    std::cout << "Triplet: ni[" << ni << "] nj[" << nj << "] nk[" << nk << "]" << std::endl;
#endif

    const Point_3 c = traits.construct_circumcenter_3_object()(CGAL::ORIGIN + ni, CGAL::ORIGIN + nj, CGAL::ORIGIN + nk);
    if(c == CGAL::ORIGIN)
    {
      // will happen if the three vectors live in the same plan, return some weighted sum
      const FT third = FT(1)/FT(3);
      return csv_3(csv_3(cslv_3(ni, third), cslv_3(nj, third)), cslv_3(nk, third));
    }

    nb = cv_3(CGAL::ORIGIN, c); // not normalized
  }

  return nb;
}

template <typename PolygonMesh, typename FaceNormalVector, typename GT>
typename GT::Vector_3
compute_most_visible_normal_2_points(std::vector<typename boost::graph_traits<PolygonMesh>::face_descriptor>& incident_faces,
                                     const FaceNormalVector& face_normals,
                                     const GT& traits)
{
  typedef typename GT::FT                                                  FT;
  typedef typename GT::Vector_3                                            Vector_3;
  typedef typename boost::property_traits<FaceNormalVector>::reference     Vector_ref;

  typename GT::Compute_scalar_product_3 sp_3 = traits.compute_scalar_product_3_object();
  typename GT::Construct_vector_3 cv_3 = traits.construct_vector_3_object();

#ifdef CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP
  std::cout << "Trying to find enclosing normal with 2 normals" << std::endl;
#endif

  FT min_sp = -1;
  Vector_3 n = cv_3(CGAL::NULL_VECTOR);

  const std::size_t nif = incident_faces.size();
  for(std::size_t i=0; i<nif; ++i)
  {
    for(std::size_t j=i+1; j<nif; ++j)
    {
      const Vector_ref ni = get(face_normals, incident_faces[i]);
      const Vector_ref nj = get(face_normals, incident_faces[j]);

      const Vector_3 nb = compute_normals_bisector(ni, nj, traits);

      // Degeneracies like ni == -nj or a numerical error in the construction of 'nb' can happen.
      if(traits.equal_3_object()(nb, CGAL::NULL_VECTOR))
        return CGAL::NULL_VECTOR;

      FT sp_bi = sp_3(nb, ni);
      sp_bi = (std::max)(FT(0), sp_bi);
      if(sp_bi <= min_sp)
        continue;

      if(!does_enclose_other_normals<PolygonMesh>(i, j, -1 /*NA*/, nb, sp_bi, incident_faces, face_normals, traits))
        continue;

      min_sp = sp_bi;
      n = nb;
    }
  }

  return n;
}

template <typename PolygonMesh, typename FaceNormalVector, typename GT>
typename GT::Vector_3
compute_most_visible_normal_3_points(const std::vector<typename boost::graph_traits<PolygonMesh>::face_descriptor>& incident_faces,
                                     const FaceNormalVector& face_normals,
                                     const GT& traits)
{
  typedef typename GT::FT                                                  FT;
  typedef typename GT::Vector_3                                            Vector_3;
  typedef typename boost::property_traits<FaceNormalVector>::reference     Vector_ref;

#ifdef CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP
  std::cout << "Trying to find enclosing normal with 3 normals" << std::endl;
#endif

  FT min_sp = -1;

  Vector_3 n = traits.construct_vector_3_object()(CGAL::NULL_VECTOR);

  const std::size_t nif = incident_faces.size();
  for(std::size_t i=0; i<nif; ++i)
  {
    for(std::size_t j=i+1; j<nif; ++j)
    {
      for(std::size_t k=j+1; k<nif; ++k)
      {
        const Vector_ref ni = get(face_normals, incident_faces[i]);
        const Vector_ref nj = get(face_normals, incident_faces[j]);
        const Vector_ref nk = get(face_normals, incident_faces[k]);

        if(ni == CGAL::NULL_VECTOR || nj == CGAL::NULL_VECTOR || nk == CGAL::NULL_VECTOR)
          continue;

        Vector_3 nb = compute_normals_bisector(ni, nj, nk, traits);
        if(traits.equal_3_object()(nb, CGAL::NULL_VECTOR))
          return nb;

        FT sp_bi = traits.compute_scalar_product_3_object()(nb, ni);
        if(sp_bi < FT(0))
        {
          nb = traits.construct_opposite_vector_3_object()(nb);
          sp_bi = - sp_bi;
        }

        if(sp_bi <= min_sp)
          continue;

        if(!does_enclose_other_normals<PolygonMesh>(i, j, k, nb, sp_bi, incident_faces, face_normals, traits))
          continue;

        min_sp = sp_bi;
        n = nb;
      }
    }
  }

#ifdef CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP
  std::cout << "Best normal from 3-normals-approach: " << n << std::endl;
#endif

  return n;
}

// Inspired by Aubry et al. On the most 'normal' normal
template <typename PolygonMesh, typename FaceNormalVector, typename GT>
typename GT::Vector_3
compute_vertex_normal_most_visible_min_circle(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
                                              const FaceNormalVector& face_normals,
                                              const PolygonMesh& pmesh,
                                              const GT& traits)
{
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor       face_descriptor;

  typedef typename GT::Vector_3                                            Vector_3;

  std::vector<face_descriptor> incident_faces;
  for(face_descriptor f : CGAL::faces_around_target(halfedge(v, pmesh), pmesh))
  {
    if(f == boost::graph_traits<PolygonMesh>::null_face())
      continue;

    incident_faces.push_back(f);
  }

  if(incident_faces.size() == 1)
    return get(face_normals, incident_faces.front());

  Vector_3 res = compute_most_visible_normal_2_points<PolygonMesh>(incident_faces, face_normals, traits);

  if(res != CGAL::NULL_VECTOR) // found a valid normal through 2 point min circle
    return res;

  // The vertex has only two incident faces with opposite normals (fold)...
  // @todo devise something based on the directions of the 2/3/4 incident edges?
  if(incident_faces.size() == 2 && res == CGAL::NULL_VECTOR)
    return res;

  CGAL_assertion(incident_faces.size() >= 2);

  return compute_most_visible_normal_3_points<PolygonMesh>(incident_faces, face_normals, traits);
}

template <typename PolygonMesh, typename FaceNormalVector, typename VertexPointMap, typename GT>
typename GT::Vector_3
compute_vertex_normal_as_sum_of_weighted_normals(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
                                                 const Vertex_normal_type& vn_type,
                                                 const FaceNormalVector& face_normals,
                                                 const VertexPointMap& vpmap,
                                                 const PolygonMesh& pmesh,
                                                 const GT& traits)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  typedef typename GT::FT                                                     FT;
  typedef typename GT::Vector_3                                               Vector_3;
  typedef typename boost::property_traits<FaceNormalVector>::reference        Vector_ref;

  typename GT::Construct_vector_3 cv_3 = traits.construct_vector_3_object();
  typename GT::Compute_squared_length_3 csl_3 = traits.compute_squared_length_3_object();

#ifdef CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP
  std::cout << "Compute normal as weighted sums; type: " << vn_type << std::endl;
#endif

  Vector_3 normal = cv_3(CGAL::NULL_VECTOR);

  halfedge_descriptor h = halfedge(v, pmesh);
  if(h == boost::graph_traits<PolygonMesh>::null_halfedge())
    return normal;

  halfedge_descriptor end = h;
  do
  {
    if(!is_border(h, pmesh))
    {
      if(vn_type == NO_WEIGHT)
      {
        const Vector_ref n = get(face_normals, face(h, pmesh));
        normal = traits.construct_sum_of_vectors_3_object()(normal, n);
      }
      else if(vn_type == SIN_WEIGHT)
      {
        const Vector_3 v1 = cv_3(get(vpmap, v), get(vpmap, source(h, pmesh)));
        const Vector_3 v2 = cv_3(get(vpmap, v), get(vpmap, target(next(h, pmesh), pmesh)));

        //v(i) and v(i+1) must be seen in ccw order, from v, so we reverse v1 and v2
        Vector_3 n = traits.construct_cross_product_vector_3_object()(v2, v1);
        const FT den = CGAL::approximate_sqrt(csl_3(v1) * csl_3(v2));

        if(den == FT(0))
        {
#ifdef CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP
          std::cout << "Null denominator, switching to no weights" << std::endl;
#endif

          return compute_vertex_normal_as_sum_of_weighted_normals(v, NO_WEIGHT, face_normals, vpmap, pmesh, traits);
        }

        n = traits.construct_scaled_vector_3_object()(n, FT(1) / den);
        normal = traits.construct_sum_of_vectors_3_object()(normal, n);
      }
      else
      {
        std::cerr << "Error: unknown vertex normal type" << std::endl;
        CGAL_assertion(false);
        return CGAL::NULL_VECTOR;
      }
    }

    h = opposite(next(h, pmesh), pmesh);
  }
  while(h != end);

  return normal;
}

} // end namespace internal

/**
* \ingroup PMP_normal_grp
* computes the unit normal at vertex `v` as the average of the normals of incident faces.
* @tparam PolygonMesh a model of `FaceGraph`
*
* @param v the vertex whose normal is computed
* @param pmesh the polygon mesh containing `v`
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* @return the computed normal. The return type is a 3D vector type. It is
* either deduced from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
* or the geometric traits class deduced from the point property map
* of `pmesh`.
*
* \warning This function involves a square root computation.
* If the field type (`FT`) of the traits does not support the `sqrt()` operation,
* the square root computation will be performed approximately.
*/
template<typename PolygonMesh, typename NamedParameters>
#ifdef DOXYGEN_RUNNING
Vector_3
#else
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Vector_3
#endif
compute_vertex_normal(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
                      const PolygonMesh& pmesh,
                      const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::is_default_parameter;
  using parameters::get_parameter;

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor          face_descriptor;

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type          GT;
  typedef typename GT::Vector_3                                               Vector_3;
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VPMap;
  VPMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                 get_const_property_map(vertex_point, pmesh));

  typedef std::map<face_descriptor, Vector_3>                                 Face_vector_map;
  typedef boost::associative_property_map<Face_vector_map>                    Default_map;

  typedef typename internal_np::Lookup_named_param_def<internal_np::face_normal_t,
                                                       NamedParameters,
                                                       Default_map>::type     Face_normal_map;
  Face_vector_map default_fvmap;
  Face_normal_map face_normals = choose_parameter(get_parameter(np, internal_np::face_normal),
                                                  Default_map(default_fvmap));
  const bool must_compute_face_normals = is_default_parameter(get_parameter(np, internal_np::face_normal));

#ifdef CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP
  std::cout << "<----- compute vertex normal at " << get(vpmap, v)
            << ", must compute face normals? " << must_compute_face_normals << std::endl;
#endif

  // handle isolated vertices
  halfedge_descriptor he = halfedge(v, pmesh);
  if(he == boost::graph_traits<PolygonMesh>::null_halfedge())
    return CGAL::NULL_VECTOR;

  if(must_compute_face_normals)
  {
    for(face_descriptor f : CGAL::faces_around_target(halfedge(v, pmesh), pmesh))
    {
      if(f == boost::graph_traits<PolygonMesh>::null_face())
        continue;

      put(face_normals, f, compute_face_normal(f, pmesh, np));
    }
  }

#ifdef CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP
  std::cout << "Incident face normals:" << std::endl;
  for(halfedge_descriptor h : CGAL::halfedges_around_target(v, pmesh))
  {
    if(!is_border(h, pmesh))
      std::cout << "get normal at f " << face(h, pmesh) << " : " << get(face_normals, face(h, pmesh)) << std::endl;
  }
#endif

  Vector_3 normal = internal::compute_vertex_normal_most_visible_min_circle(v, face_normals, pmesh, traits);
  if(traits.equal_3_object()(normal, CGAL::NULL_VECTOR)) // can't always find a most visible normal
  {
#ifdef CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP
    std::cout << "Failed to find most visible normal, use weighted sum of normals" << std::endl;
#endif
    normal = internal::compute_vertex_normal_as_sum_of_weighted_normals(
               v, internal::SIN_WEIGHT, face_normals, vpmap, pmesh, traits);
  }

  if(!traits.equal_3_object()(normal, CGAL::NULL_VECTOR))
    internal::normalize(normal, traits);

  return normal;
}

template <typename PolygonMesh>
typename GetGeomTraits<PolygonMesh>::type::Vector_3
compute_vertex_normal(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
                      const PolygonMesh& pmesh)
{
  return compute_vertex_normal(v, pmesh, CGAL::parameters::all_default());
}

/**
* \ingroup PMP_normal_grp
* computes the outward unit vector normal for all vertices of the polygon mesh.
*
* @tparam PolygonMesh a model of `FaceListGraph`
* @tparam VertexNormalMap a model of `WritablePropertyMap` with
*                         `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and
*                         the return type of `compute_vertex_normal()` as value type.
*
* @param pmesh the polygon mesh
* @param vertex_normals the property map in which the normals are written
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \warning This function involves a square root computation.
* If the field type (`FT`) of the traits does not support the `sqrt()` operation,
* the square root computation will be performed approximately.
*/
template <typename PolygonMesh, typename VertexNormalMap, typename NamedParameters>
void compute_vertex_normals(const PolygonMesh& pmesh,
                            VertexNormalMap vertex_normals,
                            const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::is_default_parameter;
  using parameters::get_parameter;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor           vertex_descriptor;

  typedef typename GetGeomTraits<PolygonMesh,NamedParameters>::type              GT;
  typedef typename GT::Vector_3                                                  Vector_3;

  typedef CGAL::dynamic_face_property_t<Vector_3>                                Face_normal_tag;
  typedef typename boost::property_map<PolygonMesh, Face_normal_tag>::const_type Face_normal_dmap;

#ifdef CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type   VPMap;
  VPMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                 get_const_property_map(vertex_point, pmesh));
#endif

  typedef typename internal_np::Lookup_named_param_def<internal_np::face_normal_t,
                                                       NamedParameters,
                                                       Face_normal_dmap>::type   Face_normal_map;
  Face_normal_map face_normals = choose_parameter(get_parameter(np, internal_np::face_normal),
                                                  get(Face_normal_tag(), pmesh));
  const bool must_compute_face_normals = is_default_parameter(get_parameter(np, internal_np::face_normal));

  if(must_compute_face_normals)
    compute_face_normals(pmesh, face_normals, np);

#ifdef CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP
  std::ofstream out("computed_normals.cgal.polylines.txt");
  const Bbox_3 bb = bbox(pmesh, np);
  const typename GT::FT bbox_diagonal = CGAL::sqrt(CGAL::square(bb.xmax() - bb.xmin()) +
                                                   CGAL::square(bb.ymax() - bb.ymin()) +
                                                   CGAL::square(bb.zmax() - bb.zmin()));
#endif

  for(vertex_descriptor v : vertices(pmesh))
  {
    const Vector_3 n = compute_vertex_normal(v, pmesh, np.face_normal_map(face_normals));
    put(vertex_normals, v, n);

#ifdef CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP
    out << "2 " << get(vpmap, v) << " "
                << get(vpmap, v) + traits.construct_scaled_vector_3_object()(n, 0.1 * bbox_diagonal) << "\n";
#endif
  }
}

template <typename PolygonMesh, typename VertexNormalMap>
void compute_vertex_normals(const PolygonMesh& pmesh, VertexNormalMap vertex_normals)
{
  compute_vertex_normals(pmesh, vertex_normals, CGAL::parameters::all_default());
}

/**
* \ingroup PMP_normal_grp
* computes the outward unit vector normal for all vertices and faces of the polygon mesh.
*
* @tparam PolygonMesh a model of `FaceListGraph`
* @tparam VertexNormalMap a model of `WritablePropertyMap` with
*    `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and
*    `Kernel::Vector_3` as value type.
* @tparam FaceNormalMap a model of `ReadWritePropertyMap` with
*    `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
*    `Kernel::Vector_3` as value type.
*
* @param pmesh the polygon mesh
* @param vertex_normals the property map in which the vertex normals are written
* @param face_normals the property map in which the face normals are written
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \warning This function involves a square root computation.
* If the field type (`FT`) of the traits does not support the `sqrt()` operation,
* the square root computation will be performed approximately.
*/
template <typename PolygonMesh,
          typename VertexNormalMap, typename FaceNormalMap,
          typename NamedParameters>
void compute_normals(const PolygonMesh& pmesh,
                     VertexNormalMap vertex_normals,
                     FaceNormalMap face_normals,
                     const NamedParameters& np)
{
  compute_face_normals(pmesh, face_normals, np);
  compute_vertex_normals(pmesh, vertex_normals, np.face_normal_map(face_normals));
}

template <typename PolygonMesh, typename VertexNormalMap, typename FaceNormalMap>
void compute_normals(const PolygonMesh& pmesh,
                     VertexNormalMap vertex_normals,
                     FaceNormalMap face_normals)
{
  compute_normals(pmesh, vertex_normals, face_normals, CGAL::parameters::all_default());
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H
