// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Jane Tournois


#ifndef CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H
#define CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H

#include <CGAL/license/Polygon_mesh_processing/Compute_normal.h>


#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/Origin.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Kernel_traits.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <boost/type_traits.hpp>

namespace CGAL{

namespace Polygon_mesh_processing{

namespace internal {

  template <class GT>
  void normalize(typename GT::Vector_3& v, const GT& traits)
  {
    typename GT::FT norm = CGAL::approximate_sqrt(
        traits.compute_squared_length_3_object()(v));
    //If the vector is small enough, approx_sqrt might return 0, and then we get nan values.
    //To avoid that, we check the resulted norm. If it is 0, we don't normalize.
    if(norm != 0)
    {
      v = traits.construct_divided_vector_3_object()(v, norm );
    }
  }

  template<typename Point
         , typename GT>
  typename GT::Vector_3
  triangle_normal(const Point& p0, const Point& p1, const Point& p2
                , const GT& traits)
  {
    typename GT::Vector_3 n = traits.construct_cross_product_vector_3_object()(
      traits.construct_vector_3_object()(p1, p2),
      traits.construct_vector_3_object()(p1, p0));

    //cross-product(AB, AC)'s norm is the area of the parallelogram
    //formed by these 2 vectors.
    //the triangle's area is half of it
    return traits.construct_scaled_vector_3_object()(n, 0.5);
  }
}

template<typename Point, typename PM, typename VertexPointMap, typename Vector
       , typename GT>
void sum_normals(const PM& pmesh,
                 typename boost::graph_traits<PM>::face_descriptor f,
                 VertexPointMap vpmap,
                 Vector& sum,
                 const GT& traits)
{
  typedef typename boost::graph_traits<PM>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor he = halfedge(f, pmesh);
  vertex_descriptor v = source(he, pmesh);
  const Point& pv = get(vpmap, v);
  while (v != target(next(he, pmesh), pmesh))
  {
    const Point& pvn  = get(vpmap, target(he, pmesh));
    const Point& pvnn = get(vpmap, target(next(he, pmesh), pmesh));

    Vector n = internal::triangle_normal(pv, pvn, pvnn, traits);
    sum = traits.construct_sum_of_vectors_3_object()(sum, n);

    he = next(he, pmesh);
  }
}


/**
* \ingroup PMP_normal_grp
* computes the outward unit vector normal to face `f`.
* @tparam PolygonMesh a model of `FaceGraph`
* @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
*
* @param f the face on which the normal is computed
* @param pmesh the polygon mesh containing `f`
* @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
*   If this parameter is omitted, an internal property map for
*   `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
*    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
* \cgalNamedParamsEnd
*
* @return the computed normal. The return type is a 3D vector type. It is
* either deduced from the `geom_traits` \ref pmp_namedparameters "Named Parameters" if provided,
* or from the geometric traits class deduced from the point property map
* of `pmesh`.
*
* \warning This function involves a square root computation.
* If `Kernel::FT` does not have a `sqrt()` operation, the square root computation
* will be done approximately.
*/
template <typename PolygonMesh, typename NamedParameters>
#ifdef DOXYGEN_RUNNING
Vector_3
#else
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Vector_3
#endif
compute_face_normal(typename boost::graph_traits<PolygonMesh>::face_descriptor f
                    , const PolygonMesh& pmesh
                    , const NamedParameters& np)
{
  using boost::choose_param;
  using boost::get_param;

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
  GT traits = choose_param(get_param(np, internal_np::geom_traits), GT());

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VPMap;
  VPMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, pmesh));

  typedef typename GT::Point_3 Point;
  typedef typename GT::Vector_3 Vector;

  Vector normal = traits.construct_vector_3_object()(CGAL::NULL_VECTOR);
  sum_normals<Point>(pmesh, f, vpmap, normal, traits);

  if (!traits.equal_3_object()(normal, CGAL::NULL_VECTOR))
    internal::normalize(normal, traits);

  return normal;
}

///\cond SKIP_IN_MANUAL

// compute_face_normal overload
template <typename PolygonMesh>
typename GetGeomTraits<PolygonMesh>::type::Vector_3
compute_face_normal(typename boost::graph_traits<PolygonMesh>::face_descriptor f,
                    const PolygonMesh& pmesh)
{
  return compute_face_normal(f, pmesh,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

/// \endcond

/**
* \ingroup PMP_normal_grp
* computes the outward unit vector normal for all faces of the polygon mesh.
* @tparam PolygonMesh a model of `FaceGraph`
* @tparam FaceNormalMap a model of `WritablePropertyMap` with
    `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
    `Kernel::Vector_3` as value type.
*
* @param pmesh the polygon mesh
* @param fnm the property map in which the normals are written
* @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
*   If this parameter is omitted, an internal property map for
*   `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
*    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
* \cgalNamedParamsEnd
*
* \warning This function involves a square root computation.
* If `Kernel::FT` does not have a `sqrt()` operation, the square root computation
* will be done approximately.
*/
template <typename PolygonMesh
          , typename FaceNormalMap
          , typename NamedParameters>
void
compute_face_normals(const PolygonMesh& pmesh
                   , FaceNormalMap fnm
                   , const NamedParameters& np)
{
  typedef typename GetGeomTraits<PolygonMesh,NamedParameters>::type Kernel;

  for(typename boost::graph_traits<PolygonMesh>::face_descriptor f : faces(pmesh)){
    typename Kernel::Vector_3 vec = compute_face_normal(f, pmesh, np);
    put(fnm, f, vec);
  }
}

///\cond SKIP_IN_MANUAL

// compute_face_normals overload
template <typename PolygonMesh, typename FaceNormalMap>
void compute_face_normals(const PolygonMesh& pmesh, FaceNormalMap fnm)
{
  compute_face_normals(pmesh, fnm, CGAL::Polygon_mesh_processing::parameters::all_default());
}

/// \endcond


// @TMP
// -----------------------------------------------------------------
template<typename K>
typename K::FT fix_sine(typename K::FT sine)
{
  if(sine >= 1)
    return 1;
  else if(sine <= -1)
    return -1;
  else
    return sine;
}

template<typename K>
typename K::FT compute_angle_rad(const typename K::Vector_3& u,
                                 const typename K::Vector_3& v)
{
  typedef typename K::FT                           NT;

  // check
  NT product = CGAL::approximate_sqrt(u * u) * CGAL::approximate_sqrt(v * v);
  if(product == 0)
    return 0;

  // cosine
  NT dot = (u * v);
  NT cosine = dot / product;

  return std::acos(fix_sine<K>(cosine));
}

//                                                      ->  ->
// Returns the angle (in radians) of (P,Q,R) corner (i.e. QP, QR angle).
template<typename K>
typename K::FT compute_angle_rad(const typename K::Point_3& P,
                                 const typename K::Point_3& Q,
                                 const typename K::Point_3& R)
{
  typedef typename K::Vector_3                     Vector_3;

  Vector_3 u = P - Q;
  Vector_3 v = R - Q;

  return compute_angle_rad<K>(u, v);
}
// -----------------------------------------------------------------

template <typename GT>
bool almost_equal(const typename GT::Vector_3 v1, const typename GT::Vector_3 v2)
{
  return (CGAL::abs(1 - GT().compute_scalar_product_3_object()(v1, v2)) < 1e-10);
}

template <typename PolygonMesh, typename FaceNormalVector, typename K>
bool
does_enclose_other_normals(const int i, const int j, const int k,
                           const typename K::Vector_3& nb,
                           const typename K::FT sp_bi,
                           std::vector<typename boost::graph_traits<PolygonMesh>::face_descriptor> incident_faces,
                           const FaceNormalVector& normalized_face_normals,
                           const K& traits)
{
  typedef typename K::FT                                                  FT;
  typedef typename K::Vector_3                                            Vector_3;

  typename K::Compute_scalar_product_3 sp = traits.compute_scalar_product_3_object();

  // check that this min circle defined by the diameter contains the other points
  const std::size_t nif = incident_faces.size();
  for(int l=0; l<nif; ++l)
  {
    if(l == i || l == j|| l == k)
      continue;

    const Vector_3& nl = normalized_face_normals[incident_faces[l]];

    const FT sp_bl = sp(nb, nl);
    if(CGAL::abs(sp_bi - sp_bl) < 1e-10) // @tolerance
      continue;

    if(sp_bl < sp_bi)
      return false;
  }

  return true;
}

template <typename GT>
typename GT::Vector_3 compute_normals_bisector(const typename GT::Vector_3 ni,
                                               const typename GT::Vector_3 nj,
                                               const GT& traits)
{
  if(ni == nj)
    return ni;

  typename GT::Vector_3 nb = traits.construct_sum_of_vectors_3_object()(ni, nj);
  CGAL::Polygon_mesh_processing::internal::normalize(nb, traits); // prob not necessary @fixme

  return nb;
}

template <typename PolygonMesh, typename FaceNormalVector, typename GT>
std::pair<typename GT::Vector_3, bool>
compute_most_visible_normal_2_points(std::vector<typename boost::graph_traits<PolygonMesh>::face_descriptor> incident_faces,
                                     const FaceNormalVector& normalized_face_normals,
                                     const GT& traits)
{
  typedef typename GT::FT                                                  FT;
  typedef typename GT::Vector_3                                            Vector_3;

  typename GT::Compute_scalar_product_3 sp = traits.compute_scalar_product_3_object();

  FT min_sp = -1;

  Vector_3 n = CGAL::NULL_VECTOR;

  const std::size_t nif = incident_faces.size();
  for(int i=0; i<nif; ++i)
  {
    for(int j=i+1; j<nif; ++j)
    {
      const Vector_3& ni = normalized_face_normals[incident_faces[i]];
      const Vector_3& nj = normalized_face_normals[incident_faces[j]];

      Vector_3 nb = compute_normals_bisector(ni, nj, traits);
      const FT sp_bi = sp(nb, ni);

      CGAL_assertion(sp_bi >= 0);
      if(sp_bi <= min_sp)
        continue;

      if(!does_enclose_other_normals<PolygonMesh>(i, j, -1 /*NA*/, nb, sp_bi, incident_faces, normalized_face_normals, traits))
        continue;

      min_sp = sp_bi;
      n = nb;
    }
  }

  return std::make_pair(n, (n != CGAL::NULL_VECTOR));
}

template <typename PolygonMesh, typename FaceNormalVector, typename GT>
typename GT::Vector_3
compute_most_visible_normal_3_points(std::vector<typename boost::graph_traits<PolygonMesh>::face_descriptor> incident_faces,
                                     const FaceNormalVector& normalized_face_normals,
                                     const GT& traits)
{
  typedef typename GT::FT                                                  FT;
  typedef typename GT::Point_3                                             Point_3;
  typedef typename GT::Vector_3                                            Vector_3;

  typename GT::Compute_scalar_product_3 sp = traits.compute_scalar_product_3_object();
  typename GT::Construct_sum_of_vectors_3 csv = traits.construct_sum_of_vectors_3_object();

  FT min_sp = -1;

  Vector_3 n = CGAL::NULL_VECTOR;

  const std::size_t nif = incident_faces.size();

  for(int i=0; i<nif; ++i)
  {
    for(int j=i+1; j<nif; ++j)
    {
      for(int k=j+1; k<nif; ++k)
      {
        const Vector_3& ni = normalized_face_normals[incident_faces[i]];
        const Vector_3& nj = normalized_face_normals[incident_faces[j]];
        const Vector_3& nk = normalized_face_normals[incident_faces[k]];

        Vector_3 nb = CGAL::NULL_VECTOR;

        if(almost_equal<GT>(ni, nj))
        {
          if(almost_equal<GT>(nj, nk))
            nb = ni;
          else // ni == nj, but nj != nk
            nb = compute_normals_bisector(nj, nk, traits);
        }
        else if(almost_equal<GT>(ni, nk)) // ni != nj
        {
          nb = compute_normals_bisector(nj, nk, traits);
        }
        else if(almost_equal<GT>(nj, nk)) // ni != nj
        {
          nb = compute_normals_bisector(ni, nk, traits);
        }
        else
        {
          CGAL_assertion(ni != nj);
          CGAL_assertion(ni != nk);
          CGAL_assertion(nj != nk);
          CGAL_assertion(!CGAL::collinear(CGAL::ORIGIN + ni, CGAL::ORIGIN + nj, CGAL::ORIGIN + nk));

          const Point_3 c = CGAL::circumcenter(CGAL::ORIGIN + ni, CGAL::ORIGIN + nj, CGAL::ORIGIN + nk);
          CGAL_assertion(c != CGAL::ORIGIN);

          nb = traits.construct_vector_3_object()(CGAL::ORIGIN, c);
          CGAL::Polygon_mesh_processing::internal::normalize(nb, traits); // prob not necessary @fixme
        }

        CGAL_assertion(nb != CGAL::NULL_VECTOR);
        FT sp_bi = sp(nb, ni);

        if(sp_bi < 0)
        {
          nb = traits.construct_opposite_vector_3_object()(nb);
          sp_bi = - sp_bi;
        }

        if(sp_bi <= min_sp)
          continue;

        if(!does_enclose_other_normals<PolygonMesh>(i, j, k, nb, sp_bi, incident_faces, normalized_face_normals, traits))
          continue;

        min_sp = sp_bi;
        n = nb;
      }
    }
  }

  CGAL_assertion(n != CGAL::NULL_VECTOR);
  return n;
}

// Complexity is high, but valence is usually low, so that's ok
template <typename FaceNormalVector, typename PolygonMesh, typename GT>
typename GT::Vector_3
compute_vertex_normal_most_visible_min_circle(typename boost::graph_traits<PolygonMesh>::vertex_descriptor vd,
                                              const FaceNormalVector& normalized_face_normals,
                                              const PolygonMesh& pmesh,
                                              const GT& traits)
{
  typedef typename GT::FT                                                  FT;
  typedef typename GT::Vector_3                                            Vector_3;

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor       face_descriptor;

#ifdef CGAL_MOST_VISIBLE_NORMAL_VERBOSE
  std::cout << std::endl << std::endl;
  std::cout << "----------------------------------------------------------------------" << std::endl;
  std::cout << "compute_vertex_normal_most_visible_min_circle at " << pmesh.point(vd) << std::endl;
#endif

  halfedge_descriptor hd = halfedge(vd, pmesh);

  std::vector<face_descriptor> incident_faces;
  for(face_descriptor fd : CGAL::faces_around_target(hd, pmesh))
  {
    if(fd == boost::graph_traits<PolygonMesh>::null_face())
      continue;

    incident_faces.push_back(fd);
  }

  if(incident_faces.size() == 1)
    return normalized_face_normals[incident_faces.front()];

  std::pair<Vector_3, bool> res = compute_most_visible_normal_2_points<PolygonMesh>(incident_faces, normalized_face_normals, traits);
  if(res.second)
    return res.first;

  CGAL_assertion(incident_faces.size() > 2);

  return compute_most_visible_normal_3_points<PolygonMesh>(incident_faces, normalized_face_normals, traits);
}

template <typename FaceNormalVector, typename PolygonMesh, typename GT>
typename GT::Vector_3
compute_vertex_normal_most_visible_with_optimization(typename boost::graph_traits<PolygonMesh>::vertex_descriptor vd,
                                                     const FaceNormalVector& normalized_face_normals,
                                                     const PolygonMesh& pmesh,
                                                     const GT& traits)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor         halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor             face_descriptor;

  typedef typename GT::FT                                                        FT;
  typedef typename GT::Vector_3                                                  Vector_3;

  typename GT::Construct_scaled_vector_3 csclv = traits.construct_scaled_vector_3_object();
  typename GT::Construct_sum_of_vectors_3 csv = traits.construct_sum_of_vectors_3_object();
  typename GT::Compute_scalar_product_3 sp = traits.compute_scalar_product_3_object();

#ifdef CGAL_MOST_VISIBLE_NORMAL_VERBOSE
  std::cout << "----------------------------------------------------------------------" << std::endl;
  std::cout << "compute_vertex_normal_most_visible_with_optimization at " << pmesh.point(vd) << std::endl;
#endif

  halfedge_descriptor hd = halfedge(vd, pmesh);

  std::map<face_descriptor, FT> weights;
  Vector_3 normal = CGAL::NULL_VECTOR;

// initial weight based on the number of incident faces
#ifdef CGAL_PMP_MOST_VISIBLE_NORMAL_INITIALIZE_WITH_EQUAL_WEIGHTS
  typedef typename CGAL::Halfedge_around_target_iterator<PolygonMesh>::difference_type   difference_type;
  difference_type incident_faces_n = 0;

  for(face_descriptor fd : CGAL::faces_around_target(hd, pmesh))
  {
    if(fd == boost::graph_traits<PolygonMesh>::null_face())
      continue;

    ++incident_faces_n;
  }

  CGAL_assertion(incident_faces_n > 0);

  // Initial normal is just weighted by 1/n
  const FT initial_weight = FT(1)/incident_faces_n;

  for(face_descriptor fd : CGAL::faces_around_target(hd, pmesh))
  {
    if(fd == boost::graph_traits<PolygonMesh>::null_face())
      continue;

    normal = csv(normal, csclv(normalized_face_normals.at(fd), initial_weight));
    weights[face] = initial_weight;
  }
#else // angle-based initial weight
  for(face_descriptor fd : CGAL::faces_around_target(hd, pmesh))
  {
    if(fd == boost::graph_traits<PolygonMesh>::null_face())
      continue;

    halfedge_descriptor can_hd = halfedge(fd, pmesh);
    while(target(can_hd, pmesh) != vd)
      can_hd = next(can_hd, pmesh);

    FT weight = compute_angle_rad<GT>(pmesh.point(target(next(can_hd, pmesh), pmesh)),
                                      pmesh.point(target(can_hd, pmesh)),
                                      pmesh.point(source(can_hd, pmesh)));

// #define CGAL_MOST_VISIBLE_NORMAL_VERBOSE
#ifdef CGAL_MOST_VISIBLE_NORMAL_VERBOSE
      std::cout << "         weight: " << weight << std::endl;
#endif

    normal = csv(normal, csclv(normalized_face_normals.at(fd), weight));
    weights[fd] = weight;
  }
#endif

  CGAL::Polygon_mesh_processing::internal::normalize(normal, traits);

  // Iterative process to find the most visible
  bool converged = false;
  int iter = 0;
  FT old_sp = std::numeric_limits<FT>::infinity();
  for(;;)
  {
#ifdef CGAL_MOST_VISIBLE_NORMAL_VERBOSE
    std::cout << "Current best normal: " << normal << " for vertex: " << pmesh.point(vd) << std::endl;
#endif

    ++iter;

    // compute scaling coefficients
    boost::container::flat_map<face_descriptor, FT> alphas;
    alphas.reserve(16);
    FT alpha_sum = 0;

    for(face_descriptor fd : CGAL::faces_around_target(hd, pmesh))
    {
      if(fd == boost::graph_traits<PolygonMesh>::null_face())
        continue;

      FT normals_sp = sp(normal, normalized_face_normals.at(fd));
      if(normals_sp > 1) // everything is normalized so that shouldn't happen in theory, but it does in practice
        normals_sp = 1;
      if(normals_sp < -1)
        normals_sp = -1;

      const FT alpha = std::acos(normals_sp);

#ifdef CGAL_MOST_VISIBLE_NORMAL_VERBOSE
      std::cout << "         alpha: " << alpha << " (" << normalized_face_normals.at(fd) << ")" << std::endl;
#endif

      CGAL_assertion(alpha >= 0);
      alphas[fd] = alpha;
      alpha_sum += alpha;
    }

    if(alpha_sum == 0) // all the sp are 1 and thus all the vectors are identical
      return normal;

    CGAL_assertion(weights.size() == alphas.size());

    // compute the new weights
    FT weighted_sum = 0;

    for(face_descriptor fd : CGAL::faces_around_target(hd, pmesh))
    {
      if(fd == boost::graph_traits<PolygonMesh>::null_face())
        continue;

      FT& weight = weights[fd];
      weight *= alphas.at(fd) / alpha_sum;

#ifdef CGAL_MOST_VISIBLE_NORMAL_VERBOSE
      std::cout << "         weights[" << fd << "] = " << weight << std::endl;
#endif

      weighted_sum += weight;
    }

    if(weighted_sum == 0) // not too sure about that one... @fixme (correct if all weights are 0, but what otherwise?)
      return normal;

    for(face_descriptor fd : CGAL::faces_around_target(hd, pmesh))
    {
      if(fd == boost::graph_traits<PolygonMesh>::null_face())
        continue;

      weights[fd] /= weighted_sum;
    }

    // compute the new normal
    Vector_3 new_normal_base = CGAL::NULL_VECTOR;
    for(face_descriptor fd : CGAL::faces_around_target(hd, pmesh))
    {
      if(fd == boost::graph_traits<PolygonMesh>::null_face())
        continue;

#ifdef CGAL_MOST_VISIBLE_NORMAL_VERBOSE
      std::cout << "add: " << normalized_face_normals.at(fd) << " with weight: " << weights[fd] << std::endl;
#endif
      new_normal_base = csv(new_normal_base, csclv(normalized_face_normals.at(fd), weights[fd]));
    }

    CGAL::Polygon_mesh_processing::internal::normalize(new_normal_base, traits);

    // check convergence
    const FT bound = 1e-10;
    const FT new_sp = sp(normal, new_normal_base);
    converged = (CGAL::abs(new_sp - old_sp) < bound);

#ifdef CGAL_MOST_VISIBLE_NORMAL_VERBOSE
    std::cout << "new normal (base): " << new_normal_base << std::endl;
    std::cout << "   convergence ? " << old_sp << " " << new_sp << " diff: " << CGAL::abs(new_sp - old_sp) << std::endl;
#endif

    if(converged || iter > 1000000) // @fixme
    {
      normal = new_normal_base;
      break;
    }

    // a bit of relaxation
    const FT beta = 0.01;
    const Vector_3 beta_scaled_normal = csclv(new_normal_base, beta);
    normal = csv(beta_scaled_normal, csclv(normal, 1 - beta));
    old_sp = new_sp;
  }

#ifdef CGAL_MOST_VISIBLE_NORMAL_VERBOSE
  std::cout << "Converged in: " << iter << " iterations" << std::endl;
#endif

  return normal;
}

template<typename PolygonMesh, typename VertexNormalMap, typename NamedParameters>
void compute_most_visible_vertex_normals(VertexNormalMap& vertex_normal_map,
                                         const PolygonMesh& pmesh,
                                         const NamedParameters& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor          vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor            face_descriptor;

  CGAL_precondition(CGAL::is_triangle_mesh(pmesh));

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type  VPMap;
//  VPMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
//                             get_const_property_map(vertex_point, pmesh));

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type            GT;
  const GT traits = choose_param(get_param(np, internal_np::geom_traits), GT());

  typedef typename GT::Vector_3                                                 Vector;

  std::cout << "compute face normals" << std::endl;

  std::vector<Vector> normalized_face_normals(num_faces(pmesh));
  for(face_descriptor fd : faces(pmesh))
  {
    Vector en = compute_face_normal(fd, pmesh);
    CGAL::Polygon_mesh_processing::internal::normalize(en, traits);
    normalized_face_normals[fd] = en;
  }

  std::ofstream out("computed_normals.cgal.polylines.txt");

  // we start by evaluating the translation normal for each vertex
  for(vertex_descriptor vd : vertices(pmesh))
  {
#if 0
    Vector n = compute_vertex_normal_most_visible_with_optimization(vd, normalized_face_normals, pmesh, traits);
#else
    Vector n = compute_vertex_normal_most_visible_min_circle(vd, normalized_face_normals, pmesh, traits);
#endif

    CGAL::Polygon_mesh_processing::internal::normalize(n, traits);
    vertex_normal_map[vd] = n;

    out << "2 " << pmesh.point(vd) << " " << pmesh.point(vd) + 0.01 * n << std::endl;
  }
}

template<typename PolygonMesh, typename VertexNormalMap>
void compute_most_visible_vertex_normals(VertexNormalMap& vertex_normal_map, const PolygonMesh& pmesh)
{
  return compute_most_visible_vertex_normals(vertex_normal_map, pmesh, CGAL::parameters::all_default());
}

/**
* \ingroup PMP_normal_grp
* computes the unit normal at vertex `v` as the average of the normals of incident faces.
* @tparam PolygonMesh a model of `FaceGraph`
*
* @param v the vertex at which the normal is computed
* @param pmesh the polygon mesh containing `v`
* @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
*   If this parameter is omitted, an internal property map for
*   `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
*    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
* \cgalNamedParamsEnd
*
* @return the computed normal. The return type is a 3D vector type. It is
* either deduced from the `geom_traits` \ref pmp_namedparameters "Named Parameters" if provided,
* or the geometric traits class deduced from the point property map
* of `pmesh`.
*
* \warning This function involves a square root computation.
* If `Kernel::FT` does not have a `sqrt()` operation, the square root computation
* will be done approximately.
*/
template<typename PolygonMesh, typename NamedParameters>
#ifdef DOXYGEN_RUNNING
Vector_3
#else
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Vector_3
#endif
compute_vertex_normal(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
                      const PolygonMesh& pmesh,
                      const NamedParameters& np
                      )
{
  using boost::choose_param;

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
  typedef typename GT::Vector_3 Vector;
  GT traits = choose_param(get_param(np, internal_np::geom_traits), GT());

  typedef typename GetFaceNormalMap<PolygonMesh, NamedParameters>::NoMap DefaultMap;
  typedef typename boost::lookup_named_param_def <
    internal_np::face_normal_t,
    NamedParameters,
    DefaultMap> ::type FaceNormalMap;
  FaceNormalMap fnmap = choose_param(get_param(np, internal_np::face_normal), DefaultMap());
  bool fnmap_valid
    = !boost::is_same<FaceNormalMap,
                      DefaultMap
                     >::value;

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  Vector normal = traits.construct_vector_3_object()(CGAL::NULL_VECTOR);

  halfedge_descriptor he = halfedge(v, pmesh);
  // handle isolated vertices
  if (he==boost::graph_traits<PolygonMesh>::null_halfedge()) return normal;
  halfedge_descriptor end = he;
  do
  {
    if (!is_border(he, pmesh))
    {
      Vector n = fnmap_valid ? get(fnmap, face(he, pmesh))
                             : compute_face_normal(face(he, pmesh), pmesh, np);
      normal = traits.construct_sum_of_vectors_3_object()(normal, n);
    }
    he = opposite(next(he, pmesh), pmesh);
  } while (he != end);

  if ( ! traits.equal_3_object()(normal, CGAL::NULL_VECTOR))
    internal::normalize(normal, traits);
  return normal;
}

///\cond SKIP_IN_MANUAL

// compute_vertex_normal overloads
template <typename PolygonMesh>
typename GetGeomTraits<PolygonMesh>::type::Vector_3
compute_vertex_normal(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
                      const PolygonMesh& pmesh)
{
  return compute_vertex_normal(v, pmesh, CGAL::Polygon_mesh_processing::parameters::all_default());
}

/// \endcond

/**
* \ingroup PMP_normal_grp
* computes the outward unit vector normal for all vertices of the polygon mesh.
* @tparam PolygonMesh a model of `FaceListGraph`
* @tparam VertexNormalMap a model of `WritablePropertyMap` with
    `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and
    the return type of `compute_vertex_normal()` as value type.
*
* @param pmesh the polygon mesh
* @param vnm the property map in which the normals are written
* @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
*   If this parameter is omitted, an internal property map for
*   `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
*    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
* \cgalNamedParamsEnd
*
* \warning This function involves a square root computation.
* If `Kernel::FT` does not have a `sqrt()` operation, the square root computation
* will be done approximately.
*/
template <typename PolygonMesh
          , typename VertexNormalMap
          , typename NamedParameters
          >
void
compute_vertex_normals(const PolygonMesh& pmesh
                      , VertexNormalMap vnm
                      , const NamedParameters& np
                      )
{
  typedef typename GetGeomTraits<PolygonMesh,NamedParameters>::type Kernel;

  for(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v : vertices(pmesh)){
    typename Kernel::Vector_3 vec = compute_vertex_normal(v, pmesh, np);
    put(vnm, v, vec);
  }
}

///\cond SKIP_IN_MANUAL

// compute_vertex_normals overloads
template <typename PolygonMesh, typename VertexNormalMap>
void compute_vertex_normals(const PolygonMesh& pmesh, VertexNormalMap vnm)
{
  compute_vertex_normals(pmesh, vnm, CGAL::Polygon_mesh_processing::parameters::all_default());
}

/// \endcond

/**
* \ingroup PMP_normal_grp
* computes the outward unit vector normal for all vertices and faces of the polygon mesh.
* @tparam PolygonMesh a model of `FaceListGraph`
* @tparam VertexNormalMap a model of `WritablePropertyMap` with
    `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and
    `Kernel::Vector_3` as value type.
* @tparam FaceNormalMap a model of `ReadWritePropertyMap` with
    `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
    `Kernel::Vector_3` as value type.
*
* @param pmesh the polygon mesh
* @param vnm the property map in which the vertex normals are written
* @param fnm the property map in which the face normals are written
* @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
*   If this parameter is omitted, an internal property map for
*   `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
*    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
* \cgalNamedParamsEnd
*
* \warning This function involves a square root computation.
* If `Kernel::FT` does not have a `sqrt()` operation, the square root computation
* will be done approximately.
*/
template <typename PolygonMesh
          , typename VertexNormalMap
          , typename FaceNormalMap
          , typename NamedParameters
          >
void
compute_normals(const PolygonMesh& pmesh
                , VertexNormalMap vnm
                , FaceNormalMap fnm
                , const NamedParameters& np
                )
{
  compute_face_normals(pmesh, fnm, np);
  compute_vertex_normals(pmesh, vnm, np.face_normal_map(fnm));
}

///\cond SKIP_IN_MANUAL

// compute_normals overload
template <typename PolygonMesh, typename VertexNormalMap, typename FaceNormalMap>
void compute_normals(const PolygonMesh& pmesh,
                     VertexNormalMap vnm,
                     FaceNormalMap fnm)
{
  compute_normals(pmesh, vnm, fnm, CGAL::Polygon_mesh_processing::parameters::all_default());
}

/// \endcond

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H
