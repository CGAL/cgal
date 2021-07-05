// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé,
//                 Julian Komaromy

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_TRIANGLE_QUADRICS_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_TRIANGLE_QUADRICS_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>

#include <Eigen/Dense>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {
// change the quadric construction for faces, but keep the standard plane quadric 
// method for (discontinous) edges
//
// this means we inherit from GarlandHeckbert_policies so we have the implementation 
// for discontinous edges and from the cost base class so that out implementation
// of face quadrics will be used
template<typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_triangle_faces
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename Eigen::Matrix<FT, 4, 4, Eigen::DontAlign> Mat_4;
  
  public: 
  template<typename VPM, typename TM> 
  Mat_4 construct_quadric_from_face(
      const VPM& point_map,
      const TM& tmesh, 
      typename boost::graph_traits<TM>::face_descriptor f, 
      const GeomTraits& gt) const
  {
    auto construct_vector = gt.construct_vector_3_object();
    auto cross_product = gt.construct_cross_product_vector_3_object();
    auto sum_vectors = gt.construct_sum_of_vectors_3_object();
    auto dot_product = gt.compute_scalar_product_3_object();

    typedef typename boost::property_traits<VPM>::reference Point_reference;
    typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;

    const halfedge_descriptor h = halfedge(f, tmesh);

    // get all points and turn them into location vectors so we can use cross product on them
    const Point_reference p = get(point_map, source(h, tmesh));
    const Point_reference q = get(point_map, target(h, tmesh));
    const Point_reference r = get(point_map, target(next(h, tmesh), tmesh));

    Vector_3 a = construct_vector(ORIGIN, p);
    Vector_3 b = construct_vector(ORIGIN, q);
    Vector_3 c = construct_vector(ORIGIN, r);

    const Vector_3 ab = cross_product(a, b);
    const Vector_3 bc = cross_product(b, c);
    const Vector_3 ca = cross_product(c, a);

    const Vector_3 sum_of_cross_product = sum_vectors(sum_vectors(ab, bc), ca);
    const FT scalar_triple_product = dot_product(ab, c);

    Eigen::Matrix<FT, 1, 4> row;
    
    row << sum_of_cross_product.x(), sum_of_cross_product.y(), 
        sum_of_cross_product.z(), -scalar_triple_product;

    // calculate the outer product of row^t*row
    return row.transpose() * row;
  }
};

} //namespace internal
} //namespace Surface_mesh_simplification
} //namespace CGAL

#endif //CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_POLICIES_H
