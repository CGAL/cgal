// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baskin Burak Senbaslar,
//                 Mael Rouxel-Labbé,
//                 Julian Komaromy

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PLANE_QUADRICS_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PLANE_QUADRICS_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_policy_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_common.h>

#include <Eigen/Dense>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

namespace plane_quadric_helpers {
  template<typename GeomTraits>
  typename Eigen::Matrix<typename GeomTraits::FT, 4, 4, Eigen::DontAlign> construct_quadric_from_normal(
      const typename GeomTraits::Vector_3& normal, 
      const typename GeomTraits::Point_3& point, 
      const GeomTraits& gt) 
  {
    typedef typename GeomTraits::FT FT;
    
    auto dot_product = gt.compute_scalar_product_3_object();
    auto construct_vector = gt.construct_vector_3_object();

    // negative dot product between the normal and the position vector
    const FT d = - dot_product(normal, construct_vector(ORIGIN, point));

    // row vector given by d appended to the normal
    const Eigen::Matrix<FT, 1, 4> row(normal.x(), normal.y(), normal.z(), d);

    // outer product
    return row.transpose() * row;
  }
}

// change the quadric construction for faces, but keep the standard plane quadric 
// method for (discontinous) edges
//
// this means we inherit from GarlandHeckbert_policies so we have the implementation 
// for discontinous edges and from the cost base class so that out implementation
// of face quadrics will be used
template<typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_plane_faces
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Point_3 Point_3;
  
  typedef typename Eigen::Matrix<FT, 4, 4, Eigen::DontAlign> Mat_4;

  public:
  template<typename VPM, typename TM> 
    Mat_4 construct_quadric_from_face(
        const VPM& point_map,
        const TM& tmesh, 
        typename boost::graph_traits<TM>::face_descriptor f, 
        const GeomTraits& gt) const
    {
      // this-> is used because we are calling a base class function
      const Vector_3 normal = common::construct_unit_normal_from_face<
        GeomTraits, VPM, TM>(point_map, tmesh, f, gt);

      // get any point of the face 
      const auto p = get(point_map, source(halfedge(f, tmesh), tmesh));

      return plane_quadric_helpers::construct_quadric_from_normal(normal, p, gt);
    }
};

template<typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_plane_edges
{
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename Eigen::Matrix<FT, 4, 4, Eigen::DontAlign> Mat_4;

  public:
  template<typename VPM, typename TM> 
  Mat_4 construct_quadric_from_edge(
      const VPM& point_map,
      const TM& tmesh, 
      typename boost::graph_traits<TM>::halfedge_descriptor he, 
      const GeomTraits& gt) const
  {
    const Vector_3 normal = construct_edge_normal<VPM, TM>(point_map, tmesh, he, gt);

    return plane_quadric_helpers::construct_quadric_from_normal(normal, 
        get(point_map, source(he, tmesh)), gt);
  }

  private:
  
  template<typename VPM, typename TM> 
  Vector_3 construct_edge_normal(
      const VPM& point_map,
      const TM& tmesh, 
      typename boost::graph_traits<TM>::halfedge_descriptor he, 
      const GeomTraits& gt) const
  {
    typedef typename GeomTraits::Vector_3 Vector_3;

    typedef typename boost::graph_traits<TM>::vertex_descriptor vertex_descriptor;
    
    // TODO we do a potentially redundant calculation here, as we have 
    // almost certainly already calculated this when constructing a quadric for the face
    const Vector_3 face_normal = common::construct_unit_normal_from_face<GeomTraits, VPM, TM>(
        point_map, tmesh, face(he, tmesh), gt);

    const vertex_descriptor vs = source(he, tmesh);
    const vertex_descriptor vt = target(he, tmesh);

    const Vector_3 edge_vector = Vector_3(get(point_map, vs), get(point_map, vt));
    const Vector_3 discontinuity_normal = cross_product(edge_vector, face_normal);

    // normalize
    const Vector_3 normal = discontinuity_normal 
      / sqrt(discontinuity_normal.squared_length());

    return normal;
  }
};
} //namespace internal
} //namespace Surface_mesh_simplification
} //namespace CGAL

#endif //CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PLANE_QUADRICS_H
