// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baskin Burak Senbaslar,
//                 Mael Rouxel-Labbé
//

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_INTERNAL_GARLAND_HECKBERT_CORE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_INTERNAL_GARLAND_HECKBERT_CORE_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_quadrics.h>


#include <Eigen/Dense>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

template<typename GT_>
struct GarlandHeckbert_matrix_type
{
  typedef typename GT_::FT                                                      FT;
  typedef Eigen::Matrix<FT, 4, 4, Eigen::DontAlign>                             type;
};

template<class TM_, class VPM_, typename GT_>
struct GarlandHeckbert_core
{
  typedef TM_                                                                   Triangle_mesh;
  typedef boost::graph_traits<Triangle_mesh>                                    GraphTraits;
  typedef typename GraphTraits::vertex_descriptor                               vertex_descriptor;
  typedef typename GraphTraits::halfedge_descriptor                             halfedge_descriptor;
  typedef typename GraphTraits::face_descriptor                                 face_descriptor;

  typedef VPM_                                                                  Vertex_point_pmap;
  typedef typename boost::property_traits<Vertex_point_pmap>::value_type        Point;
  typedef typename boost::property_traits<Vertex_point_pmap>::reference         Point_reference;

  typedef GT_                                                                   Geom_traits;
  typedef typename Geom_traits::FT                                              FT;
  typedef typename Geom_traits::Line_3                                          Line_3;
  typedef typename Geom_traits::Plane_3                                         Plane_3;
  typedef typename Geom_traits::Vector_3                                        Vector_3;
  typedef typename Geom_traits::Point_3                                         Point_3;

  typedef typename GarlandHeckbert_matrix_type<GT_>::type                       Matrix4x4;
  typedef Eigen::Matrix<FT, 1, 4>                                               Row4;
  typedef Eigen::Matrix<FT, 4, 1>                                               Col4;
  
  static Col4 point_to_homogenous_column(const Point& pt)
  {
    return Col4(pt.x(), pt.y(), pt.z(), FT(1));
  }
  /**
  * Combines two Q matrices.
  * It is simply the addition of two matrices
  */
  inline static Matrix4x4 combine_matrices(const Matrix4x4& aFirst, const Matrix4x4& aSecond)
  {
    return aFirst + aSecond;
  }

  /**
  * Returns `true` if the target of h is a discontinuity vertex in tmesh
  *
  * Currently, only checks if vertex belongs to a border.
  */
  static bool is_discontinuity_vertex(const halfedge_descriptor h, const Triangle_mesh& tmesh)
  {
    return is_border(target(h, tmesh), tmesh);
  }

  static bool is_discontinuity_edge(const halfedge_descriptor h, const Triangle_mesh& tmesh)
  {
    return is_border_edge(h, tmesh);
  }

  /*
   * TODO build selection of quadric type into the policy system
   */
  template <typename VCM>
  static void fundamental_error_quadrics(VCM& vcm, // quadrics container
                                         const Triangle_mesh& tmesh,
                                         const Vertex_point_pmap& vpm,
                                         const Geom_traits& gt,
                                         const FT discontinuity_multiplier = FT(100))
  {
    Matrix4x4 nq = Matrix4x4::Zero();
    for(vertex_descriptor v : vertices(tmesh))
      put(vcm, v, nq);

    for(face_descriptor f : faces(tmesh))
    {
      if (f == GraphTraits::null_face())
      {
        continue;
      }
      
      const halfedge_descriptor h = halfedge(f, tmesh);
      
      // get the three vertices of the triangular face
      const Point_3 p = get(vpm, source(h, tmesh));
      const Point_3 q = get(vpm, target(h, tmesh));
      const Point_3 r = get(vpm, target(next(h, tmesh), tmesh));
      
      // construtct the (4 x 4) matrix representing the plane quadric
      // TODO find good default values for std deviations
      const auto plane_matrix 
        = quadrics::construct_classical_plane_quadric<Geom_traits>(p, q, r, gt);

      for(halfedge_descriptor shd : halfedges_around_face(h, tmesh))
      {
        const vertex_descriptor vs = source(shd, tmesh);
        const vertex_descriptor vt = target(shd, tmesh);

        put(vcm, vt, combine_matrices(get(vcm, vt), plane_matrix));

        if(!is_discontinuity_edge(shd, tmesh))
          continue;

        //TODO optimize and clean up the computation of this quadric
        // some functions may be called unnecessarily often here (i.e. too 
        // much normalization and two cross products)
        //
        // also, behavious at discontinuous edges are largely untested at this stage
        const auto edge_source = get(vpm, vs);
        const auto plane_n = gt.construct_unit_normal_3_object()(r, p, q);
        // vector along the current discontinuous edge 
        const auto v_edge = gt.construct_vector_3_object()(get(vpm, vs), get(vpm, vt));
        // pstp x plane_n
        const auto disc_normal = gt.construct_cross_product_vector_3_object()(v_edge, plane_n);

        const FT disc_length = CGAL::sqrt(gt.compute_squared_length_3_object()(disc_normal));

        const auto discontinuous_matrix = discontinuity_multiplier 
          * quadrics::construct_classical_plane_quadric<Geom_traits>(
              disc_normal / disc_length, edge_source, gt);

        put(vcm, vs, combine_matrices(get(vcm, vs), discontinuous_matrix));
        put(vcm, vt, combine_matrices(get(vcm, vt), discontinuous_matrix));
      }
    }
  }


  /*
  * Return the point p that minimizes p' Q p where p is free.
  * p0, and p1 are the points that are being collapsed.
  * aQuadric is the matrix that is the combination of matrices
  * of p0 and p1.
  */
  static Col4 optimal_point(const Matrix4x4& aQuadric,
                            const Col4& p0,
                            const Col4& p1)
  {
    Matrix4x4 X;
    X << aQuadric.block(0, 0, 3, 4), 0, 0, 0, 1;

    Col4 opt_pt;

    if(X.determinant() == 0)
    {
      // not invertible
      const Col4 p1mp0 = std::move(p1 - p0);
      const FT a = (p1mp0.transpose() * aQuadric * p1mp0)(0, 0);
      const FT b = 2 * (p0.transpose() * aQuadric * p1mp0)(0, 0);

      if(a == 0)
      {
        if(b < 0)
          opt_pt = p1;
        else if(b == 0)
          opt_pt = 0.5 * (p0 + p1);
        else
          opt_pt = p0;
      }
      else
      {
        FT ext_t = -b/(2*a);
        if(ext_t < 0 || ext_t > 1 || a < 0)
        {
          // one of endpoints
          FT p0_cost = (p0.transpose() * aQuadric * p0)(0, 0);
          FT p1_cost = (p1.transpose() * aQuadric * p1)(0, 0);
          if(p0_cost > p1_cost)
            opt_pt = p1;
          else
            opt_pt = p0;
        }
        else
        {
          // extremum of the parabola
          opt_pt = p0 + ext_t * (p1 - p0);
        }
      }
    }
    else // invertible
    {
      opt_pt = X.inverse().col(3); // == X.inverse() * (0 0 0 1)
    }
    return opt_pt;
  }
};

} // namespace internal
} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_INTERNAL_GARLAND_HECKBERT_CORE_H
