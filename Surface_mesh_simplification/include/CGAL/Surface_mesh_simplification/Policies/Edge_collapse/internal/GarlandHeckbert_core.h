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
    if(is_border(target(h, tmesh), tmesh))
      return true;

    return false;
  }

  static bool is_discontinuity_edge(const halfedge_descriptor h, const Triangle_mesh& tmesh)
  {
    return is_border_edge(h, tmesh);
  }

  /*
  * fundamental error quadric for the target vertex of h in tmesh
  * Unused, but leaving it as it might still be useful somehow
  */
  static Matrix4x4 fundamental_error_quadric(const halfedge_descriptor h,
                                             const Triangle_mesh& tmesh,
                                             const Vertex_point_pmap& vpm,
                                             const Geom_traits& gt,
                                             const FT discontinuity_multiplier = FT(100))
  {
    Matrix4x4 quadric = Matrix4x4::Zero();

    const vertex_descriptor target_vd = target(h, tmesh);
    const Vector_3 target_vertex_vector = gt.construct_vector_3_object()(CGAL::ORIGIN, get(vpm, target_vd));

    // const bool discontinuity_vertex = is_discontinuity_vertex(h, tmesh);

    for(const halfedge_descriptor hd : halfedges_around_target(target_vd, tmesh))
    {
      const face_descriptor fd = face(hd, tmesh);
      if(fd == GraphTraits::null_face())
        continue;

      Plane_3 plane(get(vpm, source(hd, tmesh)),
                    get(vpm, target(hd, tmesh)),
                    get(vpm, target(next(hd, tmesh), tmesh)));

      Row4 plane_mtr;
      const FT norm = CGAL::sqrt(CGAL::square(plane.a()) +
                                 CGAL::square(plane.b()) +
                                 CGAL::square(plane.c()));
      const FT den = FT(1) / norm;

      plane_mtr << den * plane.a(),
                   den * plane.b(),
                   den * plane.c(),
                   den * plane.d();
      quadric += plane_mtr.transpose() * plane_mtr;

      if(is_discontinuity_edge(hd, tmesh))
      {
        const vertex_descriptor source_vd = source(hd, tmesh);

        const Vector_3 p1p2 = gt.construct_vector_3_object()(get(vpm, source_vd), get(vpm, target_vd));
        const Vector_3 normal = gt.construct_cross_product_vector_3_object()(
                                  p1p2, gt.construct_orthogonal_vector_3_object()(plane));

        const FT d = - normal * target_vertex_vector;
        const FT norm = CGAL::sqrt(gt.compute_squared_length_3_object()(normal));
        const FT den = FT(1) / norm;

        plane_mtr << den * normal.x(),
                     den * normal.y(),
                     den * normal.z(),
                     den * d;
        quadric += discontinuity_multiplier * plane_mtr.transpose() * plane_mtr;
      }

      const halfedge_descriptor shd = next(hd, tmesh);
      if(is_discontinuity_edge(shd, tmesh))
      {
        const Vector_3 p1p2 = gt.construct_vector_3_object()(get(vpm, target_vd), get(vpm, target(shd, tmesh)));
        const Vector_3 normal = gt.construct_cross_product_vector_3_object()(
                                  p1p2, gt.construct_orthogonal_vector_3_object()(plane));

        const FT d = - normal * target_vertex_vector;
        const FT norm = CGAL::sqrt(gt.compute_squared_length_3_object()(normal));
        const FT den = FT(1) / norm;

        plane_mtr << den * normal.x(),
                     den * normal.y(),
                     den * normal.z(),
                     den * d;
        quadric += discontinuity_multiplier * plane_mtr.transpose() * plane_mtr;
      }
    }

    return quadric;
  }

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
      const halfedge_descriptor h = halfedge(f, tmesh);

      const Point_reference p = get(vpm, source(h, tmesh));
      const Point_reference q = get(vpm, target(h, tmesh));
      const Point_reference r = get(vpm, target(next(h, tmesh), tmesh));

      const FT rpx = p.x() - r.x();
      const FT rpy = p.y() - r.y();
      const FT rpz = p.z() - r.z();
      const FT rqx = q.x() - r.x();
      const FT rqy = q.y() - r.y();
      const FT rqz = q.z() - r.z();

      // Cross product rp * rq
      const FT a = rpy*rqz - rqy*rpz;
      const FT b = rpz*rqx - rqz*rpx;
      const FT c = rpx*rqy - rqx*rpy;
      const FT d = - a*r.x() - b*r.y() - c*r.z();

      const Vector_3 plane_n = gt.construct_vector_3_object()(a, b, c);
      const FT norm = CGAL::sqrt(CGAL::square(a) + CGAL::square(b) + CGAL::square(c));
      const FT den = FT(1) / norm;

      Row4 plane_mtr_r;
      plane_mtr_r << den * a, den * b, den * c, den * d;
      const Matrix4x4 plane_mtr = plane_mtr_r.transpose() * plane_mtr_r;

      for(halfedge_descriptor shd : halfedges_around_face(h, tmesh))
      {
        const vertex_descriptor vs = source(shd, tmesh);
        const vertex_descriptor vt = target(shd, tmesh);

        put(vcm, vt, combine_matrices(get(vcm, vt), plane_mtr));

        if(!is_discontinuity_edge(shd, tmesh))
          continue;

        const Vector_3 pspt = gt.construct_vector_3_object()(get(vpm, vs), get(vpm, vt));
        const Vector_3 disc_plane_n = gt.construct_cross_product_vector_3_object()(pspt, plane_n);

        // the plane contains the edge, so taking 'vs' or 'vt' will yield the same 'd'
        const Vector_3 vvt = gt.construct_vector_3_object()(CGAL::ORIGIN, get(vpm, vt));
        const FT disc_d = - gt.compute_scalar_product_3_object()(disc_plane_n, vvt);

        const FT disc_norm = CGAL::sqrt(gt.compute_squared_length_3_object()(disc_plane_n));
        const FT disc_den = FT(1) / disc_norm;

        Row4 disc_mtr_r;
        disc_mtr_r << disc_den * disc_plane_n.x(),
                      disc_den * disc_plane_n.y(),
                      disc_den * disc_plane_n.z(),
                      disc_den * disc_d;

        const Matrix4x4 disc_mtr = discontinuity_multiplier * disc_mtr_r.transpose() * disc_mtr_r;

        put(vcm, vs, combine_matrices(get(vcm, vs), disc_mtr));
        put(vcm, vt, combine_matrices(get(vcm, vt), disc_mtr));
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
