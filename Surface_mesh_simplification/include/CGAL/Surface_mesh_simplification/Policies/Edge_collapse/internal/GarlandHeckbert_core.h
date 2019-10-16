// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
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

template<class TM_, class VPM_>
struct GarlandHeckbert_core
{
  typedef TM_                                                                   TM;
  typedef boost::graph_traits<TM>                                               GraphTraits;
  typedef typename GraphTraits::vertex_descriptor                               vertex_descriptor;
  typedef typename GraphTraits::halfedge_descriptor                             halfedge_descriptor;
  typedef typename GraphTraits::face_descriptor                                 face_descriptor;

  typedef VPM_                                                                  Vertex_point_pmap;
  typedef typename boost::property_traits<Vertex_point_pmap>::value_type        Point;
  typedef typename boost::property_traits<Vertex_point_pmap>::reference         Point_reference;

  typedef typename Kernel_traits<Point>::Kernel                                 Kernel;
  typedef typename Kernel::FT                                                   FT;
  typedef typename Kernel::RT                                                   RT;
  typedef typename Kernel::Line_3                                               Line_3;
  typedef typename Kernel::Plane_3                                              Plane_3;
  typedef typename Kernel::Vector_3                                             Vector_3;

  typedef Eigen::Matrix<FT, 4, 4>                                               Matrix4x4;
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
  * Returns `true` if the target of aHD is a discontinuity vertex in aTM
  *
  * Currently, only checks if vertex belongs to a border.
  */
  static bool is_discontinuity_vertex(const halfedge_descriptor aHD, const TM& aTM)
  {
    if(is_border(target(aHD, aTM), aTM))
      return true;

    return false;
  }

  static bool is_discontinuity_edge(const halfedge_descriptor aHD, const TM& aTM)
  {
    return is_border_edge(aHD, aTM);
  }

  /*
  * fundamental error quadric for the target vertex of aHD in aTM
  * Unused, but leaving it as it might still be useful somehow
  */
  static Matrix4x4 fundamental_error_quadric(const halfedge_descriptor aHD,
                                             const TM& aTM,
                                             const Vertex_point_pmap& aVPM,
                                             const FT aDiscontinuityMultiplier = FT(100))
  {
    Matrix4x4 quadric = Matrix4x4::Zero();

    const vertex_descriptor target_vd = target(aHD, aTM);
    const Vector_3 target_vertex_vector(CGAL::ORIGIN, get(aVPM, target_vd));

    // const bool discontinuity_vertex = is_discontinuity_vertex(aHD, aTM);

    for(const halfedge_descriptor hd : halfedges_around_target(target_vd, aTM))
    {
      const face_descriptor fd = face(hd, aTM);
      if(fd == GraphTraits::null_face())
        continue;

      Plane_3 plane(get(aVPM, source(hd, aTM)),
                    get(aVPM, target(hd, aTM)),
                    get(aVPM, target(next(hd, aTM), aTM)));

      Row4 plane_mtr;
      const FT norm = sqrt(CGAL::square(plane.a()) +
                           CGAL::square(plane.b()) +
                           CGAL::square(plane.c()));
      const FT den = FT(1) / norm;

      plane_mtr << den * plane.a(),
                   den * plane.b(),
                   den * plane.c(),
                   den * plane.d();
      quadric += plane_mtr.transpose() * plane_mtr;

      if(is_discontinuity_edge(hd, aTM))
      {
        const vertex_descriptor source_vd = source(hd, aTM);
        const Vector_3 p1p2(get(aVPM, source_vd), get(aVPM, target_vd));
        const Vector_3 normal = cross_product(p1p2, plane.orthogonal_vector());

        const FT d = - normal * target_vertex_vector;
        const FT norm = sqrt(CGAL::square(normal.x()) +
                             CGAL::square(normal.y()) +
                             CGAL::square(normal.z()));
        const FT den = FT(1) / norm;

        plane_mtr << den * normal.x(),
                     den * normal.y(),
                     den * normal.z(),
                     den * d;
        quadric += aDiscontinuityMultiplier * plane_mtr.transpose() * plane_mtr;
      }

      const halfedge_descriptor shd = next(hd, aTM);
      if(is_discontinuity_edge(shd, aTM))
      {
        const Vector_3 p1p2(get(aVPM, target_vd), get(aVPM, target(shd, aTM)));
        const Vector_3 normal = cross_product(p1p2, plane.orthogonal_vector());

        const FT d = - normal * target_vertex_vector;
        const FT norm = sqrt(CGAL::square(normal.x()) +
                             CGAL::square(normal.y()) +
                             CGAL::square(normal.z()));
        const FT den = FT(1) / norm;

        plane_mtr << den * normal.x(),
                     den * normal.y(),
                     den * normal.z(),
                     den * d;
        quadric += aDiscontinuityMultiplier * plane_mtr.transpose() * plane_mtr;
      }
    }

    return quadric;
  }

  template <typename VCM>
  static void fundamental_error_quadrics(VCM& vcm, // quadrics container
                                         const TM& aTM,
                                         const Vertex_point_pmap& aVPM,
                                         const FT aDiscontinuityMultiplier = FT(100))
  {
    Matrix4x4 nq = Eigen::Matrix<FT, 4, 4>::Zero();
    for(vertex_descriptor v : vertices(aTM))
      put(vcm, v, nq); // @todo necessary ?

    for(face_descriptor f : faces(aTM))
    {
      const halfedge_descriptor h = halfedge(f, aTM);

      const Point_reference p = get(aVPM, source(h, aTM));
      const Point_reference q = get(aVPM, target(h, aTM));
      const Point_reference r = get(aVPM, target(next(h, aTM), aTM));

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

      const Vector_3 plane_n(a, b, c);
      const FT norm = CGAL::sqrt(CGAL::square(a) + CGAL::square(b) + CGAL::square(c));
      const FT den = FT(1) / norm;

      Row4 plane_mtr_r;
      plane_mtr_r << den * a, den * b, den * c, den * d;
      const Matrix4x4 plane_mtr = plane_mtr_r.transpose() * plane_mtr_r;

      for(halfedge_descriptor shd : halfedges_around_face(h, aTM))
      {
        const vertex_descriptor vs = source(shd, aTM);
        const vertex_descriptor vt = target(shd, aTM);

        put(vcm, vt, combine_matrices(get(vcm, vt), plane_mtr));

        if(!is_discontinuity_edge(shd, aTM))
          continue;

        const Vector_3 pspt(get(aVPM, vs), get(aVPM, vt));
        const Vector_3 disc_plane_n = cross_product(pspt, plane_n);

        // the plane contains the edge, so taking 'vs' or 'vt' will yield the same 'd'
        const FT disc_d = - disc_plane_n * Vector_3(CGAL::ORIGIN, get(aVPM, vt));

        const FT disc_norm = CGAL::sqrt(CGAL::square(disc_plane_n.x()) +
                                        CGAL::square(disc_plane_n.y()) +
                                        CGAL::square(disc_plane_n.z()));
        const FT disc_den = FT(1) / disc_norm;

        Row4 disc_mtr_r;
        disc_mtr_r << disc_den * disc_plane_n.x(),
                      disc_den * disc_plane_n.y(),
                      disc_den * disc_plane_n.z(),
                      disc_den * d;

        const Matrix4x4 disc_mtr = aDiscontinuityMultiplier * disc_mtr_r.transpose() * disc_mtr_r;

        put(vcm, vs, combine_matrices(get(vcm, vs), disc_mtr));
        put(vcm, vt, combine_matrices(get(vcm, vt), disc_mtr));
      }
    }
  }


  /*
  * Return the point p that minimizes p' Q p where p is free.
  * aP0, and aP1 are the points that are being collapsed.
  * aQuadric is the matrix that is the combination of matrices
  * of aP0 and aP1.
  */
  static Col4 optimal_point(const Matrix4x4& aQuadric,
                            const Col4& aP0,
                            const Col4& aP1)
  {
    Matrix4x4 X;
    X << aQuadric.block(0, 0, 3, 4), 0, 0, 0, 1;

    Col4 opt_pt;

    if(X.determinant() == 0)
    {
      // not invertible
      Col4 p1mp0 = std::move(aP1 - aP0);
      const FT a = (p1mp0.transpose() * aQuadric * p1mp0)(0, 0);
      const FT b = 2 * (aP0.transpose() * aQuadric * p1mp0)(0, 0);

      if(a == 0)
      {
        if(b < 0)
          opt_pt = aP1;
        else if(b == 0)
          opt_pt = 0.5 * (aP0 + aP1);
        else
          opt_pt = aP0;
      }
      else
      {
        FT ext_t = -b/(2*a);
        if(ext_t < 0 || ext_t > 1 || a < 0)
        {
          // one of endpoints
          FT aP0_cost = (aP0.transpose() * aQuadric * aP0)(0, 0);
          FT aP1_cost = (aP1.transpose() * aQuadric * aP1)(0, 0);
          if(aP0_cost > aP1_cost)
            opt_pt = aP1;
          else
            opt_pt = aP0;
        }
        else
        {
          // extremum of the parabola
          opt_pt = aP0 + ext_t * (aP1 - aP0);
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
