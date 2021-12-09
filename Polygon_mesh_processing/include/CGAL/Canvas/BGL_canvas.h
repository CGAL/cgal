// Copyright (c) 2021 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_CANVAS_CANVAS_BGL_CANVAS_H
#define CGAL_CANVAS_CANVAS_BGL_CANVAS_H

#include <CGAL/Canvas/Base_canvas.h>
#include <CGAL/Canvas/BGL_canvas_point.h>
#include <CGAL/Canvas/Metric_field/Euclidean_metric_field.h>

#include <CGAL/assertions.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <Eigen/Dense>
#endif

#include <cstddef>
#include <deque>
#include <istream>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <utility>

namespace CGAL {
namespace Canvas {

template<typename TriangleMesh, typename VertexPointMap, typename Metric_field, typename GeomTraits>
class BGL_canvas
  : public Base_canvas<GeomTraits,
                       BGL_canvas_point<TriangleMesh, VertexPointMap, Metric_field, GeomTraits>,
                       Metric_field>
{
private:
  using Self = BGL_canvas<TriangleMesh, VertexPointMap, Metric_field, GeomTraits>;

public:
  using Canvas_point = BGL_canvas_point<TriangleMesh, VertexPointMap, Metric_field, GeomTraits>;
  using Canvas_point_handle = Canvas_point*;

  using VPM = VertexPointMap;
  using Geom_traits = GeomTraits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

  using Base = Base_canvas<Geom_traits, Canvas_point, Metric_field>;

  using Metric = typename Base::Metric;
  using Vector3d = typename Base::Vector3d;

  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
  using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

  using HA_Tag = CGAL::dynamic_halfedge_property_t<FT>;
  using Sector_angle_map = typename boost::property_map<TriangleMesh, HA_Tag>::const_type;
  using EN_Tag = CGAL::dynamic_edge_property_t<Vector3d>;
  using Edge_normal_map = typename boost::property_map<TriangleMesh, EN_Tag>::const_type;
  using FN_Tag = CGAL::dynamic_face_property_t<Vector3d>;
  using Face_normal_map = typename boost::property_map<TriangleMesh, FN_Tag>::const_type;

  using VIM = typename CGAL::GetInitializedVertexIndexMap<TriangleMesh>::const_type;
  using VID = typename boost::property_traits<VIM>::value_type;

protected:
  const TriangleMesh& m_tmesh;
  VPM m_vpm;
  VIM m_vim;

  Sector_angle_map m_angles;
  Edge_normal_map m_edge_normals;
  Face_normal_map m_face_normals;

public:
  BGL_canvas(const TriangleMesh& tmesh,
             const VertexPointMap vpm,
             const Metric_field* mf,
             const Geom_traits& gt = Geom_traits())
    :
      Base(mf, gt),
      m_tmesh(tmesh),
      m_vpm(vpm),
      m_vim(CGAL::get_initialized_vertex_index_map(m_tmesh)),
      m_angles(get(HA_Tag(), m_tmesh)),
      m_edge_normals(get(EN_Tag(), m_tmesh)),
      m_face_normals(get(FN_Tag(), m_tmesh))
  {
  }

public:
  const TriangleMesh& mesh() const { return m_tmesh; }
  VID index(const vertex_descriptor v) const { return get(m_vim, v); }
  FT angle(const halfedge_descriptor h) const { return get(m_angles, h); }
  Vector3d normal(const edge_descriptor e) const { return get(m_edge_normals, e); }
  Vector3d normal(const face_descriptor f) const { return get(m_face_normals, f); }

protected:
  FT compute_sec_angle(vertex_descriptor vh0, vertex_descriptor vh1, vertex_descriptor vh2) const
  {
    // compute the angle between v0v1 and v1v2
    Vector_3 v1(m_tmesh.point(vh0), m_tmesh.point(vh1));
    Vector_3 v2(m_tmesh.point(vh0), m_tmesh.point(vh2));

    return this->m_gt.compute_approximate_angle_3_object()(v1, v2);
  }

  void initialize_cumulative_angles()
  {
    for(vertex_descriptor v : vertices(m_tmesh))
    {
      const std::size_t id = get(m_vim, v);
      const Canvas_point& cp = this->m_canvas_points[id];

      FT angle = 0;

      // compute the cumulative angles in a star (the initial one is 0)
      // the entry in the map is of the form :
      // [center C of the star, vertex V on the first ring] = angle V, C, Vnext

      for(halfedge_descriptor h : halfedges_around_target(halfedge(v, m_tmesh), m_tmesh))
      {
        if(is_border(h, m_tmesh))
          continue;

        vertex_descriptor vh = source(h, m_tmesh);
        vertex_descriptor vh_next = target(next(h, m_tmesh), m_tmesh);

        // set up the angle for vh, v, vh_next in the angle memory
        put(m_angles, h, angle);

        FT loc_angle = compute_sec_angle(v, vh, vh_next);
        angle += loc_angle;
      }

      // normalize the angles (sum must be 2*pi)
      FT norm_coeff = 2. * CGAL_PI / angle;
      for(halfedge_descriptor h : halfedges_around_target(halfedge(v, m_tmesh), m_tmesh))
        put(m_angles, h, norm_coeff * get(m_angles, h));
    }
  }

  void add_to_edge_normals_map(edge_descriptor e,
                               const Vector3d& normal)
  {
    // This is fine because there is 0, 1, or 2 insertions
    // - On the 1st insertion, get is default constructed (null vector) and the normalize is just superfluous
    // - On the 2nd insertion, it gives the normalized average, since 'normal' is normalized
    Vector3d& n = get(m_edge_normals, e);
    n += normal;
    n /= n.norm();
    put(m_edge_normals, e, n);
  }

  void compute_normals()
  {
    for(face_descriptor f : faces(m_tmesh))
    {
      Vector_3 n = CGAL::Polygon_mesh_processing::compute_face_normal(f, m_tmesh);
      CGAL::Polygon_mesh_processing::internal::normalize(n, this->m_gt);
      Vector3d en(n[0], n[1], n[2]);
      put(m_face_normals, f, en);

      // add this face normal to its edges
      for(halfedge_descriptor h : halfedges_around_face(f, m_tmesh))
        add_to_edge_normals_map(edge(h, m_tmesh), en);
    }
  }

public:
  void initialize() override
  {
#if (VERBOSITY > 5)
    std::cout << "canvas initialization" << std::endl;
#endif

    // build the canvas points
    for(vertex_descriptor v : vertices(m_tmesh))
      this->m_canvas_points.emplace_back(get(m_vim, v), get(m_vpm, v), v, *this);

#if (VERBOSITY > 5)
    std::cout << this->m_canvas_points.size() << " points on canvas" << std::endl;
#endif

    initialize_cumulative_angles();
    compute_normals();
  }

  template <typename DistanceMap>
  void paint(DistanceMap distances)
  {
    Base::paint();

    for(Canvas_point& cp : this->m_canvas_points)
      put(distances, cp.vertex(), cp.distance_to_closest_seed());
  }

  template <typename VertexSet>
  void initialize_seeds(const VertexSet& seeds)
  {
#if (VERBOSITY > 0)
    std::cout << "Initialize " << seeds.size() << " seeds" << std::endl;
#endif

    int seed_id = 0;
    for(vertex_descriptor v : seeds)
    {
      const Point_3& p = m_tmesh.point(v);
      this->m_seeds.insert_new_seed(p[0], p[1], p[2]);
      Canvas_point& cp = this->m_canvas_points[get(m_vim, v)];
      Base::initialize_canvas_point(cp, 0. /*distance*/, seed_id++);
    }
  }
};

} // namespace Canvas
} // namespace CGAL

#endif // CGAL_CANVAS_CANVAS_BGL_CANVAS_H
