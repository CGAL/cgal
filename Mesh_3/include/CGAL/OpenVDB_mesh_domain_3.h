// Copyright (c) 2025 GeometryFactory Sarl (France)
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

#ifndef CGAL_MESH_3_OPENVDB_MESH_DOMAIN_3_H
#define CGAL_MESH_3_OPENVDB_MESH_DOMAIN_3_H

#include <CGAL/point_generators_3.h>
#include <CGAL/Bbox_3.h>

#include <openvdb/openvdb.h>
#include <openvdb/tools/RayIntersector.h>

#include <variant>
#include <utility>
#include <optional>

namespace CGAL
{
template<typename GeomTraits>
class OpenVDB_mesh_domain_3
{
public:
  using R = GeomTraits;
  using Point_3 = typename R::Point_3;
  using Segment_3 = typename R::Segment_3;

  using Subdomain_index = int;
  using Surface_patch_index = std::pair<Subdomain_index, Subdomain_index>;
  using Index = std::variant<Subdomain_index, Surface_patch_index>;

  using Intersection = std::tuple<Point_3, Index, int>;

  using Has_features = CGAL::Tag_false;

private:
  using Sphere_3 = typename GeomTraits::Sphere_3;
  using FT = typename GeomTraits::FT;
  using Point_creator = typename CGAL::Creator_uniform_3<FT, Point_3>;

public:
  OpenVDB_mesh_domain_3(openvdb::FloatGrid::Ptr grid) :
    _grid(grid),
    _intersector(*grid, 0.0),
    _boundingSphere(Point_3(0.0, 0.0, 0.0), 100.0 * 100.0)
  {
    m_bbox = _boundingSphere.bbox();
  }

  const CGAL::Bbox_3& bbox() const { return m_bbox; }

private:
  const Sphere_3& boundingSphere() const { return _boundingSphere; }

  static const Point_3& point(const Intersection& intersection) { return std::get<0>(intersection); }
  static const Index& index(const Intersection& intersection) { return std::get<1>(intersection); }
  static const int dimension(const Intersection& intersection) { return std::get<2>(intersection); }

  bool intersects(const Point_3 &from, const Point_3 &to, Point_3 &p) const
  {
    const openvdb::Vec3R eye(from.x(), from.y(), from.z());
    const openvdb::Vec3R direction = openvdb::Vec3R(to.x() - from.x(),
                                                    to.y() - from.y(),
                                                    to.z() - from.z());
    const openvdb::math::Ray<openvdb::Real> ray(eye, direction.unit());
    openvdb::Vec3R I;
    openvdb::Vec3R normal;
    openvdb::Real t;
    if (_intersector.intersectsWS(ray, I, normal, t) && t > 0.0)
    {
      // Only consider intersections within the segment bounds.
      const openvdb::Vec3R vI = I - eye;
      if (vI.lengthSqr() <= direction.lengthSqr())
      {
        p = Point_3(I.x(), I.y(), I.z());
        return true;
      }
    }
    return false;
  }

  template<typename Query>
  std::optional<Segment_3> clip_to_segment(const Query& query) const
  {
    const auto clipped = CGAL::intersection(query, bbox());
    if (clipped)
      if (const Segment_3* s = std::get_if<Segment_3>(&*clipped))
        return *s;

    return std::nullopt;
  }

  Surface_patch_index default_surface_patch_index() const
  {
    return Surface_patch_index(1, 1);
  }

public:

  struct Construct_initial_points
  {
    Construct_initial_points(const OpenVDB_mesh_domain_3& domain)
      : r_domain_(domain) {}

    template <typename OutputIteratorPoints>
    OutputIteratorPoints operator()(OutputIteratorPoints out, int n = 40)
    {
      using FT = typename GeomTraits::FT;

      const Sphere_3& sphere = r_domain_.boundingSphere();
      const Point_3 initial_center = GeomTraits().construct_center_3_object()(sphere);
      const FT squared_radius = GeomTraits().compute_squared_radius_3_object()(sphere);
      const FT radius = std::sqrt(squared_radius);

      typename CGAL::Random_points_on_sphere_3<Point_3, Point_creator> random_point_on_sphere(CGAL::to_double(radius));
      typename CGAL::Random_points_in_sphere_3<Point_3, Point_creator> random_point_in_sphere(CGAL::to_double(radius));
      typename GeomTraits::Construct_segment_3 segment_3 = GeomTraits().construct_segment_3_object();
      typename GeomTraits::Construct_vector_3 vector_3 = GeomTraits().construct_vector_3_object();
      typename GeomTraits::Construct_translated_point_3 translate = GeomTraits().construct_translated_point_3_object();

      auto construct_intersection = r_domain_.construct_intersection_object();

      Point_3 center = initial_center;
      while (n > 0)
      {
        const Point_3 p = translate(*random_point_on_sphere++, vector_3(CGAL::ORIGIN, initial_center));

        auto seg = segment_3(center, p);
        auto intersection = construct_intersection(seg);

        if (dimension(intersection) != 0)
        {
          *out++ = std::make_pair(point(intersection), index(intersection));
          --n;
        }
        else
        {
          center = translate(*random_point_in_sphere++, vector_3(CGAL::ORIGIN, initial_center));
        }
      }
      return out;
    }
  private:
    const OpenVDB_mesh_domain_3& r_domain_;
  };
  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points(*this);
  }

  struct Construct_intersection
  {
    Construct_intersection(const OpenVDB_mesh_domain_3& domain)
      : r_domain_(domain) {}

    template<typename Query>
    Intersection operator()(const Query& q)
    {
      const auto seg = r_domain_.clip_to_segment(q);
      if (seg != std::nullopt)
      {
        const Segment_3& s = seg.value();
        const auto point = intersect_clipped_segment(s.source(), s.target());
        if(point != std::nullopt)
        {
          const Surface_patch_index spi = r_domain_.default_surface_patch_index();
          return std::make_tuple(*point, index_from_surface_patch_index(spi), 2);
        }
      }
      return Intersection();
    }

  private:
    std::optional<Point_3> intersect_clipped_segment(const Point_3& from,
                                                     const Point_3& to) const
    {
      Point_3 p;
      if (r_domain_.intersects(from, to, p))
        return p;
      else
        return std::nullopt;
    }

  private:
    openvdb::FloatGrid::Ptr _grid;
    const OpenVDB_mesh_domain_3& r_domain_;
  };
  Construct_intersection construct_intersection_object() const
  {
    return Construct_intersection(*this);
  }

  struct Is_in_domain
  {
    template<typename Query>
    auto operator()(const Query&) const
    {
      return std::nullopt;
    }
  };
  Is_in_domain is_in_domain_object() const
  {
    return Is_in_domain();
  }

public:
  static const Index index_from_subdomain_index(const Subdomain_index& si)
  {
    return Index(si);
  }
  static const Index index_from_surface_patch_index(const Surface_patch_index& si)
  {
    return Index(si);
  }
  static const Subdomain_index subdomain_index(const Index& index)
  {
    if (const Subdomain_index* si = std::get_if<Subdomain_index>(&index))
      return *si;
    else
      return Subdomain_index();
  }
  static const Surface_patch_index surface_patch_index(const Index& index)
  {
    if (const Surface_patch_index* spi = std::get_if<Surface_patch_index>(&index))
      return *spi;
    else
      return Surface_patch_index();
  }

private:
  openvdb::FloatGrid::Ptr _grid;
  openvdb::tools::LevelSetRayIntersector<openvdb::FloatGrid> _intersector;
  Sphere_3 _boundingSphere;
  CGAL::Bbox_3 m_bbox;
};

}//end namespace CGAL

#endif //CGAL_MESH_3_OPENVDB_MESH_DOMAIN_3_H
