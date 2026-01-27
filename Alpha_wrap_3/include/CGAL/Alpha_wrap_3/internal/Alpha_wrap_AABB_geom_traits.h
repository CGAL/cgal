// Copyright (c) 2019-2022 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pierre Alliez
//                 Cedric Portaneri,
//                 Mael Rouxel-Labbé
//                 Andreas Fabri
//                 Michael Hemmer
//
#ifndef CGAL_ALPHA_WRAP_3_INTERNAL_ALPHA_WRAP_AABB_GEOM_TRAITS_H
#define CGAL_ALPHA_WRAP_3_INTERNAL_ALPHA_WRAP_AABB_GEOM_TRAITS_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/Bbox_3.h>

#include <array>
#include <iostream>
#include <bitset>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

template <typename K>
struct Tetrahedron_with_outside_info
{
  using Kernel = K;
  using Tetrahedron_3 = typename Kernel::Tetrahedron_3;
  using Triangle_3 = typename Kernel::Triangle_3;

  template <typename CellHandle>
  Tetrahedron_with_outside_info(const CellHandle c, const K& k)
  {
    typename K::Construct_bbox_3 bbox = k.construct_bbox_3_object();
    typename K::Construct_tetrahedron_3 tetrahedron = k.construct_tetrahedron_3_object();
    typename K::Construct_triangle_3 triangle = k.construct_triangle_3_object();

    m_tet = tetrahedron(c->vertex(0)->point(), c->vertex(1)->point(),
                        c->vertex(2)->point(), c->vertex(3)->point());
    m_bbox = bbox(m_tet);

    for(int i=0; i<4; ++i)
    {
      if(c->neighbor(i)->is_outside())
        m_b.set(i, true);

      m_triangles[i] = triangle(c->vertex((i+1)& 3)->point(),
                                c->vertex((i+2)& 3)->point(),
                                c->vertex((i+3)& 3)->point());
      m_tbox[i] = bbox(m_triangles[i]);
    }
  }

  Tetrahedron_with_outside_info& operator=(const Tetrahedron_with_outside_info& rhs) = default;

  friend std::ostream& operator<<(std::ostream& os, const Tetrahedron_with_outside_info& t)
  {
    os << t.m_tet;
    return os;
  }

  Tetrahedron_3 m_tet;
  Bbox_3 m_bbox;
  std::array<Bbox_3, 4> m_tbox;
  std::array<Triangle_3, 4> m_triangles;
  std::bitset<4> m_b;
};

template <typename K>
class Ball_3
  : private K::Sphere_3
{
  using FT = typename K::FT;
  using Point_3 = typename K::Point_3;
  using Sphere_3 = typename K::Sphere_3;

public:
  explicit Ball_3(const Sphere_3& s) : Sphere_3(s) { }

public:
  const Sphere_3& boundary() const { return static_cast<const Sphere_3&>(*this); }
};

template <typename GT>
class Alpha_wrap_AABB_geom_traits
  : public GT
{
public:
  using Ball_3 = internal::Ball_3<GT>;

public:
  Alpha_wrap_AABB_geom_traits(const GT& gt = GT()) : GT(gt) { }

public:
  class Construct_ball_3
  {
    using FT = typename GT::FT;
    using Point_3 = typename GT::Point_3;
    using Ball_3 = internal::Ball_3<GT>;

    const GT& m_base_traits;

public:
    Construct_ball_3(const GT& base_traits) : m_base_traits(base_traits) { }

    Ball_3 operator()(const Point_3& p, const FT sqr)
    {
      return Ball_3(m_base_traits.construct_sphere_3_object()(p, sqr));
    }
  };

  // Enrich the base's Do_intersect_3 with Tetrahedron_with_outside_info<K> and Ball_3<K> overloads
  class Do_intersect_3
    : public GT::Do_intersect_3
  {
    using Base = typename GT::Do_intersect_3;

    using Point_3 = typename GT::Point_3;
    using Segment_3 = typename GT::Segment_3;
    using Triangle_3 = typename GT::Triangle_3;
    using Ball_3 = internal::Ball_3<GT>;

    const GT& m_base_traits;

public:
    Do_intersect_3(const GT& base_traits)
      : Base(base_traits.do_intersect_3_object()),
        m_base_traits(base_traits)
    { }

    using Base::operator();

    // ======
    // Ball_3

    bool operator()(const Ball_3& ball,
                    const Point_3& p) const
    {
      return !m_base_traits.has_on_unbounded_side_3_object()(ball.boundary(), p);
    }

    bool operator()(const Ball_3& ball,
                    const Segment_3& s) const
    {
      typename GT::Construct_bbox_3 bbox = m_base_traits.construct_bbox_3_object();
      typename GT::Construct_vertex_3 vertex = m_base_traits.construct_vertex_3_object();
      typename GT::Has_on_bounded_side_3 has_on_bounded_side = m_base_traits.has_on_bounded_side_3_object();

      if(!do_overlap(bbox(ball.boundary()), bbox(s)))
         return false;

      if(Base::operator()(ball.boundary(), s))
        return true;

      return has_on_bounded_side(ball.boundary(), vertex(s, 0));
    }

    bool operator()(const Ball_3& ball,
                    const Triangle_3& tr) const
    {
      typename GT::Construct_bbox_3 bbox = m_base_traits.construct_bbox_3_object();
      typename GT::Construct_vertex_3 vertex = m_base_traits.construct_vertex_3_object();
      typename GT::Has_on_bounded_side_3 has_on_bounded_side = m_base_traits.has_on_bounded_side_3_object();

      if(!do_overlap(bbox(ball.boundary()), bbox(tr)))
         return false;

      if(Base::operator()(ball.boundary(), tr))
        return true;

      return has_on_bounded_side(ball.boundary(), vertex(tr, 0));
    }

    bool operator()(const Ball_3& ball,
                    const CGAL::Bbox_3& bb) const
    {
      typename GT::Construct_bbox_3 bbox = m_base_traits.construct_bbox_3_object();
      typename GT::Construct_point_3 point = m_base_traits.construct_point_3_object();
      typename GT::Has_on_bounded_side_3 has_on_bounded_side = m_base_traits.has_on_bounded_side_3_object();

      if(!do_overlap(bbox(ball.boundary()), bb))
         return false;

      if(Base::operator()(ball.boundary(), bb))
        return true;

      const Point_3 bbp = point(bb.xmin(), bb.ymin(), bb.zmin());

      return has_on_bounded_side(ball.boundary(), bbp);
    }

    // ======
    // Tetrahedron_with_outside_info

    template <typename K>
    bool operator()(const Tetrahedron_with_outside_info<K>& q,
                    const Point_3& p) const
    {
      typename GT::Construct_bbox_3 bbox = m_base_traits.construct_bbox_3_object();

      const CGAL::Bbox_3 pbox = bbox(p);
      if(!do_overlap(q.m_bbox, pbox))
        return false;

      for(int i=0; i<4; ++i)
      {
        if(!q.m_b.test(i) && do_overlap(q.m_tbox[i], pbox) && Base::operator()(q.m_triangles[i], p))
          return true;
      }

      return m_base_traits.has_on_bounded_side_3_object()(q.m_tet, p);
    }

    template <typename K>
    bool operator()(const Tetrahedron_with_outside_info<K>& q,
                    const Segment_3& s) const
    {
      typename GT::Construct_bbox_3 bbox = m_base_traits.construct_bbox_3_object();
      typename GT::Construct_vertex_3 vertex = m_base_traits.construct_vertex_3_object();

      const CGAL::Bbox_3 sbox = bbox(s);
      if(!do_overlap(q.m_bbox, sbox))
        return false;

      for(int i=0; i<4; ++i)
      {
        if(!q.m_b.test(i) && do_overlap(q.m_tbox[i], sbox) && Base::operator()(q.m_triangles[i], s))
          return true;
      }

      return m_base_traits.has_on_bounded_side_3_object()(q.m_tet, vertex(s, 0));
    }

    template <typename K>
    bool operator()(const Tetrahedron_with_outside_info<K>& q,
                    const Triangle_3& tr) const
    {
      typename GT::Construct_bbox_3 bbox = m_base_traits.construct_bbox_3_object();
      typename GT::Construct_vertex_3 vertex = m_base_traits.construct_vertex_3_object();

      const CGAL::Bbox_3 tbox = bbox(tr);
      if(!do_overlap(q.m_bbox, tbox))
        return false;

      for(int i=0; i<4; ++i)
      {
        if(!q.m_b.test(i) && do_overlap(q.m_tbox[i], tbox) && Base::operator()(q.m_triangles[i], tr))
          return true;
      }

      return m_base_traits.has_on_bounded_side_3_object()(q.m_tet, vertex(tr, 0));
    }

    template <typename K>
    bool operator()(const Tetrahedron_with_outside_info<K>& q,
                    const CGAL::Bbox_3& bb) const
    {
      if(!do_overlap(q.m_bbox, bb))
        return false;

      for(int i=0; i<4; ++i)
      {
        // this overload of do_intersect() must not filter based on q.m_b because
        // it is called from the AABB_tree's traversal with a node's bounding box,
        // and the fact that a facet is incident to an outside cell is irrelevant
        // for the hierarchy of bounding boxes of the primitives.
        if(do_overlap(q.m_tbox[i], bb) && Base::operator()(q.m_triangles[i], bb))
          return true;
      }

      const Point_3 bbp = m_base_traits.construct_point_3_object()(bb.xmin(), bb.ymin(), bb.zmin());
      return m_base_traits.has_on_bounded_side_3_object()(q.m_tet, bbp);
    }
  };

  Construct_ball_3 construct_ball_3_object() const
  {
    return Construct_ball_3(static_cast<const GT&>(*this));
  }

  Do_intersect_3 do_intersect_3_object() const
  {
    return Do_intersect_3(static_cast<const GT&>(*this));
  }
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_ALPHA_WRAP_AABB_GEOM_TRAITS_H
