// Copyright (c) 2019-2022 Google LLC (USA).
// Copyright (c) 2025 GeometryFactory (France).
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
#ifndef CGAL_ALPHA_WRAP_2_INTERNAL_ALPHA_WRAP_AABB_GEOM_TRAITS_H
#define CGAL_ALPHA_WRAP_2_INTERNAL_ALPHA_WRAP_AABB_GEOM_TRAITS_H

#include <CGAL/license/Alpha_wrap_2.h>

#include <CGAL/Bbox_2.h>

#include <array>
#include <iostream>
#include <bitset>

namespace CGAL {
namespace Alpha_wraps_2 {
namespace internal {

template <typename K>
struct Triangle_with_outside_info
{
  using Kernel = K;
  using Triangle_2 = typename Kernel::Triangle_2;
  using Segment_2 = typename Kernel::Segment_2;

  template <typename FaceHandle>
  Triangle_with_outside_info(const FaceHandle f, const K& k)
  {
    typename K::Construct_bbox_2 bbox = k.construct_bbox_2_object();
    typename K::Construct_triangle_2 triangle = k.construct_triangle_2_object();
    typename K::Construct_segment_2 segment = k.construct_segment_2_object();

    m_tr = triangle(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());
    m_bbox = bbox(m_tr);

    for(int i=0; i<3; ++i)
    {
      if(f->neighbor(i)->is_outside())
        m_b.set(i, true);

      m_segments[i] = segment(f->vertex((i+1)%3)->point(),
                              f->vertex((i+2)%3)->point());
      m_sbox[i] = bbox(m_segments[i]);
    }
  }

  Triangle_with_outside_info& operator=(const Triangle_with_outside_info& rhs) = default;

  friend std::ostream& operator<<(std::ostream& os, const Triangle_with_outside_info& t)
  {
    os << t.m_tr;
    return os;
  }

  Triangle_2 m_tr;
  Bbox_2 m_bbox;
  std::array<Bbox_2, 3> m_sbox;
  std::array<Segment_2, 3> m_segments;
  std::bitset<3> m_b;
};

template <typename K>
class Disk_2
  : private K::Circle_2
{
  using FT = typename K::FT;
  using Point_2 = typename K::Point_2;
  using Circle_2 = typename K::Circle_2;

public:
  explicit Disk_2(const Circle_2& c) : Circle_2(c) { }

public:
  const Circle_2& boundary() const { return static_cast<const Circle_2&>(*this); }
};

template <typename GT>
class Alpha_wrap_AABB_geom_traits
  : public GT
{
public:
  using Disk_2 = internal::Disk_2<GT>;

public:
  Alpha_wrap_AABB_geom_traits(const GT& gt = GT()) : GT(gt) { }

public:
  class Construct_disk_2
  {
    using FT = typename GT::FT;
    using Point_2 = typename GT::Point_2;
    using Disk_2 = internal::Disk_2<GT>;

    const GT& m_base_traits;

public:
    Construct_disk_2(const GT& base_traits) : m_base_traits(base_traits) { }

    Disk_2 operator()(const Point_2& p, const FT sqr)
    {
      return Disk_2(m_base_traits.construct_circle_2_object()(p, sqr));
    }
  };

  // Enrich the base's Do_intersect_2 with Triangle_with_outside_info<K> and Disk_2<K> overloads
  class Do_intersect_2
    : public GT::Do_intersect_2
  {
    using Base = typename GT::Do_intersect_2;

    using Point_2 = typename GT::Point_2;
    using Segment_2 = typename GT::Segment_2;
    using Triangle_2 = typename GT::Triangle_2;
    using Disk_2 = internal::Disk_2<GT>;

    const GT& m_base_traits;

public:
    Do_intersect_2(const GT& base_traits)
      : Base(base_traits.do_intersect_2_object()),
        m_base_traits(base_traits)
    { }

    using Base::operator();

    // ======
    // Disk_2

    bool operator()(const Disk_2& disk,
                    const Point_2& p) const
    {
      return !m_base_traits.has_on_unbounded_side_2_object()(disk.boundary(), p);
    }

    bool operator()(const Disk_2& disk,
                    const Segment_2& s) const
    {
      typename GT::Construct_bbox_2 bbox = m_base_traits.construct_bbox_2_object();
      typename GT::Construct_vertex_2 vertex = m_base_traits.construct_vertex_2_object();
      typename GT::Has_on_bounded_side_2 has_on_bounded_side = m_base_traits.has_on_bounded_side_2_object();

      if(!do_overlap(bbox(disk.boundary()), bbox(s)))
         return false;

      if(Base::operator()(disk.boundary(), s))
        return true;

      return has_on_bounded_side(disk.boundary(), vertex(s, 0));
    }

    bool operator()(const Disk_2& disk,
                    const Triangle_2& tr) const
    {
      typename GT::Construct_bbox_2 bbox = m_base_traits.construct_bbox_2_object();
      typename GT::Construct_vertex_2 vertex = m_base_traits.construct_vertex_2_object();
      typename GT::Has_on_bounded_side_2 has_on_bounded_side = m_base_traits.has_on_bounded_side_2_object();

      if(!do_overlap(bbox(disk.boundary()), bbox(tr)))
         return false;

      if(Base::operator()(disk.boundary(), tr))
        return true;

      return has_on_bounded_side(disk.boundary(), vertex(tr, 0));
    }

    bool operator()(const Disk_2& disk,
                    const CGAL::Bbox_2& bb) const
    {
      typename GT::Construct_bbox_2 bbox = m_base_traits.construct_bbox_2_object();
      typename GT::Construct_point_2 point = m_base_traits.construct_point_2_object();
      typename GT::Has_on_bounded_side_2 has_on_bounded_side = m_base_traits.has_on_bounded_side_2_object();

      if(!do_overlap(bbox(disk.boundary()), bb))
         return false;

      if(Base::operator()(disk.boundary(), bb))
        return true;

      const Point_2 bbp = point(bb.xmin(), bb.ymin());

      return has_on_bounded_side(disk.boundary(), bbp);
    }

    // ======
    // Triangle_with_outside_info

    template <typename K>
    bool operator()(const Triangle_with_outside_info<K>& q,
                    const Point_2& p) const
    {
      typename GT::Construct_bbox_2 bbox = m_base_traits.construct_bbox_2_object();

      const CGAL::Bbox_2 pbox = bbox(p);
      if(!do_overlap(q.m_bbox, pbox))
        return false;

      for(int i=0; i<3; ++i)
      {
        if(!q.m_b.test(i) && do_overlap(q.m_sbox[i], pbox) && Base::operator()(q.m_segments[i], p))
          return true;
      }

      return m_base_traits.has_on_bounded_side_2_object()(q.m_tr, p);
    }

    template <typename K>
    bool operator()(const Triangle_with_outside_info<K>& q,
                    const Segment_2& s) const
    {
      typename GT::Construct_bbox_2 bbox = m_base_traits.construct_bbox_2_object();
      typename GT::Construct_vertex_2 vertex = m_base_traits.construct_vertex_2_object();

      const CGAL::Bbox_2 sbox = bbox(s);
      if(!do_overlap(q.m_bbox, sbox))
        return false;

      for(int i=0; i<3; ++i)
      {
        if(!q.m_b.test(i) && do_overlap(q.m_sbox[i], sbox) && Base::operator()(q.m_segments[i], s))
          return true;
      }

      return m_base_traits.has_on_bounded_side_2_object()(q.m_tr, vertex(s, 0));
    }

    template <typename K>
    bool operator()(const Triangle_with_outside_info<K>& q,
                    const Triangle_2& tr) const
    {
      typename GT::Construct_bbox_2 bbox = m_base_traits.construct_bbox_2_object();
      typename GT::Construct_vertex_2 vertex = m_base_traits.construct_vertex_2_object();

      const CGAL::Bbox_2 tbox = bbox(tr);
      if(!do_overlap(q.m_bbox, tbox))
        return false;

      for(int i=0; i<3; ++i)
      {
        if(!q.m_b.test(i) && do_overlap(q.m_sbox[i], tbox) && Base::operator()(q.m_segments[i], tr))
          return true;
      }

      return m_base_traits.has_on_bounded_side_2_object()(q.m_tr, vertex(tr, 0));
    }

    template <typename K>
    bool operator()(const Triangle_with_outside_info<K>& q,
                    const CGAL::Bbox_2& bb) const
    {
      if(!do_overlap(q.m_bbox, bb))
        return false;

      for(int i=0; i<3; ++i)
      {
        // this overload of do_intersect() must not filter based on q.m_b because
        // it is called from the AABB_tree's traversal with a node's bounding box,
        // and the fact that an edge is incident to an outside face is irrelevant
        // for the hierarchy of bounding boxes of the primitives.
        if(do_overlap(q.m_sbox[i], bb) && Base::operator()(q.m_segments[i], bb))
          return true;
      }

      const Point_2 bbp = m_base_traits.construct_point_2_object()(bb.xmin(), bb.ymin());
      return m_base_traits.has_on_bounded_side_2_object()(q.m_tr, bbp);
    }
  };

  Construct_disk_2 construct_disk_2_object() const
  {
    return Construct_disk_2(static_cast<const GT&>(*this));
  }

  Do_intersect_2 do_intersect_2_object() const
  {
    return Do_intersect_2(static_cast<const GT&>(*this));
  }
};

} // namespace internal
} // namespace Alpha_wraps_2
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_2_INTERNAL_ALPHA_WRAP_AABB_GEOM_TRAITS_H
