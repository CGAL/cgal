// Copyright (c) 2015  Tel-Aviv University (Israel).
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
// Author(s): Sebastian Morr    <sebastian@morr.cc>

#ifndef CGAL_AABB_TRAITS_2_H
#define CGAL_AABB_TRAITS_2_H

namespace CGAL {

template<typename GeomTraits, typename AABB_primitive_>
class AABB_traits_2
{

public:

  typedef AABB_primitive_ Primitive;
  typedef typename Primitive::Id Id;
  typedef typename Primitive::Datum Datum;
  typedef typename Primitive::Container Container;

  typedef typename GeomTraits::Point_2 Point;
  typedef typename GeomTraits::Vector_2 Vector_2;
  typedef typename CGAL::Bbox_2 Bounding_box;

  typedef typename std::pair<Object, Id> Object_and_primitive_id;
  typedef typename std::pair<Point, Id> Point_and_primitive_id;

  // Types for AABB_tree
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_2 Point_3;
  typedef typename GeomTraits::Circle_2 Sphere_3;
  typedef typename GeomTraits::Iso_rectangle_2 Iso_cuboid_3;
  typedef typename GeomTraits::Construct_center_2 Construct_center_3;
  typedef typename GeomTraits::Construct_iso_rectangle_2 Construct_iso_cuboid_3;
  typedef typename GeomTraits::Construct_min_vertex_2 Construct_min_vertex_3;
  typedef typename GeomTraits::Construct_max_vertex_2 Construct_max_vertex_3;
  typedef typename GeomTraits::Compute_squared_radius_2 Compute_squared_radius_3;
  typedef typename GeomTraits::Cartesian_const_iterator_2
  Cartesian_const_iterator_3;
  typedef typename GeomTraits::Construct_cartesian_const_iterator_2
  Construct_cartesian_const_iterator_3;

  AABB_traits_2(const Point &translation_point): m_translation_point(
      translation_point)
  {
    m_interval_x = Interval_nt<true>(to_interval(translation_point.x()));
    m_interval_y = Interval_nt<true>(to_interval(translation_point.y()));
  };

  AABB_traits_2()
  {
    m_translation_point = Point(0, 0);
    m_interval_x = Interval_nt<true>(0);
    m_interval_y = Interval_nt<true>(0);
  };

  Interval_nt<true> get_interval_x() const
  {
    return m_interval_x;
  }

  Interval_nt<true> get_interval_y() const
  {
    return m_interval_y;
  }

  Point get_translation_point() const
  {
    return m_translation_point;
  }

  // Put the n/2 smallest primitives in the front, the n/2 largest primitives
  // in the back. They are compared along the bbox' longest axis.
  class Sort_primitives
  {
  public:
    template<typename PrimitiveIterator>
    void operator()(PrimitiveIterator first,
                    PrimitiveIterator beyond,
                    const Bounding_box &bbox) const
    {
      PrimitiveIterator middle = first + (beyond - first) / 2;

      if (bbox.xmax()-bbox.xmin() >= bbox.ymax()-bbox.ymin())
      {
        std::nth_element(first, middle, beyond, less_x); // sort along x
      }
      else
      {
        std::nth_element(first, middle, beyond, less_y); // sort along y
      }
    }
  };

  Sort_primitives sort_primitives_object() const
  {
    return Sort_primitives();
  }

  // Computes the bounding box of a set of primitives
  class Compute_bbox
  {
  public:
    template<typename ConstPrimitiveIterator>
    Bounding_box operator()(ConstPrimitiveIterator first,
                            ConstPrimitiveIterator beyond) const
    {
      Bounding_box bbox = first->datum().bbox();

      for (++first; first != beyond; ++first)
      {
        bbox = bbox + first->datum().bbox();
      }

      return bbox;
    }
  };

  Compute_bbox compute_bbox_object() const
  {
    return Compute_bbox();
  }

  class Do_intersect
  {

  private:

    AABB_traits_2 *m_traits;

  public:

    Do_intersect(AABB_traits_2 *_traits): m_traits(_traits) {}

    bool operator()(const Bounding_box &q, const Bounding_box &bbox) const
    {
      Interval_nt<true> x1 = Interval_nt<true>(q.xmin(), q.xmax());
      Interval_nt<true> y1 = Interval_nt<true>(q.ymin(), q.ymax());
      Interval_nt<true> x2 = Interval_nt<true>(bbox.xmin(),
                             bbox.xmax()) + m_traits->get_interval_x();
      Interval_nt<true> y2 = Interval_nt<true>(bbox.ymin(),
                             bbox.ymax()) + m_traits->get_interval_y();

      return x1.do_overlap(x2) && y1.do_overlap(y2);
    }

    bool operator()(const Primitive &q, const Bounding_box &bbox) const
    {
      Interval_nt<true> x1 = Interval_nt<true>(q.datum().bbox().xmin(),
                             q.datum().bbox().xmax());
      Interval_nt<true> y1 = Interval_nt<true>(q.datum().bbox().ymin(),
                             q.datum().bbox().ymax());
      Interval_nt<true> x2 = Interval_nt<true>(bbox.xmin(),
                             bbox.xmax()) + m_traits->get_interval_x();
      Interval_nt<true> y2 = Interval_nt<true>(bbox.ymin(),
                             bbox.ymax()) + m_traits->get_interval_y();

      return x1.do_overlap(x2) && y1.do_overlap(y2);
    }

    bool operator()(const Bounding_box &q, const Primitive &pr) const
    {
      Datum tr_pr = pr.datum().transform(typename GeomTraits::Aff_transformation_2(
                                         Translation(),
                                         Vector_2(ORIGIN, m_traits->get_translation_point())));

      return do_overlap(q, tr_pr.bbox());
    }

    bool operator()(const Primitive &q, const Primitive &pr) const
    {
      Datum tr_pr = pr.datum().transform(typename GeomTraits::Aff_transformation_2(
                                         Translation(), Vector_2(ORIGIN, m_traits->get_translation_point())));

      if (!do_overlap(q.datum().bbox(), tr_pr.bbox()))
      {
        return false;
      }

      return do_intersect(q.datum(), tr_pr);
    }
  };

  Do_intersect do_intersect_object()
  {
    return Do_intersect(this);
  }

private:

  Point m_translation_point;
  Interval_nt<true> m_interval_x;
  Interval_nt<true> m_interval_y;

  // Comparison functions
  static bool less_x(const Primitive &pr1, const Primitive &pr2)
  {
    return pr1.reference_point().x() < pr2.reference_point().x();
  }

  static bool less_y(const Primitive &pr1, const Primitive &pr2)
  {
    return pr1.reference_point().y() < pr2.reference_point().y();
  }
};

} // namespace CGAL

#endif
