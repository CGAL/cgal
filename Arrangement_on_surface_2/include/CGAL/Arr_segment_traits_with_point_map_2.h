// Copyright (c) 2026 Geometry Factory.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// author(s)     : Léo Valque

#ifndef CGAL_ARR_SEGMENT_TRAITS_WITH_POINT_MAP_2_H
#define CGAL_ARR_SEGMENT_TRAITS_WITH_POINT_MAP_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/Arr_segment_traits_2.h>
#include <boost/property_map/property_map.hpp>

namespace CGAL{

template<typename Kernel, typename PointPropertyMap, typename PointInserter>
class Arr_segment_traits_with_point_map_2
  : public Arr_segment_traits_2<Kernel>
{
public:
  using Base = Arr_segment_traits_2<Kernel>;
  using Base_point_2 = typename Kernel::Point_2;

  using Point_2 = typename boost::property_traits<PointPropertyMap>::key_type;
  // Using std::pair make an ambigous case between X_monotone_curve_2 and std::pair<Point_2, Multiplicity>
  struct X_monotone_curve_2{
    X_monotone_curve_2(){}
    X_monotone_curve_2(Point_2 a, Point_2 b): first(a), second(b){}
    X_monotone_curve_2 operator=(X_monotone_curve_2 cv){
      first=cv.first;
      second=cv.second;
      return *this;
    }
    Point_2 first;
    Point_2 second;
  };
  using Curve_2 = X_monotone_curve_2;

  using Comparison_result = typename Base::Comparison_result;
  using Multiplicity = typename Base::Multiplicity;

  // CGAL_static_assertion(std::is_same_v<Base_point_2, typename boost::property_traits<PointPropertyMap>::value_type>);
private:
  mutable PointPropertyMap point_map;
  PointInserter    insert; // A function that given a PointPropertyMap and a Base_point_2, insert the point in the property map and returns the corresponding Point_2
  Base             base;

public:
  Arr_segment_traits_with_point_map_2(PointPropertyMap pm, PointInserter inserter) : point_map(pm), insert(inserter), base(){}

private:
  Base_point_2 to_base(const Point_2& v) const{
    return get(point_map, v);
  }

  typename Base::X_monotone_curve_2
  to_base(const X_monotone_curve_2& cv) const{
    return typename Base::X_monotone_curve_2(to_base(cv.first), to_base(cv.second));
  }

  Point_2 from_base(const Base_point_2& p) const{
    return insert(point_map, p);
  }

public:
  class Compare_xy_2 {
    const Arr_segment_traits_with_point_map_2* traits;
  public:
    Compare_xy_2(const Arr_segment_traits_with_point_map_2* traits_) : traits(traits_) {}

    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const{
      if(p1==p2) return EQUAL;
      return traits->base.compare_xy_2_object()(traits->to_base(p1), traits->to_base(p2));
    }
  };
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(this); }

  class Construct_min_vertex_2 {
    const Arr_segment_traits_with_point_map_2* traits;
  public:
    Construct_min_vertex_2(const Arr_segment_traits_with_point_map_2* traits_) : traits(traits_) {}

    Point_2 operator()(const X_monotone_curve_2& cv) const{
      auto cmp = traits->compare_xy_2_object();
      return (cmp(cv.first, cv.second) == SMALLER) ? cv.first : cv.second;
    }
  };
  Construct_min_vertex_2 construct_min_vertex_2_object() const { return Construct_min_vertex_2(this); }

  class Construct_max_vertex_2 {
    const Arr_segment_traits_with_point_map_2* traits;
  public:
    Construct_max_vertex_2(const Arr_segment_traits_with_point_map_2* traits_) : traits(traits_) {}

    Point_2 operator()(const X_monotone_curve_2& cv) const{
      auto cmp = traits->compare_xy_2_object();
      return (cmp(cv.first, cv.second) == LARGER) ? cv.first : cv.second;
    }
  };
  Construct_max_vertex_2 construct_max_vertex_2_object() const { return Construct_max_vertex_2(this); }

  class Is_vertical_2 {
    const Arr_segment_traits_with_point_map_2* traits;
  public:
    Is_vertical_2(const Arr_segment_traits_with_point_map_2* traits_) : traits(traits_) {}

    bool operator()(const X_monotone_curve_2& cv) const{
      return traits->base.is_vertical_2_object()(traits->to_base(cv));
    }
  };
  Is_vertical_2 is_vertical_2_object() const { return Is_vertical_2(this); }

  class Intersect_2 {
    const Arr_segment_traits_with_point_map_2* traits;
  public:
    Intersect_2(const Arr_segment_traits_with_point_map_2* traits_) : traits(traits_) {}

    template<typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              OutputIterator oi) const
    {
      using Intersection_point = std::pair<Point_2, Multiplicity>;

      auto b1 = traits->to_base(cv1);
      auto b2 = traits->to_base(cv2);

      if (! traits->base.do_intersect_2_object()(b1, b2)) return oi;

      // An intersection is guaranteed.

      // Intersect the two supporting lines.
      auto res = Kernel().intersect_2_object()(b1.line(), b2.line());
      CGAL_assertion(bool(res));

      // Check if we have a single intersection point.
      const Base_point_2* ip = std::get_if<Base_point_2>(&*res);
      if (ip != nullptr) {
        Intersection_point ip_mult;
        if(*ip == b1.source())
          ip_mult = Intersection_point(cv1.first, 1);
        else if(*ip == b1.target())
          ip_mult = Intersection_point(cv1.second, 1);
        else if(*ip == b2.source())
          ip_mult = Intersection_point(cv2.first, 1);
        else if(*ip == b2.target())
          ip_mult = Intersection_point(cv2.second, 1);
        else {
          auto new_p = traits->from_base(*ip);
          ip_mult = Intersection_point(new_p, 1);
        }
        *oi++ = ip_mult;
        return oi;
      }

      // In this case, the two supporting lines overlap.
      // The overlapping segment is therefore [p_l,p_r], where p_l is the
      // rightmost of the two left endpoints and p_r is the leftmost of the
      // two right endpoints.
      auto compare_xy = traits->base.compare_xy_2_object();
      auto left = traits->construct_min_vertex_2_object();
      auto right = traits->construct_max_vertex_2_object();
      Point_2 p_l = (compare_xy(b1.left(), b2.left()) == SMALLER) ?
        left(cv2) : left(cv1);
      Point_2 p_r = (compare_xy(b1.right(), b2.right()) == SMALLER) ?
        right(cv1) : right(cv2);

      // Examine the resulting segment.
      if (p_l == p_r) {
        // The two segment have the same supporting line, but they just share
        // a common endpoint. Thus we have an intersection point, but we leave
        // the multiplicity of this point undefined.
        Intersection_point ip_mult(p_r, 0);
        *oi++ = ip_mult;
        return oi;
      }
      X_monotone_curve_2 overlap_seg(p_l, p_r);
      *oi++ = overlap_seg;
      return oi;
    }
  };
  Intersect_2 intersect_2_object() const { return Intersect_2(this); }

  class Split_2 {
    const Arr_segment_traits_with_point_map_2* traits;
  public:
    Split_2(const Arr_segment_traits_with_point_map_2* traits_) : traits(traits_){}

    void operator()(const X_monotone_curve_2& cv,
                    const Point_2& p,
                    X_monotone_curve_2& cv1,
                    X_monotone_curve_2& cv2) const
    {
      cv1.first = traits->construct_min_vertex_2_object()(cv);
      cv1.second = p;
      cv2.first = p;
      cv2.second = traits->construct_max_vertex_2_object()(cv);
    }
  };
  Split_2 split_2_object() const { return Split_2(this); }

  class Make_x_monotone_2 {
  public:
    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const {
      // Wrap the segment with a variant.
      *oi++ = cv;
      return oi;
    }
  };
  Make_x_monotone_2 make_x_monotone_2_object() const { return Make_x_monotone_2(); }

  class Compare_y_at_x_2 {
    const Arr_segment_traits_with_point_map_2* traits;
  public:
    Compare_y_at_x_2(const Arr_segment_traits_with_point_map_2* traits_) : traits(traits_){}

    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& cv) const {
      return traits->base.compare_y_at_x_2_object()(traits->to_base(p), traits->to_base(cv));
    }
  };
  Compare_y_at_x_2 compare_y_at_x_2_object() const { return Compare_y_at_x_2(this); }

  class Compare_y_at_x_right_2 {
    const Arr_segment_traits_with_point_map_2* traits;
  public:
    Compare_y_at_x_right_2(const Arr_segment_traits_with_point_map_2* traits_) : traits(traits_){}

    Comparison_result operator()(const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2,
                                 const Point_2& p) const {
      return traits->base.compare_y_at_x_right_2_object()(traits->to_base(cv1), traits->to_base(cv2), traits->to_base(p));
    }
  };
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const { return Compare_y_at_x_right_2(this); }

  class Compare_y_at_x_left_2 {
    const Arr_segment_traits_with_point_map_2* traits;
  public:
    Compare_y_at_x_left_2(const Arr_segment_traits_with_point_map_2* traits_) : traits(traits_){}

    Comparison_result operator()(const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2,
                                 const Point_2& p) const {
      return traits->base.compare_y_at_x_left_2_object()(traits->to_base(cv1), traits->to_base(cv2), traits->to_base(p));
    }
  };
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const { return Compare_y_at_x_left_2(this); }

  class Equal_2 {
    const Arr_segment_traits_with_point_map_2* traits;
  public:
    Equal_2(const Arr_segment_traits_with_point_map_2* traits_) : traits(traits_){}

    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const {
      return ((*this)(cv1.first, cv2.first) && (*this)(cv1.second, cv2.second)) ||
             ((*this)(cv1.first, cv2.second) && (*this)(cv1.second, cv2.first));
    }

    bool operator()(const Point_2& p1, const Point_2& p2) const {
      CGAL_assertion_code(auto equal = traits->base.equal_2_object();)
      CGAL_assertion( (p1==p2) == equal(traits->to_base(p1), traits->to_base(p2)));
      return (p1==p2);
    }
  };
  Equal_2 equal_2_object() const { return Equal_2(this); }
};

template<typename PointPropertyMap, typename PointInserter>
Arr_segment_traits_with_point_map_2<typename Kernel_traits<typename PointPropertyMap::value_type>::Kernel, PointPropertyMap, PointInserter>
make_arr_segments_traits_with_point_map(PointPropertyMap pm, PointInserter inserter){
  using Kernel = typename Kernel_traits<typename PointPropertyMap::value_type>::Kernel;
  return Arr_segment_traits_with_point_map_2<Kernel,PointPropertyMap, PointInserter>(pm, inserter);
}

} // end of namespace CGAL

#endif // CGAL_ARR_SEGMENT_TRAITS_WITH_POINT_MAP_2_H