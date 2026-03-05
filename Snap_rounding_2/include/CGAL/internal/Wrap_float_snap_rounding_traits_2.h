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

#include <CGAL/license/Snap_rounding_2.h>

#ifndef CGAL_WRAP_FLOAT_SNAP_ROUNDING_TRAITS_2_H
#define CGAL_WRAP_FLOAT_SNAP_ROUNDING_TRAITS_2_H

#include <CGAL/license/Snap_rounding_2.h>

#include <boost/property_map/property_map.hpp>

namespace CGAL::internal {

template<typename FloatSnapRoundingTraits, typename PointPropertyMap, typename PointInserter = void*>
class Wrap_float_snap_rounding_traits_2 : public FloatSnapRoundingTraits
{
public:
  using Base = FloatSnapRoundingTraits;

  using Base_point_2   = typename Base::Point_2;
  using Base_segment_2 = typename Base::Segment_2;
  using FT             = typename Base::FT;

  using Point_2 = typename boost::property_traits<PointPropertyMap>::key_type;
  using Segment_2 = std::pair<Point_2, Point_2>;
  using X_monotone_curve_2 = Segment_2;
  using Curve_2 = X_monotone_curve_2;

private:
  Base base;
  mutable PointPropertyMap point_map;
  PointInserter insert;

public:
  Wrap_float_snap_rounding_traits_2(FloatSnapRoundingTraits &traits, PointPropertyMap pm, PointInserter inserter): base(traits), point_map(pm), insert(inserter){}

private:
  Base_point_2 to_base(const Point_2& p) const {
    return get(point_map, p);
  }

  Base_segment_2 to_base(const Segment_2& s) const {
    return Base_segment_2(to_base(s.first), to_base(s.second));
  }

  Point_2 from_base(const Base_point_2& p) const {
    return put(point_map, p);
  }

public:
  class Compare_xy_2 {
    const Wrap_float_snap_rounding_traits_2* traits;
  public:
    Compare_xy_2(const Wrap_float_snap_rounding_traits_2* traits_) : traits(traits_) {}

    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const{
      if(p1==p2) return EQUAL;
      return traits->base.compare_xy_2_object()(traits->to_base(p1), traits->to_base(p2));
    }
  };
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(this); }

  class Construct_min_vertex_2 {
    const Wrap_float_snap_rounding_traits_2* traits;
  public:
    Construct_min_vertex_2(const Wrap_float_snap_rounding_traits_2* traits_) : traits(traits_) {}

    Point_2 operator()(const X_monotone_curve_2& cv) const{
      auto cmp = traits->compare_xy_2_object();
      return (cmp(cv.first, cv.second) == SMALLER) ? cv.first : cv.second;
    }
  };
  Construct_min_vertex_2 construct_min_vertex_2_object() const { return Construct_min_vertex_2(this); }

  class Construct_max_vertex_2 {
    const Wrap_float_snap_rounding_traits_2* traits;
  public:
    Construct_max_vertex_2(const Wrap_float_snap_rounding_traits_2* traits_) : traits(traits_) {}

    Point_2 operator()(const X_monotone_curve_2& cv) const{
      auto cmp = traits->compare_xy_2_object();
      return (cmp(cv.first, cv.second) == LARGER) ? cv.first : cv.second;
    }
  };
  Construct_max_vertex_2 construct_max_vertex_2_object() const { return Construct_max_vertex_2(this); }

  class Is_vertical_2 {
    const Wrap_float_snap_rounding_traits_2* traits;
  public:
    Is_vertical_2(const Wrap_float_snap_rounding_traits_2* traits_) : traits(traits_) {}

    bool operator()(const X_monotone_curve_2& cv) const{
      return traits->base.is_vertical_2_object()(traits->to_base(cv));
    }
  };
  Is_vertical_2 is_vertical_2_object() const { return Is_vertical_2(this); }

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
    const Wrap_float_snap_rounding_traits_2* traits;
  public:
    Compare_y_at_x_2(const Wrap_float_snap_rounding_traits_2* traits_) : traits(traits_){}

    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& cv) const {
      return traits->base.compare_y_at_x_2_object()(traits->to_base(p), traits->to_base(cv));
    }
  };
  Compare_y_at_x_2 compare_y_at_x_2_object() const { return Compare_y_at_x_2(this); }

  class Compare_y_at_x_right_2 {
    const Wrap_float_snap_rounding_traits_2* traits;
  public:
    Compare_y_at_x_right_2(const Wrap_float_snap_rounding_traits_2* traits_) : traits(traits_){}

    Comparison_result operator()(const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2,
                                 const Point_2& CGAL_assertion_code(p)) const {
      return traits->base.compare_y_at_x_right_2_object()(traits->to_base(cv1), traits->to_base(cv2), traits->to_base(p));
    }
  };
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const { return Compare_y_at_x_right_2(this); }

  class Compare_y_at_x_left_2 {
    const Wrap_float_snap_rounding_traits_2* traits;
  public:
    Compare_y_at_x_left_2(const Wrap_float_snap_rounding_traits_2* traits_) : traits(traits_){}

    Comparison_result operator()(const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2,
                                 const Point_2& CGAL_assertion_code(p)) const {
      return traits->base.compare_y_at_x_left_2_object()(traits->to_base(cv1), traits->to_base(cv2), traits->to_base(p));
    }
  };
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const { return Compare_y_at_x_left_2(this); }

  class Equal_2 {
    const Wrap_float_snap_rounding_traits_2* traits;
  public:
    Equal_2(const Wrap_float_snap_rounding_traits_2* traits_) : traits(traits_){}

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

  class Split_2 {
    const Wrap_float_snap_rounding_traits_2* traits;
  public:
    Split_2(const Wrap_float_snap_rounding_traits_2* traits_) : traits(traits_){}

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
};

template<typename FloatSnapRoundingTraits, typename PointPropertyMap, typename PointInserter>
Wrap_float_snap_rounding_traits_2<FloatSnapRoundingTraits, PointPropertyMap, PointInserter>
make_wrap_float_snap_rounding_traits_2(FloatSnapRoundingTraits& traits, PointPropertyMap pm, PointInserter inserter){
  return Wrap_float_snap_rounding_traits_2<FloatSnapRoundingTraits, PointPropertyMap, PointInserter>(traits, pm, inserter);
}

} // namespace CGAL::internal

#endif