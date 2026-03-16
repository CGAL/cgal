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
#include <boost/container/small_vector.hpp>

namespace CGAL::internal {

template<typename FloatSnapRoundingTraits, typename PointPropertyMap>
class Wrap_float_snap_rounding_traits_2
{
public:
  using Base = FloatSnapRoundingTraits;

  using Base_point_2   = typename Base::Point_2;
  using Base_segment_2 = typename Base::Segment_2;
  using FT             = typename Base::FT;

  using Evaluation_tag = typename Base::Evaluation_tag;

  using Multiplicity = typename Base::Multiplicity;
  using Has_left_category = typename Base::Has_left_category;
  using Has_merge_category = typename Base::Has_merge_category;

  using Left_side_category = typename Base::Left_side_category;
  using Bottom_side_category = typename Base::Bottom_side_category;
  using Top_side_category = typename Base::Top_side_category;
  using Right_side_category = typename Base::Right_side_category;

  using Point_2 = typename boost::property_traits<PointPropertyMap>::key_type;
  // using Segment_2 = std::pair<Point_2, Point_2>;
  struct Segment_2 {
    Segment_2(){}
    Segment_2(Point_2 a, Point_2 b): src(a), trg(b){}
    Segment_2(const Segment_2 &s) : src(s.src), trg(s.trg), polyline_indices(s.polyline_indices){}
    Segment_2& operator=(const Segment_2 &s){
      src = s.src;
      trg = s.trg;
      polyline_indices = s.polyline_indices;
      return *this;
    }
    bool operator<(const Segment_2 &s) const {
      return (src == s.src) ? trg < s.trg: src < s.src;
    }
    void add_polyline(std::size_t pl_idx) const {
      polyline_indices.emplace_back(pl_idx);
    }

    Point_2 src;
    Point_2 trg;
    mutable boost::container::small_vector< std::size_t, 1> polyline_indices; // Indices of all overlapping polylines on the segment
                                                                              // Made mutable to be fullfill by a set
  };
  using X_monotone_curve_2 = Segment_2;
  using Curve_2 = X_monotone_curve_2;

private:
  const Base &base;
  mutable PointPropertyMap point_map;

public:
  Wrap_float_snap_rounding_traits_2(const FloatSnapRoundingTraits &traits, PointPropertyMap pm): base(traits), point_map(pm){}

private:
  Base_point_2 to_base(const Point_2& p) const {
    return get(point_map, p);
  }

  Base_segment_2 to_base(const Segment_2& s) const {
    return Base_segment_2(to_base(s.src), to_base(s.trg));
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
      return (cmp(cv.src, cv.trg) == SMALLER) ? cv.src : cv.trg;
    }
  };
  Construct_min_vertex_2 construct_min_vertex_2_object() const { return Construct_min_vertex_2(this); }

  class Construct_max_vertex_2 {
    const Wrap_float_snap_rounding_traits_2* traits;
  public:
    Construct_max_vertex_2(const Wrap_float_snap_rounding_traits_2* traits_) : traits(traits_) {}

    Point_2 operator()(const X_monotone_curve_2& cv) const{
      auto cmp = traits->compare_xy_2_object();
      return (cmp(cv.src, cv.trg) == LARGER) ? cv.src : cv.trg;
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
                                 const Point_2& p) const {
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
                                 const Point_2& p) const {
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
      return ((*this)(cv1.src, cv2.src) && (*this)(cv1.trg, cv2.trg)) ||
             ((*this)(cv1.src, cv2.trg) && (*this)(cv1.trg, cv2.src));
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
      cv1.src = traits->construct_min_vertex_2_object()(cv);
      cv1.trg= p;
      cv2.src = p;
      cv2.trg = traits->construct_max_vertex_2_object()(cv);
    }
  };
  Split_2 split_2_object() const { return Split_2(this); }

  using Evaluate                          = typename Base::Evaluate;
  using Converter_to_exact                = typename Base::Converter_to_exact ;
  using Converter_from_exact              = typename Base::Converter_from_exact;
  using Construct_point_at_x_on_segment_2 = typename Base::Construct_point_at_x_on_segment_2;
  using Compute_squared_round_bound_2     = typename Base::Compute_squared_round_bound_2;
  using Construct_rounded_point_2         = typename Base::Construct_rounded_point_2;
  using Compare_squared_distance_2        = typename Base::Compare_squared_distance_2;

  Evaluate evaluate_object() const{ return base.evaluate_object(); }
  Converter_to_exact converter_to_exact_object() const{ return base.converter_to_exact_object(); }
  Converter_from_exact converter_from_exact_object() const{ return base.converter_from_exact_object(); }
  Construct_point_at_x_on_segment_2 construct_point_at_x_on_segment_2_object() const{ return base.construct_point_at_x_on_segment_2_object(); }
  Compute_squared_round_bound_2 compute_squared_round_bound_2_object() const{ return base.compute_squared_round_bound_2_object(); }
  Construct_rounded_point_2 construct_rounded_point_2_object() const{ return base.construct_rounded_point_2_object(); }
  Compare_squared_distance_2 compare_squared_distance_2_object() const{ return base.compare_squared_distance_2_object(); }
};

template<typename FloatSnapRoundingTraits, typename PointPropertyMap>
Wrap_float_snap_rounding_traits_2<FloatSnapRoundingTraits, PointPropertyMap>
make_wrap_float_snap_rounding_traits_2(FloatSnapRoundingTraits& traits, PointPropertyMap pm){
  return Wrap_float_snap_rounding_traits_2<FloatSnapRoundingTraits, PointPropertyMap>(traits, pm);
}

} // namespace CGAL::internal

#endif