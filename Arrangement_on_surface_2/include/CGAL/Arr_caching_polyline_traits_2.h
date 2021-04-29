#ifndef CGAL_ARR_CACHING_POLYLINE_TRAITS_2_H
#define CGAL_ARR_CACHING_POLYLINE_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

#include <fstream>

#include <boost/variant.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/tags.h>
#include <CGAL/intersections.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_polycurve_basic_traits_2.h>
#include <CGAL/Arr_caching_polyline_subtraits_2.h>

namespace CGAL {

template <typename Kernel_, typename Range_>
class Arr_caching_polyline_traits_2
  : public Arr_polycurve_basic_traits_2
    <Arr_caching_polyline_subtraits_2<Kernel_, Range_>,
     internal::X_monotone_caching_polyline_2<Kernel_, Range_> >
{
public:

  using Kernel = Kernel_;
  using Range = Range_;
  using Subcurve_traits_2 = Arr_caching_polyline_subtraits_2<Kernel, Range>;
  using Curve_2 = internal::Caching_polyline_2<Kernel, Range>;
  using X_monotone_curve_2 = internal::X_monotone_caching_polyline_2<Kernel, Range>;

private:
  using Base = Arr_polycurve_basic_traits_2<Subcurve_traits_2, X_monotone_curve_2>;
  using Self = Arr_caching_polyline_traits_2<Kernel, Range>;
  using Extreme_point = typename Curve_2::Extreme_point;

public:

  using Has_left_category = typename Base::Has_left_category;
  using Has_do_intersect_category = typename Base::Has_do_intersect_category;

  using Left_side_category = typename Base::Left_side_category;
  using Bottom_side_category = typename Base::Bottom_side_category;
  using Top_side_category = typename Base::Top_side_category;
  using Right_side_category = typename Base::Right_side_category;

  using Are_all_sides_oblivious_tag = typename Base::Are_all_sides_oblivious_tag;

  using X_monotone_subcurve_2 = typename Base::X_monotone_subcurve_2;
  using Size = typename Base::Size;
  using size_type = typename Base::size_type;

  using Point_2 = typename Base::Point_2;

  using Compare_x_2 = typename Base::Compare_x_2;
  using Compare_xy_2 = typename Base::Compare_xy_2;
  using Construct_min_vertex_2 = typename Base::Construct_min_vertex_2;
  using Construct_max_vertex_2 = typename Base::Construct_max_vertex_2;
  using Is_vertical_2 = typename Base::Is_vertical_2;
  using Compare_y_at_x_2 = typename Base::Compare_y_at_x_2;
  using Compare_y_at_x_left_2 = typename Base::Compare_y_at_x_left_2;
  using Compare_y_at_x_right_2 = typename Base::Compare_y_at_x_right_2;
  using Equal_2 = typename Base::Equal_2;
  using Compare_endpoints_xy_2 = typename Base::Compare_endpoints_xy_2;
  // using Construct_opposite_2 = typename Base::Construct_opposite_2;
  using Approximate_2 = typename Base::Approximate_2;
  using Construct_x_monotone_curve_2 = typename Base::Construct_x_monotone_curve_2;
  using Parameter_space_in_x_2 = typename Base::Parameter_space_in_x_2;
  using Parameter_space_in_y_2 = typename Base::Parameter_space_in_y_2;
  using Compare_x_on_boundary_2 = typename Base::Compare_x_on_boundary_2;
  using Compare_x_at_limit_2 = typename Base::Compare_x_at_limit_2;
  using Compare_x_near_boundary_2 = typename Base::Compare_x_near_boundary_2;
  using Compare_x_near_limit_2 = typename Base::Compare_x_near_limit_2;
  using Compare_y_on_boundary_2 = typename Base::Compare_y_on_boundary_2;
  using Compare_y_near_boundary_2 = typename Base::Compare_y_near_boundary_2;
  using Is_on_y_identification_2 = typename Base::Is_on_y_identification_2;
  using Is_on_x_identification_2 = typename Base::Is_on_x_identification_2;
  using Trim_2 = typename Base::Trim_2;

  using Has_merge_category = typename Subcurve_traits_2::Has_merge_category;
  using Multiplicity = typename Subcurve_traits_2::Multiplicity;
  using Subcurve_2 = typename Subcurve_traits_2::Curve_2;


  Arr_caching_polyline_traits_2() : Base() { }
  Arr_caching_polyline_traits_2(const Subcurve_traits_2* geom_traits) : Base(geom_traits) { }

  class Construct_opposite_2 {
  public:
    X_monotone_curve_2 operator()(const X_monotone_curve_2& xcv) const
    {
      return xcv.opposite();
    }
  };
  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(); }

  class Make_x_monotone_2
  {
  protected:
    using Traits = Arr_caching_polyline_traits_2<Kernel, Range>;
    const Traits& m_traits;
    Make_x_monotone_2(const Traits& traits) : m_traits(traits) {}
    friend class Arr_caching_polyline_traits_2<Kernel, Range>;
  public:
    using Curve_iterator = typename Curve_2::iterator;

    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const
    {
      using Make_x_monotone_result = boost::variant<Point_2, X_monotone_curve_2>;

      const Kernel& kernel = *m_traits.subcurve_traits_2();

      auto compare_x_2 = m_traits.compare_x_2_object();

      Curve_iterator start = cv.points_begin();
      Curve_iterator it_prev = start + 1;
      Curve_iterator it_curr = it_prev + 1;

      Comparison_result previous_comp = compare_x_2(*start, *it_prev);
      for (; it_curr != cv.points_end(); ++ it_prev, ++ it_curr)
      {
        Comparison_result current_comp = compare_x_2(*it_prev, *it_curr);
        if (current_comp != previous_comp)
        {
          previous_comp = current_comp;
          *oi ++ = Make_x_monotone_result(X_monotone_curve_2(kernel, start, it_curr));
          start = it_prev;
        }
      }
      *oi ++ = Make_x_monotone_result(X_monotone_curve_2(kernel, start, cv.points_end()));

      return oi;
    }
  };
  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(*this); }

  class Split_2
  {
  protected:
    using Traits = Arr_caching_polyline_traits_2<Kernel, Range>;
    const Traits& m_traits;
    Split_2(const Traits& traits) : m_traits(traits) {}
    friend class Arr_caching_polyline_traits_2<Kernel, Range>;
  public:

    void operator()(const X_monotone_curve_2& xcv, const Point_2& p,
                    X_monotone_curve_2& xcv1, X_monotone_curve_2& xcv2) const
    {
      const Subcurve_traits_2* geom_traits = m_traits.subcurve_traits_2();
      const Kernel& kernel = *geom_traits;
      auto min_vertex = geom_traits->construct_min_vertex_2_object();
      auto max_vertex = geom_traits->construct_max_vertex_2_object();
      auto equal = geom_traits->equal_2_object();
      auto cmp_seg_endpts = geom_traits->compare_endpoints_xy_2_object();

      // Make sure the split point is not one of the curve endpoints.
      CGAL_precondition((! equal(m_traits.
                                 construct_min_vertex_2_object()(xcv), p)));
      CGAL_precondition((! equal(m_traits.
                                 construct_max_vertex_2_object()(xcv), p)));

      CGAL_precondition_msg(xcv.number_of_subcurves() > 0,
                            "Cannot split a polycurve of length zero.");

      Comparison_result dir = cmp_seg_endpts(xcv[0]);

      // Locate the subcurve on the polycurve xcv that contains p.
      std::size_t i = m_traits.locate(xcv, p);

      CGAL_precondition(i != Traits::INVALID_INDEX);

      /*
          i
      A B C D E F

      1 = A B C D (0 -> i+2)
      2 = D E F (i+1 -> end)

      1 = A B C (0 -> i+1)
      2 = C D E F (i-> end)
      */

      if (equal(max_vertex(xcv[i]), p)) {
        // The entire i'th subcurve belongs to xcv1:
        xcv1 = X_monotone_curve_2(kernel, xcv.points_begin(), xcv[i+1]);
        xcv2 = X_monotone_curve_2(kernel, xcv[i], xcv.points_end());
      }
      else if (equal(min_vertex(xcv[i]), p)) {
        // The entire i'th subcurves belongs to xcv2:
        xcv1 = X_monotone_curve_2(kernel, xcv.points_begin(), xcv[i]);
        xcv2 = X_monotone_curve_2(kernel, xcv[i-1], xcv.points_end());
      }
      else {
        // The i'th subcurve should be split: The left part(seg1)
        // goes to xcv1, and the right part(seg2) goes to xcv2.
        auto p_ptr = xcv.extreme_point(p, i);
        xcv1 = X_monotone_curve_2(kernel, Extreme_point(), xcv.points_begin(), xcv[i+1], p_ptr);
        xcv2 = X_monotone_curve_2(kernel, p_ptr, xcv[i+1], xcv.points_end(), Extreme_point());
      }

      if (dir != SMALLER) std::swap(xcv1, xcv2);
    }

  };
  Split_2 split_2_object() const
  { return Split_2(*this); }

  class Intersect_2
  {
  protected:
    using Traits = Arr_caching_polyline_traits_2<Kernel, Range>;
    const Traits& m_traits;
    Intersect_2(const Traits& traits) : m_traits(traits) {}
    friend class Arr_caching_polyline_traits_2<Kernel, Range>;
  public:
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              OutputIterator oi) const
    {
      typedef std::pair<Point_2, Multiplicity>        Intersection_point;
      typedef boost::variant<Intersection_point, X_monotone_curve_2>
                                                      Intersection_result;

      const Subcurve_traits_2* geom_traits = m_traits.subcurve_traits_2();
      auto cmp_y_at_x = m_traits.compare_y_at_x_2_object();
      auto equal = geom_traits->equal_2_object();
      auto min_vertex = geom_traits->construct_min_vertex_2_object();
      auto max_vertex = geom_traits->construct_max_vertex_2_object();
      auto intersect = geom_traits->intersect_2_object();
      auto cmp_seg_endpts = geom_traits->compare_endpoints_xy_2_object();
      auto construct_opposite = geom_traits->construct_opposite_2_object();

      X_monotone_subcurve_2 cv10 = cv1[0];
      X_monotone_subcurve_2 cv20 = cv2[0];
      Comparison_result dir1 = cmp_seg_endpts(cv10);
      Comparison_result dir2 = cmp_seg_endpts(cv20);

      const std::size_t n1 = cv1.number_of_subcurves();
      const std::size_t n2 = cv2.number_of_subcurves();

      std::size_t i1 = (dir1 == SMALLER) ? 0 : n1-1;
      std::size_t i2 = (dir2 == SMALLER) ? 0 : n2-1;

      X_monotone_subcurve_2 cv1i1 = cv1[i1];
      X_monotone_subcurve_2 cv2i2 = cv2[i2];
      auto compare_xy = m_traits.compare_xy_2_object();
      Comparison_result left_res =
        compare_xy(cv1i1, ARR_MIN_END, cv2i2, ARR_MIN_END);

      if (left_res == SMALLER) {
        // cv1's left endpoint is to the left of cv2's left endpoint:
        // Locate the index i1 of the subcurve in cv1 which contains cv2's
        // left endpoint.
        i1 = m_traits.locate_impl(cv1, cv2i2, ARR_MIN_END,
                                       Are_all_sides_oblivious_tag());
        if (i1 == Traits::INVALID_INDEX) return oi;
        cv1i1 = cv1[i1];

        if (equal(max_vertex(cv1i1), min_vertex(cv2i2))) {
          if (((dir1 == SMALLER) && (i1 == n1-1)) ||
              ((dir1 == LARGER) && (i1 == 0))){
            // cv1's right endpoint equals cv2's left endpoint
            // Thus we can return this single(!) intersection point
            Intersection_point p(max_vertex(cv1i1), 0);
            *oi++ = Intersection_result(p);
            return oi;
          }
          dir1 == SMALLER ?
            ++i1 :
            (i1 != 0) ? --i1 : (std::size_t) Traits::INVALID_INDEX;
          left_res = EQUAL;
          cv1i1 = cv1[i1];
        }
      }
      else if (left_res == LARGER) {
        // cv1's left endpoint is to the right of cv2's left endpoint:
        // Locate the index i2 of the subcurve in cv2 which contains cv1's
        // left endpoint.
        i2 = m_traits.locate_impl(cv2, cv1i1, ARR_MIN_END,
                                       Are_all_sides_oblivious_tag());
        if (i2 == Traits::INVALID_INDEX) return oi;
        cv2i2 = cv2[i2];

        if (equal(max_vertex(cv2i2), min_vertex(cv1i1))) {
          if (((dir2 == SMALLER) && (i2 == n2-1)) ||
              ((dir2 == LARGER) && (i2 == 0))){
            // cv2's right endpoint equals cv1's left endpoint
            // Thus we can return this single(!) intersection point
            Intersection_point p(max_vertex(cv2i2), 0);
            *oi++ = Intersection_result(p);
            return oi;
          }

          dir2 == SMALLER ?
            ++i2 :
            (i2 != 0) ? --i2 : (std::size_t) Traits::INVALID_INDEX;
          left_res = EQUAL;
          cv2i2 = cv2[i2];
        }
      }

      // Check if the the left endpoint lies on the other polycurve.
      bool left_coincides = (left_res == EQUAL);
      bool left_overlap = false;

      if (left_res == SMALLER)
        left_coincides = (cmp_y_at_x(cv2i2, ARR_MIN_END, cv1i1) == EQUAL);
      else if (left_res == LARGER)
        left_coincides = (cmp_y_at_x(cv1i1, ARR_MIN_END, cv2i2) == EQUAL);

      // NOTE: overlaps are handled in a very inefficient way: each
      // overlapping segment instantiates a new Polyline (with only 2
      // vertices), whereas it should be possible to create a unique
      // Polyline with >2 vertices. However, considering that the
      // caching implementation could only handle that efficiently
      // if curves overlap *with common vertices only*, it seems
      // simpler to keep this current implementation for now.

      // The main loop: Go simultaneously over both polycurves.
      Comparison_result right_res = left_res;
      bool right_coincides = left_coincides;
      bool right_overlap = false;

      while (((dir1 == SMALLER) && (dir2 == SMALLER) &&
              (i1 < n1) && (i2 < n2)) ||
             ((dir1 != SMALLER) && (dir2 == SMALLER) &&
              (i1 != Traits::INVALID_INDEX) && (i2 < n2)) ||
             ((dir1 == SMALLER) && (dir2 != SMALLER) && (i1 < n1) &&
              (i2 != Traits::INVALID_INDEX)) ||
             ((dir1 != SMALLER) && (dir2 != SMALLER) &&
              (i1 != Traits::INVALID_INDEX) &&
              (i2 != Traits::INVALID_INDEX)))
      {
        cv1i1 = cv1[i1];
        cv2i2 = cv2[i2];

        right_res = compare_xy(cv1i1, ARR_MAX_END, cv2i2, ARR_MAX_END);

        right_coincides = (right_res == EQUAL);
        if (right_res == SMALLER)
          right_coincides =
            (cmp_y_at_x(cv1i1, ARR_MAX_END, cv2i2) == EQUAL);
        else if (right_res == LARGER)
          right_coincides =
            (cmp_y_at_x(cv2i2, ARR_MAX_END, cv1i1) == EQUAL);

        right_overlap = false;

        if (! right_coincides && ! left_coincides) {
          oi = intersect(cv1i1, cv2i2, oi);
        }
        else if (right_coincides && left_coincides) {
          right_overlap = true;
          oi = intersect(cv1i1, cv2i2, oi);
        }
        else if (left_coincides && ! right_coincides) {
          if (left_overlap) {
            oi = intersect(cv1i1, cv2i2, oi);
          }
          else {
            if (left_res == SMALLER) {
              Intersection_point p(min_vertex(cv2i2), 0);
              *oi++ = Intersection_result(p);
            }
            else {
              Intersection_point p(min_vertex(cv1i1), 0);
              *oi++ = Intersection_result(p);
            }
          }
        }

        // Proceed forward.
        if (right_res != SMALLER) {
          if (dir2 == SMALLER) ++i2;
          else {
            if (i2 == 0) i2 = Traits::INVALID_INDEX;
            else --i2;
          }
        }
        if (right_res != LARGER) {
          if (dir1 == SMALLER)
            ++i1;
          else {
            if (i1 == 0) i1 = Traits::INVALID_INDEX;
            else --i1;
          }
        }
        left_res = (right_res == SMALLER) ? LARGER :
          (right_res == LARGER) ? SMALLER : EQUAL;

        left_coincides = right_coincides;
      } // END of while loop


      if (right_coincides) {
        typedef std::pair<Point_2,Multiplicity> return_point;
        return_point ip;
        if (right_res == SMALLER) {
          ip = (dir1 == SMALLER) ?
            return_point(max_vertex(cv1[i1-1]), 0) :
            (i1 != Traits::INVALID_INDEX) ?
            return_point(max_vertex(cv1[i1+1]), 0) :
            return_point(max_vertex(cv10), 0);
          *oi++ = Intersection_result(ip);
        }
        else if (right_res == LARGER) {
          ip = (dir2 == SMALLER) ?
            return_point(max_vertex(cv2[i2-1]), 0) :
            (i2 != Traits::INVALID_INDEX) ?
            return_point(max_vertex(cv2[i2+1]), 0) :
            return_point(max_vertex(cv20), 0);
          *oi++ = Intersection_result(ip);
        }
        else if (((i1 > 0) && (dir1 == SMALLER)) ||
                 ((i1 < n1) && (dir1 != SMALLER)) ||
                 ((i1 == Traits::INVALID_INDEX) &&
                  (dir1 != SMALLER)))
        {
          ip = (dir1 == SMALLER) ?
            return_point(max_vertex(cv1[i1-1]), 0) :
            (i1 != Traits::INVALID_INDEX) ?
            return_point(max_vertex(cv1[i1+1]), 0) :
            return_point(max_vertex(cv10), 0);
          *oi++ = Intersection_result(ip);
        }
        else {
          CGAL_assertion_msg((dir2 == SMALLER && i2 > 0) ||
                             (dir2 != SMALLER && i2 < n2) ||
                             (dir2 != SMALLER &&
                              ((i1 == Traits::INVALID_INDEX) ||
                               (i2 == Traits::INVALID_INDEX))),
                             "Wrong index for xcv2 in Intersect_2 of "
                             "polycurves.");
          ip = (dir2 == SMALLER) ?
            return_point(max_vertex(cv2[i2-1]), 0) :
            (i2 != Traits::INVALID_INDEX) ?
            return_point(max_vertex(cv2[i2+1]), 0) :
            return_point(max_vertex(cv20), 0);
          *oi++ = Intersection_result(ip);
        }
      }

      return oi;
    }
  };
  Intersect_2 intersect_2_object() const
  { return Intersect_2(*this); }

  class Are_mergeable_2
  {
  public:
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const
    {
      CGAL_assertion_msg(false, "Are_mergeable_2 not implemented");
      return true;
    }
  };
  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(); }

  class Merge_2
  {
  public:
    void operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2,
                    X_monotone_curve_2& c) const
    {
      CGAL_assertion_msg(false, "Merge_2 not implemented");
    }
  };
  Merge_2 merge_2_object() const
  { return Merge_2(); }

  class Construct_curve_2
  {
  protected:
    using Traits = Arr_caching_polyline_traits_2<Kernel, Range>;
    const Traits& m_traits;
    Construct_curve_2(const Traits& traits) : m_traits(traits) {}
    friend class Arr_caching_polyline_traits_2<Kernel, Range>;
  public:

    Curve_2 operator()(const Range& range, bool duplicate_first = true) const
    {
      const Kernel& kernel = *m_traits.subcurve_traits_2();
      return Curve_2(kernel, range, duplicate_first);
    }
  };

  /*! Obtain a Construct_curve_2 functor object. */
  Construct_curve_2 construct_curve_2_object() const
  { return Construct_curve_2(*this); }
};


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
