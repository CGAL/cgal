#ifndef CGAL_ARR_CACHING_POLYLINE_SUBTRAITS_2_H
#define CGAL_ARR_CACHING_POLYLINE_SUBTRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Arr_geometry_traits/Caching_polyline_2.h>

#include <boost/variant.hpp>

namespace CGAL {

template <typename Kernel_, typename PointIterator>
class Arr_caching_polyline_subtraits_2 : public Kernel_
{
public:
  using Kernel = Kernel_;
  using Base_iterator = PointIterator;
  using Self = Arr_caching_polyline_subtraits_2<Kernel, Base_iterator>;

  using Polyline = internal::X_monotone_caching_polyline_2<Kernel, Base_iterator>;
  using Polyline_iterator = typename Polyline::iterator;

  using FT = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Line_2 = typename Kernel::Line_2;

  using Has_exact_division = typename Algebraic_structure_traits<FT>::Is_exact;

  using Has_left_category = Tag_true;
  using Has_merge_category = Tag_true;
  using Has_do_intersect_category = Tag_false;

  using Left_side_category = Arr_oblivious_side_tag;
  using Bottom_side_category = Arr_oblivious_side_tag;
  using Top_side_category = Arr_oblivious_side_tag;
  using Right_side_category = Arr_oblivious_side_tag;

  using X_monotone_curve_2 = Polyline_iterator;
  using Curve_2 = Polyline_iterator;
  using Multiplicity = unsigned int;


public:
  Arr_caching_polyline_subtraits_2() {}

  class Compare_x_2 {
  protected:
    using Traits = Self;
    const Traits& m_traits;
    Compare_x_2(const Traits& traits) : m_traits(traits) {}
    friend Traits;
  public:
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      const Kernel& kernel = m_traits;
      return (kernel.compare_x_2_object()(p1, p2));
    }
  };
  Compare_x_2 compare_x_2_object() const { return Compare_x_2(*this); }

  class Compare_xy_2
  {
  protected:
    using Traits = Self;
    const Traits& m_traits;
    Compare_xy_2(const Traits& traits) : m_traits(traits) {}
    friend Traits;
  public:
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      const Kernel& kernel = m_traits;
      return (kernel.compare_xy_2_object()(p1, p2));
    }
  };
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(*this); }

  class Construct_min_vertex_2
  {
  public:
    const Point_2& operator()(const X_monotone_curve_2& cv) const
    { return cv.left(); }
  };
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(); }

  class Construct_max_vertex_2
  {
  public:
    const Point_2& operator()(const X_monotone_curve_2& cv) const
    { return cv.right(); }
  };
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(); }

  class Is_vertical_2
  {
  public:
    bool operator()(const X_monotone_curve_2& cv) const
    { return cv.is_vertical(); }
  };
  Is_vertical_2 is_vertical_2_object () const { return Is_vertical_2(); }

  class Compare_y_at_x_2
  {
  protected:
    using Traits = Self;
    const Traits& m_traits;
    Compare_y_at_x_2(const Traits& traits) : m_traits(traits) {}
    friend Traits;
  public:
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& cv) const
    {
      CGAL_precondition(m_traits.is_in_x_range_2_object()(cv, p));

      const Kernel& kernel = m_traits;

      if (! cv.is_vertical()) {
        // Compare p with the segment supporting line.
        CGAL_assertion_code(auto cmp_x = kernel.compare_x_2_object());
        CGAL_assertion(cmp_x(cv.left(), cv.right()) == SMALLER);
        return kernel.orientation_2_object()(cv.left(), cv.right(), p);
      }

      // Compare with the vertical segment endpoints.
      auto compare_y = kernel.compare_y_2_object();
      Comparison_result res1 = compare_y(p, cv.left());
      Comparison_result res2 = compare_y(p, cv.right());
      return (res1 == res2) ? res1 : EQUAL;
    }
  };
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(*this); }

  class Compare_y_at_x_left_2
  {
  protected:
    using Traits = Self;
    const Traits& m_traits;
    Compare_y_at_x_left_2(const Traits& traits) : m_traits(traits) {}
    friend Traits;
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2,
                                 const Point_2& CGAL_assertion_code(p)) const
    {
      const Kernel& kernel = m_traits;

      // Make sure that p lies on both curves, and that both are defined to its
      // left (so their left endpoint is lexicographically smaller than p).
      CGAL_precondition_code(auto compare_xy = kernel.compare_xy_2_object());

      CGAL_precondition((m_traits.compare_y_at_x_2_object()(p, cv1) == EQUAL) &&
                        (m_traits.compare_y_at_x_2_object()(p, cv2) == EQUAL));

      CGAL_precondition(compare_xy(cv1.left(), p) == SMALLER &&
                        compare_xy(cv2.left(), p) == SMALLER);

      const Point_2& src1 = cv1.source();
      const Point_2& tgt1 = cv1.target();
      const Point_2& src2 = cv2.source();
      const Point_2& tgt2 = cv2.target();

      return compare_slopesC2(src2.x(), src2.y(), tgt2.x(), tgt2.y(),
                              src1.x(), src1.y(), tgt1.x(), tgt1.y());

//      return (kernel.compare_slope_2_object()(cv2.line(), cv1.line()));
    }
  };
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(*this); }

  class Compare_y_at_x_right_2
  {
  protected:
    using Traits = Self;
    const Traits& m_traits;
    Compare_y_at_x_right_2(const Traits& traits) : m_traits(traits) {}
    friend Traits;
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2,
                                 const Point_2& CGAL_assertion_code(p)) const
    {
      const Kernel& kernel = m_traits;

      // Make sure that p lies on both curves, and that both are defined to its
      // right (so their right endpoint is lexicographically larger than p).
      CGAL_precondition_code(auto compare_xy = kernel.compare_xy_2_object());

      CGAL_precondition((m_traits.compare_y_at_x_2_object()(p, cv1) == EQUAL) &&
                        (m_traits.compare_y_at_x_2_object()(p, cv2) == EQUAL));

      CGAL_precondition(compare_xy(cv1.right(), p) == LARGER &&
                        compare_xy(cv2.right(), p) == LARGER);

      // Compare slopes without constructing the lines.
      const Point_2& src1 = cv1.source();
      const Point_2& tgt1 = cv1.target();
      const Point_2& src2 = cv2.source();
      const Point_2& tgt2 = cv2.target();

      return compare_slopesC2(src1.x(), src1.y(), tgt1.x(), tgt1.y(),
                              src2.x(), src2.y(), tgt2.x(), tgt2.y());

//      return (kernel.compare_slope_2_object()(cv1.line(), cv2.line()));
    }
  };
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(*this); }

  class Equal_2
  {
  protected:
    using Traits = Self;
    const Traits& m_traits;
    Equal_2(const Traits& traits) : m_traits(traits) {}
    friend Traits;
  public:
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const
    {
      const Kernel& kernel = m_traits;
      typename Kernel::Equal_2  equal = kernel.equal_2_object();

      return (equal(cv1.left(), cv2.left()) &&
              equal(cv1.right(), cv2.right()));
    }

    bool operator()(const Point_2& p1, const Point_2& p2) const
    {
      const Kernel& kernel = m_traits;
      return (kernel.equal_2_object()(p1, p2));
    }
  };
  Equal_2 equal_2_object() const { return Equal_2(*this); }

  class Make_x_monotone_2
  {
  public:
    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const
    {
      // Wrap the segment with a variant.
      typedef boost::variant<Point_2, X_monotone_curve_2>
        Make_x_monotone_result;
      *oi++ = Make_x_monotone_result(cv);
      return oi;
    }
  };
  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(); }

  class Split_2
  {
  protected:
    using Traits = Self;
    const Traits& m_traits;
    Split_2(const Traits& traits) : m_traits(traits) {}
    friend Traits;

  public:
    void operator()(const X_monotone_curve_2& cv, const Point_2& p,
                    X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      CGAL_error_msg ("Split_2 not implemented");
    }
  };
  Split_2 split_2_object() const { return Split_2(*this); }

  class Intersect_2
  {
  protected:
    using Traits = Self;
    const Traits& m_traits;
    Intersect_2(const Traits& traits) : m_traits(traits) {}
    friend Traits;

    bool do_intersect(const Point_2& A1, const Point_2& A2,
                      const Point_2& B1, const Point_2& B2) const
    {
      const Kernel& kernel = m_traits;
      auto compare_xy = kernel.compare_xy_2_object();
      namespace interx = CGAL::Intersections::internal;

      switch(make_certain(compare_xy(A1,B1))) {
       case SMALLER:
        switch(make_certain(compare_xy(A2,B1))) {
         case SMALLER: return false;
         case EQUAL: return true;
         default: // LARGER
          switch(make_certain(compare_xy(A2,B2))) {
           case SMALLER:
            return interx::seg_seg_do_intersect_crossing(A1,A2,B1,B2, kernel);
           case EQUAL: return true;
           default: // LARGER
            return interx::seg_seg_do_intersect_contained(A1,A2,B1,B2, kernel);
          }
        }
       case EQUAL: return true;
       default: // LARGER
        switch(make_certain(compare_xy(B2,A1))) {
         case SMALLER: return false;
         case EQUAL: return true;
         default: // LARGER
          switch(make_certain(compare_xy(B2,A2))) {
           case SMALLER:
            return interx::seg_seg_do_intersect_crossing(B1,B2,A1,A2, kernel);
           case EQUAL: return true;
           default: // LARGER
            return interx::seg_seg_do_intersect_contained(B1,B2,A1,A2, kernel);
          }
        }
      }
      CGAL_error();     // never reached
      return false;
    }

    bool do_bboxes_overlap(const X_monotone_curve_2& cv1,
                           const X_monotone_curve_2& cv2) const
    {
      const Kernel& kernel = m_traits;
      auto construct_bbox = kernel.construct_bbox_2_object();
      auto bbox1 = construct_bbox(cv1.source()) + construct_bbox(cv1.target());
      auto bbox2 = construct_bbox(cv2.source()) + construct_bbox(cv2.target());
      return CGAL::do_overlap(bbox1, bbox2);
    }

  public:
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              OutputIterator oi) const
    {
      typedef std::pair<Point_2, Multiplicity>          Intersection_point;
      typedef boost::variant<Intersection_point, Polyline>
                                                        Intersection_result;

      // Early ending with Bbox overlapping test
      if (! do_bboxes_overlap(cv1, cv2)) return oi;

      // Early ending with specialized do_intersect
      const Kernel& kernel = m_traits;
      if (! do_intersect(cv1.left(), cv1.right(), cv2.left(), cv2.right()))
        return oi;

      // An intersection is guaranteed.

      // Intersect the two supporting lines.
      auto res = kernel.intersect_2_object()(cv1.line(), cv2.line());
      CGAL_assertion(bool(res));

      // Check if we have a single intersection point.
      const Point_2* ip = boost::get<Point_2>(&*res);
      if (ip != nullptr) {
        CGAL_assertion(cv1.is_vertical() ?
                       m_traits.is_in_y_range_2_object()(cv1, *ip) :
                       m_traits.is_in_x_range_2_object()(cv1, *ip));
        CGAL_assertion(cv2.is_vertical() ?
                       m_traits.is_in_y_range_2_object()(cv2, *ip) :
                       m_traits.is_in_x_range_2_object()(cv2, *ip));
        Intersection_point ip_mult(*ip, 1);
        *oi++ = Intersection_result(ip_mult);
        return oi;
      }

      // In this case, the two supporting lines overlap.
      // The overlapping segment is therefore [p_l,p_r], where p_l is the
      // rightmost of the two left endpoints and p_r is the leftmost of the
      // two right endpoints.
      auto compare_xy = kernel.compare_xy_2_object();
      const Point_2& p_l = (compare_xy(cv1.left(), cv2.left()) == SMALLER) ?
        cv2.left() : cv1.left();
      const Point_2& p_r = (compare_xy(cv1.right(), cv2.right()) == SMALLER) ?
        cv1.right() : cv2.right();

      // Examine the resulting segment.
      const Comparison_result cmp_res = compare_xy(p_l, p_r);
      if (cmp_res == EQUAL) {
        // The two segment have the same supporting line, but they just share
        // a common endpoint. Thus we have an intersection point, but we leave
        // the multiplicity of this point undefined.
        Intersection_point ip_mult(p_r, 0);
        *oi++ = Intersection_result(ip_mult);
        return oi;
      }

      CGAL_assertion(cmp_res == SMALLER);

      // We have discovered an overlapping segment, we simply create a
      // new polyline with 2 hanging vertices only:
      if (!cv1.is_directed_right() && !cv2.is_directed_right())
        *oi ++ = Intersection_result(Polyline(kernel, p_r, p_l));
      else
        *oi ++ = Intersection_result(Polyline(kernel, p_l, p_r));
      return oi;
    }
  };
  Intersect_2 intersect_2_object() const { return Intersect_2(*this); }

  class Are_mergeable_2
  {
  protected:
    using Traits = Self;
    const Traits& m_traits;
    Are_mergeable_2(const Traits& traits) : m_traits(traits) {}
    friend Traits;

  public:
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const
    {
      const Kernel& kernel = m_traits;
      typename Kernel::Equal_2 equal = kernel.equal_2_object();
      if (! equal(cv1.right(), cv2.left()) &&
          ! equal(cv2.right(), cv1.left()))
        return false;

      // Check whether the two curves have the same supporting line.
      return (equal(cv1.line(), cv2.line()) ||
              equal(cv1.line(),
                    kernel.construct_opposite_line_2_object()(cv2.line())));
    }
  };
  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(*this); }

  class Merge_2
  {
  protected:
    using Traits = Self;
    const Traits& m_traits;
    Merge_2(const Traits& traits) : m_traits(traits) {}
    friend Traits;

  public:
    void operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2,
                    X_monotone_curve_2& c) const
    {
      CGAL_precondition(m_traits.are_mergeable_2_object()(cv1, cv2));

      const Kernel& kernel = m_traits;
      auto equal = kernel.equal_2_object();

      // Check which curve extends to the right of the other.
      if (equal(cv1.right(), cv2.left())) {
        // cv2 extends cv1 to the right.
        c = cv1;
        c.set_right(cv2.right());
        return;
      }

      CGAL_precondition(equal(cv2.right(), cv1.left()));

      // cv1 extends cv2 to the right.
      c = cv2;
      c.set_right(cv1.right());
    }
  };
  Merge_2 merge_2_object() const { return Merge_2(*this); }

  // using Approximate_number_type = double;
  // class Approximate_2
  // {
  // public:
  //   Approximate_number_type operator()(const Point_2& p, int i) const
  //   {
  //     CGAL_precondition((i == 0) || (i == 1));
  //     return (i == 0) ? (CGAL::to_double(p.x())) : (CGAL::to_double(p.y()));
  //   }
  // };
  // Approximate_2 approximate_2_object() const { return Approximate_2(); }

  // class Construct_x_monotone_curve_2
  // {
  // protected:
  //   using Traits = Self;
  //   const Traits& m_traits;
  //   Construct_x_monotone_curve_2(const Traits& traits) : m_traits(traits) {}
  //   friend Traits;

  // public:
  //   using Segment_2 = typename Kernel::Segment_2;

  //   X_monotone_curve_2 operator()(const Point_2& source,
  //                                 const Point_2& target) const
  //   {
  //     const Kernel& kernel = m_traits;
  //     auto line = kernel.construct_line_2_object()(source, target);
  //     Comparison_result res = kernel.compare_xy_2_object()(source, target);
  //     auto is_degen = (res == EQUAL);
  //     auto is_directed_right = (res == SMALLER);
  //     CGAL_precondition_msg(! is_degen,
  //                           "Cannot construct a degenerate segment.");
  //     auto is_vert = kernel.is_vertical_2_object()(line);
  //     return X_monotone_curve_2(line, source, target,
  //                               is_directed_right, is_vert, is_degen);
  //   }

  //   X_monotone_curve_2 operator()(const Segment_2& seg) const
  //   {
  //     const Kernel& kernel = m_traits;
  //     auto line = kernel.construct_line_2_object()(seg);
  //     auto vertex_ctr = kernel.construct_vertex_2_object();
  //     auto source = vertex_ctr(seg, 0);
  //     auto target = vertex_ctr(seg, 1);
  //     Comparison_result res = kernel.compare_xy_2_object()(source, target);
  //     auto is_degen = (res == EQUAL);
  //     auto is_directed_right = (res == SMALLER);
  //     CGAL_precondition_msg(! is_degen,
  //                           "Cannot construct a degenerate segment.");
  //     auto is_vert = kernel.is_vertical_2_object()(seg);
  //     return X_monotone_curve_2(line, source, target,
  //                               is_directed_right, is_vert, is_degen);
  //   }

  //   X_monotone_curve_2 operator()(const Line_2& line,
  //                                 const Point_2& source,
  //                                 const Point_2& target) const
  //   {
  //     const Kernel& kernel = m_traits;
  //     auto is_vert = kernel.is_vertical_2_object()(line);
  //     Comparison_result res = kernel.compare_xy_2_object()(source, target);
  //     auto is_degen = (res == EQUAL);
  //     auto is_directed_right = (res == SMALLER);
  //     CGAL_precondition_msg(! is_degen,
  //                           "Cannot construct a degenerate segment.");
  //     return X_monotone_curve_2(line, source, target,
  //                               is_directed_right, is_vert, is_degen);
  //   }
  // };
  // Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const
  // { return Construct_x_monotone_curve_2(*this); }

  // class Trim_2
  // {
  //  protected:
  //   using Traits = Self;
  //   const Traits& m_traits;
  //   Trim_2(const Traits& traits) : m_traits(traits) {}
  //   friend Traits;

  // public:
  //   X_monotone_curve_2 operator()(const X_monotone_curve_2& xcv,
  //                                 const Point_2& src,
  //                                 const Point_2& tgt)const
  //   {
  //     CGAL_precondition_code(Equal_2 equal = m_traits.equal_2_object());
  //     CGAL_precondition_code(Compare_y_at_x_2 compare_y_at_x =
  //                            m_traits.compare_y_at_x_2_object());
  //     Compare_x_2 compare_x_2 = m_traits.compare_x_2_object();

  //     // check whether source and taget are two distinct points and they lie
  //     // on the line.
  //     CGAL_precondition(!equal(src, tgt));
  //     CGAL_precondition(compare_y_at_x(src, xcv) == EQUAL);
  //     CGAL_precondition(compare_y_at_x(tgt, xcv) == EQUAL);

  //     // exchange src and tgt IF they do not conform with the direction
  //     X_monotone_curve_2 trimmed_segment;
  //     if (xcv.is_directed_right() && compare_x_2(src, tgt) == LARGER)
  //       trimmed_segment = X_monotone_curve_2(tgt, src);
  //     else if (! xcv.is_directed_right() && (compare_x_2(src, tgt) == SMALLER))
  //       trimmed_segment = X_monotone_curve_2(tgt, src);
  //     else trimmed_segment = X_monotone_curve_2(src, tgt);
  //     return trimmed_segment;
  //   }
  // };
  // Trim_2 trim_2_object() const { return Trim_2(*this); }

  class Compare_endpoints_xy_2
  {
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv) const
    { return (cv.is_directed_right()) ? (SMALLER) : (LARGER); }
  };
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(); }

  class Construct_opposite_2
  {
  public:
    X_monotone_curve_2 operator()(const X_monotone_curve_2& cv) const
    { return cv.opposite(); }
  };
  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(); }

  class Is_in_x_range_2
  {
  protected:
    using Traits = Self;
    const Traits& m_traits;
    Is_in_x_range_2(const Traits& traits) : m_traits(traits) {}
    friend Traits;

  public:
    bool operator()(const X_monotone_curve_2& cv, const Point_2& p) const
    {
      const Kernel& kernel = m_traits;
      auto compare_x = kernel.compare_x_2_object();
      Comparison_result res1 = compare_x(p, cv.left());

      if (res1 == SMALLER) return false;
      else if (res1 == EQUAL) return true;

      Comparison_result res2 = compare_x(p, cv.right());
      return (res2 != LARGER);
    }
  };
  Is_in_x_range_2 is_in_x_range_2_object() const
  { return Is_in_x_range_2(*this); }

  class Is_in_y_range_2
  {
  protected:
    using Traits = Self;
    const Traits& m_traits;
    Is_in_y_range_2(const Traits& traits) : m_traits(traits) {}
    friend Traits;

  public:
    bool operator()(const X_monotone_curve_2& cv, const Point_2& p) const
    {
      const Kernel& kernel = m_traits;
      auto compare_y = kernel.compare_y_2_object();
      Comparison_result res1 = compare_y(p, cv.left());

      if (res1 == SMALLER) return false;
      else if (res1 == EQUAL) return true;

      Comparison_result res2 = compare_y(p, cv.right());
      return (res2 != LARGER);
    }
  };
  Is_in_y_range_2 is_in_y_range_2_object() const
  { return Is_in_y_range_2(*this); }
};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
