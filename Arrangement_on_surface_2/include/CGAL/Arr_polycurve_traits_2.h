// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University(Israel).
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
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>
//                 Ron Wein  <wein@post.tau.ac.il>
//                 Dror Atariah <dror.atariah@fu-berlin.de>
//                 Waqar Khan <wkhan@mpi-inf.mpg.de>

#ifndef CGAL_ARR_POLYCURVE_TRAITS_2_H
#define CGAL_ARR_POLYCURVE_TRAITS_2_H

/*! \file
 * The traits-class for the general piece-wise (polycurve) type of curves of the
 * arrangement package.
 */

#include <iterator>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polycurve_basic_traits_2.h>
#include <CGAL/Arr_geometry_traits/Polycurve_2.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>

namespace CGAL {

template <typename SubcurveTraits_2 = Arr_segment_traits_2<> >
class Arr_polycurve_traits_2 :
    public Arr_polycurve_basic_traits_2<SubcurveTraits_2>
{
public:
  typedef SubcurveTraits_2                                Subcurve_traits_2;

private:
  typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2> Base;

public:
  /// \name Types inherited from the polycurve basic traits class.
  //@{
  typedef typename Base::Has_left_category            Has_left_category;
  typedef typename Base::Has_do_intersect_category
    Has_do_intersect_category;

  typedef typename Base::Left_side_category           Left_side_category;
  typedef typename Base::Bottom_side_category         Bottom_side_category;
  typedef typename Base::Top_side_category            Top_side_category;
  typedef typename Base::Right_side_category          Right_side_category;

  typedef typename Base::Are_all_sides_oblivious_tag
    Are_all_sides_oblivious_tag;

  typedef typename Base::X_monotone_subcurve_2        X_monotone_subcurve_2;
  typedef typename Base::Size                         Size;
  typedef typename Base::size_type                    size_type;

  typedef typename Base::Point_2                      Point_2;
  typedef typename Base::X_monotone_curve_2           X_monotone_curve_2;

  typedef typename Base::Compare_x_2                  Compare_x_2;
  typedef typename Base::Compare_xy_2                 Compare_xy_2;
  typedef typename Base::Construct_min_vertex_2       Construct_min_vertex_2;
  typedef typename Base::Construct_max_vertex_2       Construct_max_vertex_2;
  typedef typename Base::Is_vertical_2                Is_vertical_2;
  typedef typename Base::Compare_y_at_x_2             Compare_y_at_x_2;
  typedef typename Base::Compare_y_at_x_left_2        Compare_y_at_x_left_2;
  typedef typename Base::Compare_y_at_x_right_2       Compare_y_at_x_right_2;
  typedef typename Base::Equal_2                      Equal_2;
  typedef typename Base::Compare_endpoints_xy_2       Compare_endpoints_xy_2;
  typedef typename Base::Construct_opposite_2         Construct_opposite_2;
  typedef typename Base::Approximate_2                Approximate_2;
  typedef typename Base::Construct_x_monotone_curve_2
    Construct_x_monotone_curve_2;
  typedef typename Base::Parameter_space_in_x_2       Parameter_space_in_x_2;
  typedef typename Base::Parameter_space_in_y_2       Parameter_space_in_y_2;
  typedef typename Base::Compare_x_on_boundary_2      Compare_x_on_boundary_2;
  typedef typename Base::Compare_x_at_limit_2         Compare_x_at_limit_2;
  typedef typename Base::Compare_x_near_limit_2       Compare_x_near_limit_2;
  typedef typename Base::Compare_y_on_boundary_2      Compare_y_on_boundary_2;
  typedef typename Base::Compare_y_near_boundary_2    Compare_y_near_boundary_2;
  typedef typename Base::Is_on_y_identification_2     Is_on_y_identification_2;
  typedef typename Base::Is_on_x_identification_2     Is_on_x_identification_2;

  typedef typename Base::Trim_2                       Trim_2;

  //@}

  /// \name Types and functors inherited from the subcurve geometry traits.
  //@{

  typedef typename Subcurve_traits_2::Has_merge_category  Has_merge_category;
  typedef typename Subcurve_traits_2::Multiplicity        Multiplicity;
  typedef typename Subcurve_traits_2::Curve_2             Subcurve_2;

  //@}

  // Backward compatibility:
  typedef Subcurve_2                                      Segment_2;

private:
  typedef Arr_polycurve_traits_2<Subcurve_traits_2>       Self;

public:
  /*! Default constructor */
  Arr_polycurve_traits_2() : Base() {}

  /*! Constructor with given subcurve traits
   * \param seg_traits an already existing subcurve tarits which is passed will
   *        be used by the class.
   */
  Arr_polycurve_traits_2(const Subcurve_traits_2* geom_traits) :
    Base(geom_traits)
  {}

  /*! A polycurve represents a general continuous piecewise-linear
   * curve, without degenerated subcurves.
   */
  typedef internal::Polycurve_2<Subcurve_2, Point_2>      Curve_2;

  /// \name Basic predicate functors(based on the subcurve traits).
  //@{

  /*! \class
   * A functor that obtains the number of points of a polycurve.
   */
  class Number_of_points_2 : public Base::Number_of_points_2 {
  public:
    size_type operator()(const Curve_2& cv) const
    {
      size_type num_seg = cv.number_of_subcurves();
      return (num_seg == 0) ? 0 : num_seg + 1;
    }
  };

  /*! Obtain a number_of_points_2 functor object. */
  Number_of_points_2 number_of_points_2_object() const
  { return Number_of_points_2(); }

  ///@}

  /// \name Construction functors(based on the subcurve traits).
  //@{

  /*! \class
   * A functor that divides an arc into x-monotone arcs. That are, arcs that
   * do not cross the identification arc.
   */
  class Make_x_monotone_2 {
  protected:
    typedef Arr_polycurve_traits_2<Subcurve_traits_2>     Polycurve_traits_2;
    /*! The traits (in case it has state) */
    const Polycurve_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Make_x_monotone_2(const Polycurve_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Cut the given curve into x-monotone sub-curves and insert them into the
     * given output iterator.
     *
     * \pre if `cv` is not empty then it must be continuous and well-oriented.
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. The output
     *           object is a wrapper of a X_monotone_curve_2.
     * \return The past-the-end iterator.
     */
  private:
    template <typename OutputIterator>
    OutputIterator operator_impl(const Curve_2& cv, OutputIterator oi,
                                 Arr_all_sides_oblivious_tag) const
    {
       typedef typename Curve_2::Subcurve_const_iterator const_seg_iterator;

      // If the polycurve is empty, return.
      if (cv.number_of_subcurves() == 0) return oi;

      Construct_x_monotone_curve_2 ctr_x_curve =
        m_poly_traits.construct_x_monotone_curve_2_object();

      typename Subcurve_traits_2::Make_x_monotone_2 make_seg_x_monotone =
        m_poly_traits.subcurve_traits_2()->make_x_monotone_2_object();

      typename Subcurve_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
        m_poly_traits.subcurve_traits_2()->compare_endpoints_xy_2_object();

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      typename Subcurve_traits_2::Construct_opposite_2 ctr_seg_opposite =
        m_poly_traits.subcurve_traits_2()->construct_opposite_2_object();
#endif

      // Convert the input polycurve to a sequence of CGAL objects, such
      // that each Object wraps an x-monotone subcurve.
      std::vector<Object> x_seg_objects;
      const_seg_iterator it_segs;
      for (it_segs = cv.begin_subcurves(); it_segs != cv.end_subcurves();
           ++it_segs)
        make_seg_x_monotone(*it_segs, std::back_inserter(x_seg_objects));
      typename std::vector<Object>::iterator it = x_seg_objects.begin();
      X_monotone_subcurve_2 x_seg;
#if defined (CGAL_NO_ASSERTIONS)
      CGAL::assign(x_seg, *it);
#else
      bool assign_res = CGAL::assign(x_seg, *it);
      CGAL_assertion(assign_res);
#endif

      // If the polycurve consists of a single x-monotone subcurve, return.
      if (x_seg_objects.size() == 1) {
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
        if (cmp_seg_endpts(x_seg) == LARGER)
          x_seg = ctr_seg_opposite(x_seg);
#endif
        *oi++ = make_object(ctr_x_curve(x_seg));
        x_seg_objects.clear();
        return oi;
      }

      CGAL_precondition_code
        (
         // To be used in order to verify continuity and well-orientedness
         // of the input curve cv.
         typename Subcurve_traits_2::Construct_min_vertex_2 min_seg_v =
           m_poly_traits.subcurve_traits_2()->construct_min_vertex_2_object();
         typename Subcurve_traits_2::Construct_max_vertex_2 max_seg_v =
           m_poly_traits.subcurve_traits_2()->construct_max_vertex_2_object();
         typename Subcurve_traits_2::Equal_2 equal =
           m_poly_traits.subcurve_traits_2()->equal_2_object();
         Point_2 last_target = (cmp_seg_endpts(x_seg) == SMALLER) ?
           max_seg_v(x_seg) : min_seg_v(x_seg);
         Point_2 next_src;
         );

      // The polycurve consists of at least 2 x-monotone subcurves:
      Push_back_2 push_back = m_poly_traits.push_back_2_object();
      typename Subcurve_traits_2::Is_vertical_2 is_seg_vertical =
        m_poly_traits.subcurve_traits_2()->is_vertical_2_object();

      bool is_start_vertical = is_seg_vertical(x_seg);
      Comparison_result start_dir = cmp_seg_endpts(x_seg);

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      Push_front_2 push_front = m_poly_traits.push_front_2_object();
      if (cmp_seg_endpts(x_seg) == LARGER) x_seg = ctr_seg_opposite(x_seg);
#endif
      X_monotone_curve_2 x_polycurve = ctr_x_curve(x_seg);

      for (++it; it != x_seg_objects.end(); ++it){
        X_monotone_subcurve_2 x_seg;
#if defined (CGAL_NO_ASSERTIONS)
        CGAL::assign(x_seg, *it);
#else
        bool assign_res = CGAL::assign(x_seg, *it);
        CGAL_assertion(assign_res);
#endif

        // Test that cv is continuous and well-oriented.
        CGAL_precondition_code
          (
           next_src = (cmp_seg_endpts(x_seg) == SMALLER) ?
             min_seg_v(x_seg) : max_seg_v(x_seg);
           );
        CGAL_precondition_msg
          (
           equal(last_target, next_src),
             "cv must form a continuous and well oriented curve."
           );
        CGAL_precondition_code
          (
           last_target = (cmp_seg_endpts(x_seg) == SMALLER) ?
             max_seg_v(x_seg) : min_seg_v(x_seg);
           );

        if ((cmp_seg_endpts(x_seg) != start_dir) ||
            (is_seg_vertical(x_seg) != is_start_vertical))
        {
            // Construct an x-monotone curve from the sub-range which was found
          *oi++ = make_object(x_polycurve);
          is_start_vertical = is_seg_vertical(x_seg);
          start_dir = cmp_seg_endpts(x_seg);
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
          if (cmp_seg_endpts(x_seg) == LARGER) x_seg = ctr_seg_opposite(x_seg);
#endif
          x_polycurve = ctr_x_curve(x_seg);
        }
        else {
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
          if (cmp_seg_endpts(x_seg) == LARGER) {
            x_seg = ctr_seg_opposite(x_seg);
              push_front(x_polycurve, x_seg);
          }
          else
            push_back(x_polycurve, x_seg);
#else
          push_back(x_polycurve, x_seg);
#endif
        }

      } // for loop
      if (x_polycurve.number_of_subcurves() != 0)
        *oi++ = make_object(x_polycurve);
      x_seg_objects.clear();
      return oi;
    }
    template <typename OutputIterator>
    OutputIterator operator_impl(const Curve_2& cv, OutputIterator oi,
                                 Arr_not_all_sides_oblivious_tag) const
    {
      typedef typename Curve_2::Subcurve_const_iterator const_seg_iterator;

      // If the polycurve is empty, return.
      if (cv.number_of_subcurves() == 0) return oi;

      Construct_x_monotone_curve_2 ctr_x_curve =
        m_poly_traits.construct_x_monotone_curve_2_object();

      typename Subcurve_traits_2::Make_x_monotone_2 make_seg_x_monotone =
        m_poly_traits.subcurve_traits_2()->make_x_monotone_2_object();

      typename Subcurve_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
        m_poly_traits.subcurve_traits_2()->compare_endpoints_xy_2_object();

      typename Subcurve_traits_2::Parameter_space_in_x_2 ps_x =
           m_poly_traits.subcurve_traits_2()->parameter_space_in_x_2_object();
      typename Subcurve_traits_2::Parameter_space_in_y_2 ps_y =
           m_poly_traits.subcurve_traits_2()->parameter_space_in_y_2_object();

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      typename Subcurve_traits_2::Construct_opposite_2 ctr_seg_opposite =
        m_poly_traits.subcurve_traits_2()->construct_opposite_2_object();
#endif

      // Convert the input polycurve to a sequence of CGAL objects, such
      // that each Object wraps an x-monotone subcurve.
      std::vector<Object> x_seg_objects;
      const_seg_iterator it_segs;
      for (it_segs = cv.begin_subcurves(); it_segs != cv.end_subcurves();
           ++it_segs)
        make_seg_x_monotone(*it_segs, std::back_inserter(x_seg_objects));
      typename std::vector<Object>::iterator it = x_seg_objects.begin();
      X_monotone_subcurve_2 x_seg;
#if defined (CGAL_NO_ASSERTIONS)
      CGAL::assign(x_seg, *it);
#else
      bool assign_res = CGAL::assign(x_seg, *it);
      CGAL_assertion(assign_res);
#endif

      // If the polycurve consists of a single x-monotone subcurve, return.
      if (x_seg_objects.size() == 1) {
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
        if (cmp_seg_endpts(x_seg) == LARGER)
          x_seg = ctr_seg_opposite(x_seg);
#endif
        *oi++ = make_object(ctr_x_curve(x_seg));
        x_seg_objects.clear();
        return oi;
      }

      CGAL_precondition_code
        (
         // To be used in order to verify continuity and well-orientedness
         // of the input curve cv.
         typename Subcurve_traits_2::Construct_min_vertex_2 min_seg_v =
           m_poly_traits.subcurve_traits_2()->construct_min_vertex_2_object();
         typename Subcurve_traits_2::Construct_max_vertex_2 max_seg_v =
           m_poly_traits.subcurve_traits_2()->construct_max_vertex_2_object();
         typename Subcurve_traits_2::Equal_2 equal =
           m_poly_traits.subcurve_traits_2()->equal_2_object();
         Point_2 last_target = (cmp_seg_endpts(x_seg) == SMALLER) ?
           max_seg_v(x_seg) : min_seg_v(x_seg);
         Point_2 next_src;
         );

      // The polycurve consists of at least 2 x-monotone subcurves:
      Push_back_2 push_back = m_poly_traits.push_back_2_object();
      typename Subcurve_traits_2::Is_vertical_2 is_seg_vertical =
        m_poly_traits.subcurve_traits_2()->is_vertical_2_object();

      bool is_start_vertical = is_seg_vertical(x_seg);
      Comparison_result start_dir = cmp_seg_endpts(x_seg);

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      Push_front_2 push_front = m_poly_traits.push_front_2_object();
      if (cmp_seg_endpts(x_seg) == LARGER) x_seg = ctr_seg_opposite(x_seg);
#endif
      X_monotone_curve_2 x_polycurve = ctr_x_curve(x_seg);

      for (++it; it != x_seg_objects.end(); ++it){
        X_monotone_subcurve_2 x_seg;
#if defined (CGAL_NO_ASSERTIONS)
        CGAL::assign(x_seg, *it);
#else
        bool assign_res = CGAL::assign(x_seg, *it);
        CGAL_assertion(assign_res);
#endif

        // Test that cv is continuous and well-oriented.
        CGAL_precondition_code
          (
           next_src = (cmp_seg_endpts(x_seg) == SMALLER) ?
             min_seg_v(x_seg) : max_seg_v(x_seg);
           );
        CGAL_precondition_msg
          (
           equal(last_target, next_src),
             "cv must form a continuous and well oriented curve."
           );
        CGAL_precondition_code
          (
           last_target = (cmp_seg_endpts(x_seg) == SMALLER) ?
             max_seg_v(x_seg) : min_seg_v(x_seg);
           );

          Arr_curve_end polycurve_target =
            (cmp_seg_endpts(x_polycurve[0]) == SMALLER) ?
            ARR_MAX_END : ARR_MIN_END;
          Arr_curve_end seg_source = (cmp_seg_endpts(x_seg) == SMALLER) ?
            ARR_MIN_END : ARR_MAX_END;
          unsigned int num_segs = x_polycurve.number_of_subcurves();

        if ((cmp_seg_endpts(x_seg) != start_dir) ||
            (is_seg_vertical(x_seg) != is_start_vertical))
        {
            // Construct an x-monotone curve from the sub-range which was found
          *oi++ = make_object(x_polycurve);
          is_start_vertical = is_seg_vertical(x_seg);
          start_dir = cmp_seg_endpts(x_seg);
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
          if (cmp_seg_endpts(x_seg) == LARGER) x_seg = ctr_seg_opposite(x_seg);
#endif
          x_polycurve = ctr_x_curve(x_seg);
        }
        else if (ps_x(x_polycurve[num_segs-1], polycurve_target) !=
                 ARR_INTERIOR ||
                 ps_x(x_seg, seg_source) != ARR_INTERIOR)
        {
          *oi++ = make_object(x_polycurve);
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
          if (cmp_seg_endpts(x_seg) == LARGER) x_seg = ctr_seg_opposite(x_seg);
#endif
            x_polycurve = ctr_x_curve(x_seg);
        }

        else {
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
          if (cmp_seg_endpts(x_seg) == LARGER) {
            x_seg = ctr_seg_opposite(x_seg);
            push_front(x_polycurve, x_seg);
          }
          else
            push_back(x_polycurve, x_seg);
#else
          push_back(x_polycurve, x_seg);
#endif
        }
      } // for loop
      if (x_polycurve.number_of_subcurves() != 0)
        *oi++ = make_object(x_polycurve);
      x_seg_objects.clear();
      return oi;
    }
public:
    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const
    { return operator_impl(cv, oi, Are_all_sides_oblivious_tag()); }
  };

  /*! Obtain a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(*this); }

  /* Functor to augment a polycurve by either adding a vertex or a subcurve
   * at the back.
   * TODO: Test all the operator()'s. (Don't forget vertical cases!)
   */
  class Push_back_2 : public Base::Push_back_2 {
  protected:
    typedef Arr_polycurve_traits_2<Subcurve_traits_2>   Polycurve_traits_2;

  public:
    /*! Constructor. */
    Push_back_2(const Polycurve_traits_2& traits) :
      Base::Push_back_2(traits)
    {}

    // Normally, the moment the compiler finds a name, it stops looking. In
    // other words, the compiler first finds the operator() in the current
    // class and stops looking, never finding the one in the base class.
    // Explicitly bring the base operator() into scope, unnecesitating the
    // code below.
    using Base::Push_back_2::operator();

    // /*! Append a subcurve to an existing x-monotone polycurve at the back.
    //  */
    // void operator()(X_monotone_curve_2& xcv,
    //                 const X_monotone_subcurve_2& seg)
    //   const
    // { Base::Push_back_2::operator()(xcv, seg); }

    /* Append a subcurve to an existing polycurve at the back.
     * If the polycurve is empty, the subcurve will be its only subcurve.
     */
    void operator()(Curve_2& cv, const Subcurve_2& seg) const
    { cv.push_back(seg); }
  };

  /*! Obtain a Push_back_2 functor object. */
  Push_back_2 push_back_2_object() const { return Push_back_2(*this); }

  /* Functor to augment a polycurve by either adding a vertex or a subcurve
   * at the front.
   * TODO: Test all the operator()'s. (Don't forget vertical cases!)
   */
  class Push_front_2 : public Base::Push_front_2 {
  protected:
    typedef Arr_polycurve_traits_2<Subcurve_traits_2>     Polycurve_traits_2;

  public:
    /*! Constructor. */
    Push_front_2(const Polycurve_traits_2& traits) :
      Base::Push_front_2(traits)
    {}

    // Normally, the moment the compiler finds a name, it stops looking. In
    // other words, the compiler first finds the operator() in the current
    // class and stops looking, never finding the one in the base class.
    // Explicitly bring the base operator() into scope, unnecesitating the
    // code below.
    using Base::Push_front_2::operator();

    // /*! Append a subcurve to an existing x-monotone polycurve at the front.
    //  */
    // void operator()(X_monotone_curve_2& xcv,
    //                 const X_monotone_subcurve_2& seg)
    //   const
    // { Base::Push_front_2::operator()(xcv, seg); }

    /* Append a subcurve to an existing polycurve at the front. */
    void operator()(Curve_2& cv, const Subcurve_2& seg) const
    { cv.push_front(seg); }
  };

  /*! Obtain a Push_front_2 functor object. */
  Push_front_2 push_front_2_object() const { return Push_front_2(*this); }

  class Split_2 {
  protected:
    typedef Arr_polycurve_traits_2<Subcurve_traits_2>       Polycurve_traits_2;
    /*! The polycurve traits (in case it has state) */
    const Polycurve_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Split_2(const Polycurve_traits_2& traits) : m_poly_traits(traits) {}

  public:
    /*! Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve(p is its right endpoint).
     * \param c2 Output: The right resulting subcurve(p is its left endpoint).
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator()(const X_monotone_curve_2& xcv, const Point_2& p,
                    X_monotone_curve_2& xcv1, X_monotone_curve_2& xcv2) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      typename Subcurve_traits_2::Construct_min_vertex_2 min_vertex =
        geom_traits->construct_min_vertex_2_object();
      typename Subcurve_traits_2::Construct_max_vertex_2 max_vertex =
        geom_traits->construct_max_vertex_2_object();
      typename Subcurve_traits_2::Equal_2 equal =
        geom_traits->equal_2_object();
      typename Subcurve_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
        geom_traits->compare_endpoints_xy_2_object();

      // Make sure the split point is not one of the curve endpoints.
      CGAL_precondition((!equal(m_poly_traits.
                                construct_min_vertex_2_object()(xcv), p)));
      CGAL_precondition((!equal(m_poly_traits.
                                construct_max_vertex_2_object()(xcv), p)));

      CGAL_precondition_msg(xcv.number_of_subcurves() > 0,
                            "Cannot split a polycurve of length zero.");

      Comparison_result dir = cmp_seg_endpts(xcv[0]);

      // Locate the subcurve on the polycurve xcv that contains p.
      std::size_t i = m_poly_traits.locate(xcv, p);

      CGAL_precondition(i != Polycurve_traits_2::INVALID_INDEX);

      // Clear the output curves.
      xcv1.clear();
      xcv2.clear();

      // Push all subcurves labeled(0, 1, ... , i-1) into xcv1.
      for (std::size_t j = 0; j < i; ++j) xcv1.push_back(xcv[j]);

      if (dir == SMALLER){
        // Check whether the split point is xcv[i]'s source or target.
        if (equal(max_vertex(xcv[i]), p)) {
          // The entire i'th subcurve belongs to xcv1:
          xcv1.push_back(xcv[i]);
        }
        else if (equal(min_vertex(xcv[i]), p)) {
          // The entire i'th subcurves belongs to xcv2:
          xcv2.push_back(xcv[i]);
        }
        else {
          // The i'th subcurve should be split: The left part(seg1)
          // goes to xcv1, and the right part(seg2) goes to xcv2.
          X_monotone_subcurve_2 seg1, seg2;
          m_poly_traits.subcurve_traits_2()->split_2_object()(xcv[i], p,
                                                              seg1, seg2);

          xcv1.push_back(seg1);
          xcv2.push_back(seg2);
        }
      }
      else {
        if (equal(min_vertex(xcv[i]), p)) {
          xcv1.push_back(xcv[i]);
        }
        else if (equal(max_vertex(xcv[i]), p)) {
          xcv2.push_back(xcv[i]);
        }
        else {
          X_monotone_subcurve_2 seg1, seg2;
          m_poly_traits.subcurve_traits_2()->
            split_2_object()(xcv[i], p, seg1, seg2);

          if (cmp_seg_endpts(seg2) == LARGER){
            xcv1.push_back(seg2);
          }
          else {
            // seg2 has to be reversed
            seg2 = m_poly_traits.subcurve_traits_2()->
              construct_opposite_2_object()(seg2);
            xcv1.push_back(seg2);
          }

          if (cmp_seg_endpts(seg1) == LARGER){
            xcv2.push_back(seg1);
          }
          else {
            // seg2 has to be reversed
            seg1 = m_poly_traits.subcurve_traits_2()->
              construct_opposite_2_object()(seg1);
            xcv1.push_back(seg1);
          }
        }
      }

      // Push all subcurves labeled(i+1, i+2, ... , n-1) into xcv1.
      std::size_t n = xcv.number_of_subcurves();

      for (std::size_t j = i+1; j < n; ++j) xcv2.push_back(xcv[j]);

      if (dir != SMALLER) std::swap(xcv1, xcv2);
    }
  };

  /*! Obtain a Split_2 functor object. */
  Split_2 split_2_object() const { return Split_2(*this); }

  class Intersect_2 {
  protected:
    typedef Arr_polycurve_traits_2<Subcurve_traits_2>       Polycurve_traits_2;
    /*! The polycurve traits (in case it has state) */
    const Polycurve_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Intersect_2(const Polycurve_traits_2& traits) : m_poly_traits(traits) {}

    /*! Find the intersections of the two given curves and insert them into the
     * given output iterator. As two subcurves may itersect only once, only a
     * single intersection will be contained in the iterator.
     * Note: If the intersection yields an overlap then it will be oriented
     *       from left-to-right.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator
    operator()(const X_monotone_curve_2& cv1,
               const X_monotone_curve_2& cv2,
               OutputIterator oi) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      Compare_y_at_x_2 cmp_y_at_x = m_poly_traits.compare_y_at_x_2_object();
      typename Subcurve_traits_2::Equal_2 equal =
        geom_traits->equal_2_object();
      typename Subcurve_traits_2::Construct_min_vertex_2 min_vertex =
        geom_traits->construct_min_vertex_2_object();
      typename Subcurve_traits_2::Construct_max_vertex_2 max_vertex =
        geom_traits->construct_max_vertex_2_object();
      typename Subcurve_traits_2::Intersect_2 intersect =
        geom_traits->intersect_2_object();
      typename Subcurve_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
        geom_traits->compare_endpoints_xy_2_object();
      typename Subcurve_traits_2::Construct_opposite_2 construct_opposite =
        geom_traits->construct_opposite_2_object();

      typedef std::pair<Point_2,Multiplicity> Point_2_pair;

      Comparison_result dir1 = cmp_seg_endpts(cv1[0]);
      Comparison_result dir2 = cmp_seg_endpts(cv2[0]);

      const std::size_t n1 = cv1.number_of_subcurves();
      const std::size_t n2 = cv2.number_of_subcurves();

      std::size_t i1 = (dir1 == SMALLER) ? 0 : n1-1;
      std::size_t i2 = (dir2 == SMALLER) ? 0 : n2-1;

      X_monotone_curve_2 ocv;           // Used to represent overlaps.

      Compare_xy_2 compare_xy = m_poly_traits.compare_xy_2_object();
      Comparison_result left_res =
        compare_xy(cv1[i1], ARR_MIN_END, cv2[i2], ARR_MIN_END);

      if (left_res == SMALLER) {
        // cv1's left endpoint is to the left of cv2's left endpoint:
        // Locate the index i1 of the subcurve in cv1 which contains cv2's
        // left endpoint.
        i1 = m_poly_traits.locate_impl(cv1, cv2[i2], ARR_MIN_END,
                                       Are_all_sides_oblivious_tag());
        if (i1 == Polycurve_traits_2::INVALID_INDEX) return oi;

        if (equal(max_vertex(cv1[i1]), min_vertex(cv2[i2]))) {
          if (((dir1 == SMALLER) && (i1 == n1-1)) ||
              ((dir1 == LARGER) && (i1 == 0))){
            // cv1's right endpoint equals cv2's left endpoint
            // Thus we can return this single(!) intersection point
            std::pair<Point_2, Multiplicity>  p(max_vertex(cv1[i1]), 0);
            *oi++ = make_object(p);
            return oi;
          }
          dir1 == SMALLER ?
            ++i1 : (i1 != 0) ? --i1 : (std::size_t) Polycurve_traits_2::INVALID_INDEX;
          left_res = EQUAL;
        }
      }
      else if (left_res == LARGER) {
        // cv1's left endpoint is to the right of cv2's left endpoint:
        // Locate the index i2 of the subcurve in cv2 which contains cv1's
        // left endpoint.
        i2 = m_poly_traits.locate_impl(cv2, cv1[i1], ARR_MIN_END,
                                       Are_all_sides_oblivious_tag());
        if (i2 == Polycurve_traits_2::INVALID_INDEX) return oi;

        if (equal(max_vertex(cv2[i2]), min_vertex(cv1[i1]))) {
          if (((dir2 == SMALLER) && (i2 == n2-1)) ||
              ((dir2 == LARGER) && (i2 == 0))){
            // cv2's right endpoint equals cv1's left endpoint
            // Thus we can return this single(!) intersection point
            std::pair<Point_2, Multiplicity>  p(max_vertex(cv2[i2]), 0);
            *oi++ = make_object(p);
            return oi;
          }

          dir2 == SMALLER ?
            ++i2 : (i2 != 0) ? --i2 : (std::size_t) Polycurve_traits_2::INVALID_INDEX;
          left_res = EQUAL;
        }
      }

      // Check if the the left endpoint lies on the other polycurve.
      bool left_coincides = (left_res == EQUAL);
      bool left_overlap = false;

      if (left_res == SMALLER)
        left_coincides = (cmp_y_at_x(cv2[i2], ARR_MIN_END, cv1[i1]) == EQUAL);
      else if (left_res == LARGER)
        left_coincides = (cmp_y_at_x(cv1[i1], ARR_MIN_END, cv2[i2]) == EQUAL);

      // The main loop: Go simultaneously over both polycurves.
      Comparison_result right_res = left_res;
      bool right_coincides = left_coincides;
      bool right_overlap = false;

      while (((dir1 == SMALLER) && (dir2 == SMALLER) &&
              (i1 < n1) && (i2 < n2)) ||
             ((dir1 != SMALLER) && (dir2 == SMALLER) &&
              (i1 != Polycurve_traits_2::INVALID_INDEX) && (i2 < n2)) ||
             ((dir1 == SMALLER) && (dir2 != SMALLER) && (i1 < n1) &&
              (i2 != Polycurve_traits_2::INVALID_INDEX)) ||
             ((dir1 != SMALLER) && (dir2 != SMALLER) &&
              (i1 != Polycurve_traits_2::INVALID_INDEX) &&
              (i2 != Polycurve_traits_2::INVALID_INDEX)))
      {
        right_res = compare_xy(cv1[i1], ARR_MAX_END, cv2[i2], ARR_MAX_END);

        right_coincides = (right_res == EQUAL);
        if (right_res == SMALLER)
          right_coincides =
            (cmp_y_at_x(cv1[i1], ARR_MAX_END, cv2[i2]) == EQUAL);
        else if (right_res == LARGER)
          right_coincides =
            (cmp_y_at_x(cv2[i2], ARR_MAX_END, cv1[i1]) == EQUAL);

        right_overlap = false;

        if (!right_coincides && !left_coincides) {
          // Non of the endpoints of the current subcurve of one polycurve
          // coincides with the curent subcurve of the other polycurve:
          // Output the intersection if exists.
          oi = intersect(cv1[i1], cv2[i2], oi);
        }
        else if (right_coincides && left_coincides) {
          // An overlap exists between the current subcurves of the
          // polycurves: Output the overlapping subcurve.
          right_overlap = true;

          std::vector<CGAL::Object> int_seg;
          intersect(cv1[i1], cv2[i2], std::back_inserter(int_seg));

          for (size_t i = 0; i < int_seg.size(); ++i) {
            const X_monotone_subcurve_2* x_seg =
              CGAL::object_cast<X_monotone_subcurve_2> (&(int_seg[i]));
            if (x_seg != NULL) {
              X_monotone_subcurve_2 seg = *x_seg;

              // If for some reason the subcurve intersection
              // results in left oriented curve.
              if ( cmp_seg_endpts(seg) == LARGER)
                seg = construct_opposite(seg);
              ocv.push_back(seg);
            }

            const Point_2_pair* p_ptr =
              CGAL::object_cast<Point_2_pair>(&(int_seg[i]));
            if (p_ptr != NULL) {
              // Any point that is not equal to the max_vertex of the
              // subcurve should be inserted into oi.
              // The max_vertex of the current subcurve (if intersecting)
              // will be taken care of as the min_vertex of in the next
              // iteration.
              if (!equal(p_ptr->first, max_vertex(cv1[i1])))
                *oi++ = make_object(*p_ptr);
            }
          }
        }

        else if (left_coincides && !right_coincides) {
          // std::cout << "Left is coinciding but right is not." << std::endl;
          // The left point of the current subcurve of one polycurve
          // coincides with the current subcurve of the other polycurve.
          if (left_overlap) {
            // An overlap occured at the previous iteration:
            // Output the overlapping polycurve.
            CGAL_assertion(ocv.number_of_subcurves() > 0);
            *oi++ = make_object(ocv);
            ocv.clear();
          }
          else {
            // The left point of the current subcurve of one
            // polycurve coincides with the current subcurve of the
            // other polycurve, and no overlap occured at the
            // previous iteration: Output the intersection
            // point. The derivative of at least one of the
            // polycurves is not defined at this point, so we give
            // it multiplicity 0.
            if (left_res == SMALLER) {
              std::pair<Point_2, Multiplicity> p(min_vertex(cv2[i2]), 0);
              *oi++ = make_object(p);
            }
            else {
              std::pair<Point_2, Multiplicity> p(min_vertex(cv1[i1]), 0);
              *oi++ = make_object(p);
            }
          }
        }

        // Proceed forward.
        if (right_res != SMALLER) {
          if (dir2 == SMALLER) ++i2;
          else {
            if (i2 == 0) i2 = Polycurve_traits_2::INVALID_INDEX;
            else --i2;
          }
        }
        if (right_res != LARGER) {
          if (dir1 == SMALLER)
            ++i1;
          else {
            if (i1 == 0) i1 = Polycurve_traits_2::INVALID_INDEX;
            else --i1;
          }
        }
        left_res = (right_res == SMALLER) ? LARGER :
          (right_res == LARGER) ? SMALLER : EQUAL;

        left_coincides = right_coincides;
        left_overlap = right_overlap;
      } // END of while loop

        // Output the remaining overlapping polycurve, if necessary.
      if (ocv.number_of_subcurves() > 0) {
        *oi++ = make_object(ocv);
      }
      else if (right_coincides) {
        typedef std::pair<Point_2,Multiplicity> return_point;
        return_point ip;
        if (right_res == SMALLER) {
          ip = (dir1 == SMALLER) ?
            return_point(max_vertex(cv1[i1-1]), 0) :
            (i1 != Polycurve_traits_2::INVALID_INDEX) ?
            return_point(max_vertex(cv1[i1+1]), 0) :
            return_point(max_vertex(cv1[0]), 0);
          *oi++ = make_object(ip);
        }
        else if (right_res == LARGER) {
          ip = (dir2 == SMALLER) ?
            return_point(max_vertex(cv2[i2-1]), 0) :
            (i2 != Polycurve_traits_2::INVALID_INDEX) ?
            return_point(max_vertex(cv2[i2+1]), 0) :
            return_point(max_vertex(cv2[0]), 0);
          *oi++ = make_object(ip);
        }
        else if (((i1 > 0) && (dir1 == SMALLER)) ||
                 ((i1 < n1) && (dir1 != SMALLER)) ||
                 ((i1 == Polycurve_traits_2::INVALID_INDEX) &&
                  (dir1 != SMALLER)))
        {
          ip = (dir1 == SMALLER) ?
            return_point(max_vertex(cv1[i1-1]), 0) :
            (i1 != Polycurve_traits_2::INVALID_INDEX) ?
            return_point(max_vertex(cv1[i1+1]), 0) :
            return_point(max_vertex(cv1[0]), 0);
          *oi++ = make_object(ip);
        }
        else {
          CGAL_assertion_msg((dir2 == SMALLER && i2 > 0) ||
                             (dir2 != SMALLER && i2 < n2) ||
                             (dir2 != SMALLER &&
                              ((i1 == Polycurve_traits_2::INVALID_INDEX) ||
                               (i2 == Polycurve_traits_2::INVALID_INDEX))),
                             "Wrong index for xcv2 in Intersect_2 of "
                             "polycurves.");
          ip = (dir2 == SMALLER) ?
            return_point(max_vertex(cv2[i2-1]), 0) :
            (i2 != Polycurve_traits_2::INVALID_INDEX) ?
            return_point(max_vertex(cv2[i2+1]), 0) :
            return_point(max_vertex(cv2[0]), 0);
          *oi++ = make_object(ip);
        }
      }

      return oi;
    }
  };

  /*! Obtain an Intersect_2 functor object. */
  Intersect_2 intersect_2_object() const
  { return Intersect_2(*this); }

  class Are_mergeable_2 {
  protected:
    typedef Arr_polycurve_traits_2<Subcurve_traits_2>       Polycurve_traits_2;
    /*! The polycurve traits (in case it has state) */
    const Polycurve_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Are_mergeable_2(const Polycurve_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return(true) if the two curves are mergeable, that is, they share a
     * common endpoint and the same orientation;(false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      Construct_min_vertex_2 min_vertex =
        m_poly_traits.construct_min_vertex_2_object();
      Construct_max_vertex_2 max_vertex =
        m_poly_traits.construct_max_vertex_2_object();
      typename Subcurve_traits_2::Equal_2 equal =
        geom_traits->equal_2_object();
      typename Subcurve_traits_2::Is_vertical_2 is_seg_vertical =
        geom_traits->is_vertical_2_object();

      Comparison_result dir1 =
        m_poly_traits.compare_endpoints_xy_2_object()(cv1);
      Comparison_result dir2 =
        m_poly_traits.compare_endpoints_xy_2_object()(cv2);

      if (dir1 != dir2)
        return false;

      bool ver1 = is_seg_vertical(cv1[0]);
      bool ver2 = is_seg_vertical(cv2[0]);

      return (((// Both are directed from left-to-right
                (dir1 == SMALLER) &&
                ((equal(max_vertex(cv1),min_vertex(cv2))) ||
                 (equal(max_vertex(cv2),min_vertex(cv1))))) ||
               (// Both are directed from right-to-left
                (dir1 == LARGER) &&
                ((equal(min_vertex(cv1),max_vertex(cv2))) ||
                 (equal(max_vertex(cv1),min_vertex(cv2))))
                )) &&
              (// Either both should be vertical or both should
               // be NOT vertical.
               (ver1 && ver2) || (!ver1 && !ver2)));
    }
  };

  /*! Obtain an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(*this); }

  /*! \class Merge_2
   * A functor that merges two x-monotone curves into one.
   */
  /* Roadmap: Allow merging of overlapping polycurves. This means also
   *          changing the subcurve traits class.
   */
  class Merge_2 {
  protected:
    typedef Arr_polycurve_traits_2<Subcurve_traits_2>     Geometry_traits;
    /*! The traits (in case it has state) */
    const Geometry_traits& m_poly_traits;

  public:
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Merge_2(const Geometry_traits& traits) : m_poly_traits(traits) {}

    /*! Merge two given x-monotone curves into a single curve(segment).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable.
     */
    void operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2,
                    X_monotone_curve_2& c) const
    {
      CGAL_precondition(m_poly_traits.are_mergeable_2_object()(cv1, cv2));

      Construct_min_vertex_2 get_min_v =
        m_poly_traits.construct_min_vertex_2_object();
      Construct_max_vertex_2 get_max_v =
        m_poly_traits.construct_max_vertex_2_object();
      Compare_endpoints_xy_2 cmp_seg_endpts =
        m_poly_traits.compare_endpoints_xy_2_object();
      Equal_2 equal = m_poly_traits.equal_2_object();

      c.clear();
      if (// Either both are left-to-right and cv2 is to the right of cv1
          ((cmp_seg_endpts(cv1)==SMALLER) &&
           (equal(get_max_v(cv1),get_min_v(cv2)))) ||
          // or both are right-to-left and cv2 is to the left of cv1
          ((cmp_seg_endpts(cv1)==LARGER) &&
           (equal(get_min_v(cv1), get_max_v(cv2)))))
      {
        const std::size_t n1 = cv1.number_of_subcurves();
        const std::size_t n2 = cv2.number_of_subcurves();
        std::size_t i;

        // cv2 extends cv1 to the right:
        for (i = 0; i < n1 - 1; ++i) c.push_back(cv1[i]);

        // Try to merge the to contiguous line subcurves:
        if (m_poly_traits.subcurve_traits_2()->
            are_mergeable_2_object()(cv1[n1 - 1], cv2[0]))
        {
          X_monotone_subcurve_2 seg;
          m_poly_traits.subcurve_traits_2()->
            merge_2_object()(cv1[n1 - 1], cv2[0], seg);
          c.push_back(seg);
        }
        else {
          c.push_back(cv1[n1 - 1]);
          c.push_back(cv2[0]);
        }

        for (i = 1; i < n2; ++i) c.push_back(cv2[i]);
      }
      else
        return this->operator()(cv2,cv1,c);
    }
  };

  /*! Obtain a Merge_2 functor object. */
  Merge_2 merge_2_object() const { return Merge_2(*this); }
  ///@}

  /*! \class
   * A functor that constructs a (general) polycurve.
   */
  class Construct_curve_2 {
  protected:
    typedef Arr_polycurve_traits_2<Subcurve_traits_2>       Polycurve_traits_2;
    /*! The polycurve traits (in case it has state) */
    const Polycurve_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Construct_curve_2(const Polycurve_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Obtain a polycurve that consists of one given subcurve. */
    Curve_2 operator()(const Subcurve_2& seg) const { return Curve_2(seg); }

    /* Construct a well-oriented polycurve from a range of either
     * `SubcurveTraits::Point_2` or `SubcurveTraits::Subcurve_2`.
     */
    template <typename ForwardIterator>
    Curve_2 operator()(ForwardIterator begin, ForwardIterator end) const
    {
      typedef typename std::iterator_traits<ForwardIterator>::value_type VT;
      typedef typename boost::is_same<VT, Point_2>::type Is_point;
      // Dispatch the range to the appropriate implementation.
      return constructor_impl(begin, end, Is_point());
    }

    /*! Construction of a polycurve from a range of points.
     * \pre The range contains at least two points
     * \pre Consecutive points are disjoint.
     * \return Well-oriented polycurve connecting the given
     *         points. The order of the vertices is determined by
     *         their order in the range.  Furthermore, the
     *         orientation of the polycurve is induced by their
     *         order.
     */
    template <typename ForwardIterator>
    Curve_2 constructor_impl(ForwardIterator /* begin */,
                             ForwardIterator /* end */,
                             boost::true_type) const
    {  CGAL_error_msg("Cannot construct a polycurve from a range of points!"); }

    /*! Construction implementation from a range of subcurves.
     *  Note that the subcurves in the range are NOT necessarily x-monotone,
     *  thus it is impossible to test (even in precondition) whether the input
     *  forms a continuous and well oriented polycurve.
     *  \pre Range should contain at least one subcurve.
     */
    template <typename ForwardIterator>
    Curve_2 constructor_impl(ForwardIterator begin, ForwardIterator end,
                             boost::false_type) const
    {
      // Range has to contain at least one subcurve
      CGAL_precondition(begin != end);
      return Curve_2(begin, end);
    }
  };

  /*! Obtain a Construct_curve_2 functor object. */
  Construct_curve_2 construct_curve_2_object() const
  { return Construct_curve_2(*this); }
};

} //namespace CGAL

#endif
