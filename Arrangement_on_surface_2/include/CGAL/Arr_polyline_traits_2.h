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
// $URL$
// $Id$
//
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>
//                 Ron Wein  <wein@post.tau.ac.il>
//                 Dror Atariah <dror.atariah@fu-berlin.de>

#ifndef CGAL_ARR_POLYLINE_TRAITS_2_H
#define CGAL_ARR_POLYLINE_TRAITS_2_H

/*! \file
 * The traits-class for the linear piece-wiese(polyline) type of curves of the
 * arrangement package.
 */

#include <iterator>
#include <boost/type_traits/is_same.hpp>

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_geometry_traits/Polyline_2.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>

namespace CGAL {

  template <typename SegmentTraits_2 = Arr_segment_traits_2<> >
  class Arr_polyline_traits_2 {
  public:
    typedef SegmentTraits_2                          Segment_traits_2;

    // Tag definitions:
    typedef Tag_true                                 Has_left_category;
    typedef Tag_true                                 Has_merge_category;
    typedef Tag_false                                Has_do_intersect_category;

    typedef typename Segment_traits_2::Left_side_category   Left_side_category;
    typedef typename Segment_traits_2::Bottom_side_category Bottom_side_category;
    typedef typename Segment_traits_2::Top_side_category    Top_side_category;
    typedef typename Segment_traits_2::Right_side_category  Right_side_category;

    typedef typename Arr_are_all_sides_oblivious_tag
    <Left_side_category, Bottom_side_category,
     Top_side_category, Right_side_category>::result
    Are_all_sides_oblivious_tag;

  private:
    typedef Arr_polyline_traits_2<Segment_traits_2>  Self;

    // Data members:
    const Segment_traits_2* m_seg_traits;    // The base segment-traits class.
    bool m_own_traits;

  private:
    enum { INVALID_INDEX = 0xffffffff };

  public:
    /*! Default constructor */
    Arr_polyline_traits_2() :
      m_seg_traits(new Segment_traits_2()), m_own_traits(true) {}

    /*! Constructor with given segment traits
     * \param seg_traits an already existing segment tarits which is passed will
     *        be used by the class.
     */
    Arr_polyline_traits_2(const Segment_traits_2* seg_traits) :
      m_seg_traits(seg_traits), m_own_traits(false){ }

    /* Destructor
     * Deletes the segment tarits class in case it was constructed during the
     * construction of this.
     */
    ~Arr_polyline_traits_2(){
      if (m_own_traits)
        delete m_seg_traits;
    }

    /*! Obtain the segment traits.
     * \return the segment traits.
     */
    const Segment_traits_2* segment_traits_2() const { return m_seg_traits; }

    /// \name Types and functors inherited from the base segment traits.
    //@{

    // Traits types:
    typedef typename Segment_traits_2::Point_2            Point_2;
    typedef typename Segment_traits_2::X_monotone_curve_2 X_monotone_segment_2;
    typedef typename Segment_traits_2::Curve_2            Segment_2;

    /*!
     * A polyline represents a general continuous piecewise-linear
     * curve, without degenerated segments.
     */
    typedef polyline::Polyline_2<Segment_2, Point_2>      Curve_2;
    /*!
     * An x monotone polyline represents a continuous piecewise-linear
     * curve which is either strongly x-monotone or vertical. Again,
     * the polyline is without degenerated segments.
     */
    typedef polyline::X_monotone_polyline_2<X_monotone_segment_2, Point_2>
      X_monotone_curve_2;

    typedef typename Segment_traits_2::Multiplicity       Multiplicity;

    /*! Compare the x-coordinates of two points. */
    class Compare_x_2 {
    private:
      // Oblivious implementation
      Comparison_result operator()(const X_monotone_segment_2& xs1,
                                   Arr_curve_end ce1,
                                   const X_monotone_segment_2& xs2,
                                   Arr_curve_end ce2,
                                   Arr_all_sides_oblivious_tag) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        const Point_2& p1 = (ce1 == ARR_MAX_END) ?
          seg_traits->construct_max_vertex_2_object()(xs1) :
          seg_traits->construct_min_vertex_2_object()(xs1);
        const Point_2& p2 = (ce2 == ARR_MAX_END) ?
          seg_traits->construct_max_vertex_2_object()(xs2) :
          seg_traits->construct_min_vertex_2_object()(xs2);
        return seg_traits->compare_x_2_object()(p1, p2);
      }

      // Boundary implementation
      Comparison_result operator()(const X_monotone_segment_2& xs1,
                                   Arr_curve_end ce1,
                                   const X_monotone_segment_2& xs2,
                                   Arr_curve_end ce2,
                                   Arr_not_all_sides_oblivious_tag) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        typename Segment_traits_2::Parameter_space_in_x_2 ps_x =
          seg_traits->parameter_space_in_x_2_object();
        const Arr_parameter_space ps_x1 = ps_x(xs1, ce1);
        const Arr_parameter_space ps_x2 = ps_x(xs2, ce2);

        if (ps_x1 != ps_x2) {
          if (ps_x1 == ARR_LEFT_BOUNDARY) return LARGER;
          if (ps_x1 == ARR_RIGHT_BOUNDARY) return SMALLER;
          if (ps_x2 == ARR_LEFT_BOUNDARY) return SMALLER;
          if (ps_x2 == ARR_RIGHT_BOUNDARY) return LARGER;
        }

        // ps_x1 == ps_x2
        if (ps_x1 != ARR_INTERIOR) return EQUAL;

        typename Segment_traits_2::Parameter_space_in_y_2 ps_y =
          seg_traits->parameter_space_in_y_2_object();
        const Arr_parameter_space ps_y1 = ps_y(xs1, ce1);
        const Arr_parameter_space ps_y2 = ps_y(xs2, ce2);
        if (ps_y1 == ARR_INTERIOR) {
          const Point_2& p1 = (ce1 == ARR_MAX_END) ?
            seg_traits->construct_max_vertex_2_object()(xs1) :
            seg_traits->construct_min_vertex_2_object()(xs1);
          if (ps_y2 == ARR_INTERIOR) {
            const Point_2& p2 = (ce2 == ARR_MAX_END) ?
              seg_traits->construct_max_vertex_2_object()(xs2) :
              seg_traits->construct_min_vertex_2_object()(xs2);
            return seg_traits->compare_x_2_object()(p1, p2);
          }
          typename Segment_traits_2::Compare_x_on_boundary_2 cmp_x_on_bnd =
            seg_traits->compare_x_on_boundary_2_object();
          return cmp_x_on_bnd(p1, xs2, ce2);
        }
        if (ps_y2 == ARR_INTERIOR) {
          const Point_2& p2 = (ce2 == ARR_MAX_END) ?
            seg_traits->construct_max_vertex_2_object()(xs2) :
            seg_traits->construct_min_vertex_2_object()(xs2);
          typename Segment_traits_2::Compare_x_on_boundary_2 cmp_x_on_bnd =
            seg_traits->compare_x_on_boundary_2_object();
          return opposite(cmp_x_on_bnd(p2, xs1, ce1));
        }
        return seg_traits->compare_x_on_boundary_2_object()(xs1, ce1, xs2, ce2);
      }

    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Compare_x_2(const Polyline_traits_2& traits) : m_poly_traits(traits) {}

      /*! Compare the x-coordinates of two directional points.
       * \param p1 the first directional point.
       * \param p2 the second directional point.
       * \return SMALLER - x(p1) < x(p2);
       *         EQUAL   - x(p1) = x(p2);
       *         LARGER  - x(p1) > x(p2).
       * \pre p1 does not lie on the boundary.
       * \pre p2 does not lie on the boundary.
       */
      Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
      { return m_poly_traits.segment_traits_2()->compare_x_2_object()(p1, p2); }

      /*! Compare two ends of x-monotone curves in x.
       * \param xs1 the first curve.
       * \param ce1 the curve-end indicator of the first x-monotone curve xs1:
       *            ARR_MIN_END - the minimal end of xs1 or
       *            ARR_MAX_END - the maximal end of xs1.
       * \param xs2 the second curve.
       * \param ce2 the curve-end indicator of the second x-monoton curve xs2:
       *            ARR_MIN_END - the minimal end of xs2 or
       *            ARR_MAX_END - the maximal end of xs2.
       */
      Comparison_result operator()(const X_monotone_segment_2& xs1,
                                   Arr_curve_end ce1,
                                   const X_monotone_segment_2& xs2,
                                   Arr_curve_end ce2)
      { return operator()(xs1, ce1, xs2, ce2, Are_all_sides_oblivious_tag()); }
    };

    /*! Get a Compare_x_2 functor object. */
    Compare_x_2 compare_x_2_object() const
    { return Compare_x_2(*this); }

    /*! Compare two curve-ends or points lexigoraphically: by x, then by y. */
    class Compare_xy_2 {
    private:
      // Oblivious implementation
      Comparison_result operator()(const X_monotone_segment_2& xs1,
                                   Arr_curve_end ce1,
                                   const X_monotone_segment_2& xs2,
                                   Arr_curve_end ce2,
                                   Arr_all_sides_oblivious_tag) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        const Point_2& p1 = (ce1 == ARR_MAX_END) ?
          seg_traits->construct_max_vertex_2_object()(xs1) :
          seg_traits->construct_min_vertex_2_object()(xs1);
        const Point_2& p2 = (ce2 == ARR_MAX_END) ?
          seg_traits->construct_max_vertex_2_object()(xs2) :
          seg_traits->construct_min_vertex_2_object()(xs2);
        return seg_traits->compare_xy_2_object()(p1, p2);
      }

      // Boundary implementation
      Comparison_result operator()(const X_monotone_segment_2& xs1,
                                   Arr_curve_end ce1,
                                   const X_monotone_segment_2& xs2,
                                   Arr_curve_end ce2,
                                   Arr_not_all_sides_oblivious_tag) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        typename Segment_traits_2::Parameter_space_in_x_2 ps_x =
          seg_traits->parameter_space_in_x_2_object();
        typename Segment_traits_2::Parameter_space_in_y_2 ps_y =
          seg_traits->parameter_space_in_y_2_object();
        const Arr_parameter_space ps_x1 = ps_x(xs1, ce1);
        const Arr_parameter_space ps_y1 = ps_y(xs1, ce1);
        const Arr_parameter_space ps_x2 = ps_x(xs2, ce2);
        const Arr_parameter_space ps_y2 = ps_y(xs2, ce2);

        if (ps_x1 != ps_x2) {
          if (ps_x1 == ARR_LEFT_BOUNDARY) return LARGER;
          if (ps_x1 == ARR_RIGHT_BOUNDARY) return SMALLER;
          if (ps_x2 == ARR_LEFT_BOUNDARY) return SMALLER;
          if (ps_x2 == ARR_RIGHT_BOUNDARY) return LARGER;
        }

        if ((ps_x1 == ARR_INTERIOR) && (ps_y1 == ARR_INTERIOR)) {
          const Point_2& p1 = (ce1 == ARR_MAX_END) ?
            seg_traits->construct_max_vertex_2_object()(xs1) :
            seg_traits->construct_min_vertex_2_object()(xs1);
          // ps1 == ARR_INTERIOR

          if ((ps_x2 == ARR_INTERIOR) && (ps_y2 == ARR_INTERIOR)) {
            const Point_2& p2 = (ce2 == ARR_MAX_END) ?
              seg_traits->construct_max_vertex_2_object()(xs2) :
              seg_traits->construct_min_vertex_2_object()(xs2);

            // ps1 == ARR_INTERIOR
            // ps2 == ARR_INTERIOR
            return seg_traits->compare_xy_2_object()(p1, p2);
          }

          // The cases ps_x2 == ARR_{LEFT,RIGHT}_BOUNDARY are handled above

          // ps1 == ARR_INTERIOR
          // ps_x2 == ARR_INTERIOR
          // ps_y2 != ARR_INTERIOR
          CGAL_assertion(ps_x2 == ARR_INTERIOR);
          // EFEF: missing implementation for open boundary.
          typename Segment_traits_2::Compare_x_on_boundary_2 cmp_x_on_bnd =
            seg_traits->compare_x_on_boundary_2_object();
          Comparison_result res = cmp_x_on_bnd(p1, xs2, ce2);
          if (res != EQUAL) return res;
          if (ps_y2 == ARR_TOP_BOUNDARY) return SMALLER;
          CGAL_assertion(ps_y2 == ARR_BOTTOM_BOUNDARY);
          return LARGER;
        }

        // ps1 != ARR_INTERIOR
        if ((ps_x2 == ARR_INTERIOR) && (ps_y2 == ARR_INTERIOR)) {
          const Point_2& p2 = (ce2 == ARR_MAX_END) ?
            seg_traits->construct_max_vertex_2_object()(xs2) :
            seg_traits->construct_min_vertex_2_object()(xs2);

          // The cases ps_x1 == ARR_{LEFT,RIGHT}_BOUNDARY are handled above

          // ps_x1 == ARR_INTERIOR
          // ps_y1 != ARR_INTERIOR
          // ps2 == ARR_INTERIOR
          CGAL_assertion(ps_x1 == ARR_INTERIOR);
          typename Segment_traits_2::Compare_x_on_boundary_2 cmp_x_on_bnd =
            seg_traits->compare_x_on_boundary_2_object();
          Comparison_result res = cmp_x_on_bnd(p2, xs1, ce1);
          if (res != EQUAL) return opposite(res);
          if (ps_y1 == ARR_TOP_BOUNDARY) return LARGER;
          CGAL_assertion(ps_y1 == ARR_BOTTOM_BOUNDARY);
          return SMALLER;
        }

        // ps1 != ARR_INTERIOR
        // ps2 != ARR_INTERIOR
        // ps_x1 == ps_x2
        if (ps_x1 == ARR_INTERIOR) {
          // ps_y1 != ARR_INTERIOR
          // ps_y2 != ARR_INTERIOR
          Comparison_result res =
            seg_traits->compare_x_on_boundary_2_object()(xs1, ce1, xs2, ce2);
          if (res != EQUAL) return res;
          if (ps_y1 == ps_y2) return EQUAL;
          return (ps_y1 == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER;
        }

        CGAL_assertion(ce1 == ce2);
        const Point_2& p1 = (ce1 == ARR_MAX_END) ?
          seg_traits->construct_max_vertex_2_object()(xs1) :
          seg_traits->construct_min_vertex_2_object()(xs1);
        const Point_2& p2 = (ce2 == ARR_MAX_END) ?
          seg_traits->construct_max_vertex_2_object()(xs2) :
          seg_traits->construct_min_vertex_2_object()(xs2);
        return seg_traits->compare_y_on_boundary_2_object()(p1, p2);
      }

    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Compare_xy_2(const Polyline_traits_2& traits) : m_poly_traits(traits) {}

      /*! Compare two directional points lexigoraphically: by x, then by y.
       * \param p1 the first enpoint directional point.
       * \param p2 the second endpoint directional point.
       * \return SMALLER - x(p1) < x(p2);
       *         SMALLER - x(p1) = x(p2) and y(p1) < y(p2);
       *         EQUAL   - x(p1) = x(p2) and y(p1) = y(p2);
       *         LARGER  - x(p1) = x(p2) and y(p1) > y(p2);
       *         LARGER  - x(p1) > x(p2).
       * \pre p1 does not lie on the boundary.
       * \pre p2 does not lie on the boundary.
       */
      Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
      { return m_poly_traits.segment_traits_2()->compare_xy_2_object()(p1, p2); }

      /*! Compare two ends of x-monotone curves lexicographically.
       * \param xs1 the first curve.
       * \param ce1 the curve-end indicator of the first x-monotone curve xs1:
       *            ARR_MIN_END - the minimal end of xs1 or
       *            ARR_MAX_END - the maximal end of xs1.
       * \param xs2 the second curve.
       * \param ce2 the curve-end indicator of the second x-monoton curve xs2:
       *            ARR_MIN_END - the minimal end of xs2 or
       *            ARR_MAX_END - the maximal end of xs2.
       */
      Comparison_result operator()(const X_monotone_segment_2& xs1,
                                   Arr_curve_end ce1,
                                   const X_monotone_segment_2& xs2,
                                   Arr_curve_end ce2)
      { return operator()(xs1, ce1, xs2, ce2, Are_all_sides_oblivious_tag()); }
    };

    /*! Get a Compare_xy_2 functor object. */
    Compare_xy_2 compare_xy_2_object() const
    { return Compare_xy_2(*this); }

    ///@}

    /// \name Basic predicate functors(based on the segment traits).
    //@{

    class Number_of_points_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /* The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /* Constructor. */
      Number_of_points_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits)
      {}

      int operator()(const Curve_2& cv) const
      {
        int num_seg = cv.number_of_segments();
        return (num_seg == 0) ? 0 : num_seg + 1;
      }
    };

    Number_of_points_2 number_of_points_2_object() const
    { return Number_of_points_2(this); }


    class Construct_min_vertex_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /* Constructor. */
      Construct_min_vertex_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits)
      {}

      /*!
       * Get the left endpoint of the x-monotone curve(segment).
       * \param cv The polyline curve.
       * \return The left endpoint.
       */
      const Point_2& operator()(const X_monotone_curve_2& cv) const
      {
        CGAL_assertion(cv.number_of_segments() > 0);

        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();

        if (seg_traits->compare_endpoints_xy_2_object()(cv[0]) == SMALLER)
          return seg_traits->construct_min_vertex_2_object()(cv[0]);
        else
          return seg_traits->
            construct_min_vertex_2_object()(cv[cv.number_of_segments()-1]);
      }
    };

    /*! Get a Construct_min_vertex_2 functor object. */
    Construct_min_vertex_2 construct_min_vertex_2_object() const
    { return Construct_min_vertex_2(*this); }

    class Construct_max_vertex_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Construct_max_vertex_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits)
      {}

      /*!
       * Get the right endpoint of the x-monotone curve(segment).
       * \param cv The polyline.
       * \return The right endpoint.
       */
      const Point_2& operator()(const X_monotone_curve_2& cv) const
      {
        CGAL_assertion(cv.number_of_segments() > 0);

        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();

        if (seg_traits->compare_endpoints_xy_2_object()(cv[0]) == SMALLER)
          return seg_traits->
            construct_max_vertex_2_object()(cv[cv.number_of_segments()-1]);
        else
          return seg_traits->construct_max_vertex_2_object()(cv[0]);
      }
    };

    /*! Get a Construct_max_vertex_2 functor object. */
    Construct_max_vertex_2 construct_max_vertex_2_object() const
    { return Construct_max_vertex_2(*this); }

    class Is_vertical_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Is_vertical_2(const Polyline_traits_2& traits) : m_poly_traits(traits) {}

      /*!
       * Check whether the given x-monotone curve is a vertical segment.
       * \param cv The curve.
       * \return (true) if the curve is a vertical segment;(false) otherwise.
       */
      bool operator()(const X_monotone_curve_2& cv) const
      {
        // An x-monotone polyline can represent a vertical segment only if it
        // is comprised of vertical segments. If the first segment is vertical,
        // all segments are vertical in an x-monotone polyline
        return m_poly_traits.segment_traits_2()->is_vertical_2_object()(cv[0]);
      }
    };

    /*! Get an Is_vertical_2 functor object. */
    Is_vertical_2 is_vertical_2_object() const
    { return Is_vertical_2(*this); }

    class Compare_y_at_x_2 {
    private:
      // Oblivious implementation
      Comparison_result operator()(const X_monotone_segment_2& xs1,
                                   Arr_curve_end ce1,
                                   const X_monotone_segment_2& xs2,
                                   Arr_all_sides_oblivious_tag) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        const Point_2& p = (ce1 == ARR_MAX_END) ?
          seg_traits->construct_max_vertex_2_object()(xs1) :
          seg_traits->construct_min_vertex_2_object()(xs1);
        return seg_traits->compare_y_at_x_2_object()(p, xs2);
      }

      // Boundary implementation
      Comparison_result operator()(const X_monotone_segment_2& xs1,
                                   Arr_curve_end ce1,
                                   const X_monotone_segment_2& xs2,
                                   Arr_not_all_sides_oblivious_tag) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        typename Segment_traits_2::Parameter_space_in_x_2 ps_x =
          seg_traits->parameter_space_in_x_2_object();
        typename Segment_traits_2::Parameter_space_in_y_2 ps_y =
          seg_traits->parameter_space_in_y_2_object();
        typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
          seg_traits->construct_min_vertex_2_object();
        typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
          seg_traits->construct_max_vertex_2_object();

        const Arr_parameter_space ps_x1 = ps_x(xs1, ce1);
        const Arr_parameter_space ps_y1 = ps_y(xs1, ce1);

        CGAL_assertion(((ce1 == ARR_MAX_END) && (ps_x1 != ARR_LEFT_BOUNDARY)) ||
                       ((ce1 == ARR_MIN_END) && (ps_x1 != ARR_RIGHT_BOUNDARY)));

        if (ps_x1 == ARR_INTERIOR) {
          if (ps_y1 == ARR_TOP_BOUNDARY) {
            typename Segment_traits_2::Equal_2 equal =
              seg_traits->equal_2_object();
            const Point_2& p = (ce1 == ARR_MAX_END) ?
              max_vertex(xs1) : min_vertex(xs1);
            if (equal(p, max_vertex(xs2))) return EQUAL;
            if (equal(p, min_vertex(xs2))) return EQUAL;
            return LARGER;
          }
          if (ps_y1 == ARR_BOTTOM_BOUNDARY) {
            typename Segment_traits_2::Equal_2 equal =
              seg_traits->equal_2_object();
            const Point_2& p = (ce1 == ARR_MAX_END) ?
              max_vertex(xs1) : min_vertex(xs1);
            if (equal(p, max_vertex(xs2))) return EQUAL;
            if (equal(p, min_vertex(xs2))) return EQUAL;
            return SMALLER;
          }
          // ps_y1 == ARR_INTERIOR
          const Point_2& p = (ce1 == ARR_MAX_END) ?
            max_vertex(xs1) : min_vertex(xs1);
          return seg_traits->compare_y_at_x_2_object()(p, xs2);
        }
        // ps_x1 == ARR_RIGHT_BOUNDARY || ARR_LEFT_BOUNDARY
        const Point_2& p1 = (ce1 == ARR_MAX_END) ?
          max_vertex(xs1) : min_vertex(xs1);
        const Point_2& p2 = (ce1 == ARR_MAX_END) ?
          max_vertex(xs2) : min_vertex(xs2);
        return seg_traits->compare_y_on_boundary_2_object()(p1, p2);
      }

    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>       Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Compare_y_at_x_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits) {}

      /*!
       * Return the location of the given point with respect to the input curve.
       * \param p The point.
       * \param xcv The polyline curve.
       * \pre p is in the x-range of cv.
       * \return SMALLER if y(p) < cv(x(p)), i.e. the point is below the curve;
       *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
       *         EQUAL if p lies on the curve.
       */
      Comparison_result operator()(const Point_2& p,
                                   const X_monotone_curve_2& xcv) const
      {
        // Get the index of the segment in xcv containing p.
        std::size_t i = m_poly_traits.locate(xcv, p);
        CGAL_precondition(i != INVALID_INDEX);

        // Compare the segment xcv[i] and p.
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        return seg_traits->compare_y_at_x_2_object()(p, xcv[i]);
      }

      /*!
       * Return the location of the given curve_end with respect to the input
       * curve.
       * \param xcv The polyline curve.
       * \param ce the curve-end indicator of the x-monotone segment xl:
       *            ARR_MIN_END - the minimal end of xl or
       *            ARR_MAX_END - the maximal end of xl.
       * \param xcv The polyline curve.
       * \pre the curve-end is in the x-range of xcv.
       * \return SMALLER if if y(xs, ce) < cv(x(xs, ce)), i.e.  the curve-end
       *           is below the curve xcv;
       *         LARGER if y(xs, ce) > cv(x(xs, ce)), i.e. the curve-end is
       *           above the curve xcv;
       *         EQUAL if the curve-end lies on the curve xcv.
       */
      Comparison_result operator()(const X_monotone_segment_2& xs1,
                                   Arr_curve_end ce1,
                                   const X_monotone_segment_2& xs2) const
      { return operator()(xs1, ce1, xs2, Are_all_sides_oblivious_tag()); }
    };

    /*! Get a Compare_y_at_x_2 functor object. */
    Compare_y_at_x_2 compare_y_at_x_2_object() const
    { return Compare_y_at_x_2(*this); }

    class Compare_y_at_x_left_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>       Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Compare_y_at_x_left_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits)
      {}

      /*!
       * Compare the y value of two x-monotone curves immediately to the left
       * of their intersection point.
       * \param cv1 The first polyline curve.
       * \param cv2 The second polyline curve.
       * \param p The intersection point.
       * \pre The point p lies on both curves, and both of them must be also be
       *      defined(lexicographically) to its left.
       * \return The relative position of cv1 with respect to cv2 immdiately to
       *         the left of p: SMALLER, LARGER or EQUAL.
       */
      Comparison_result operator()(const X_monotone_curve_2& cv1,
                                   const X_monotone_curve_2& cv2,
                                   const Point_2& p) const
      {
        // Get the indices of the segments in cv1 and cv2 containing p and
        // defined to its left.
        std::size_t i1 = m_poly_traits.locate_side(cv1, p, false);
        std::size_t i2 = m_poly_traits.locate_side(cv2, p, false);

        CGAL_precondition(i1 != INVALID_INDEX);
        CGAL_precondition(i2 != INVALID_INDEX);

        // Compare cv1[i1] and cv2[i2] at p's left.
        return m_poly_traits.segment_traits_2()->
          compare_y_at_x_left_2_object()(cv1[i1], cv2[i2], p);
      }
    };

    /*! Get a Compare_y_at_x_left_2 functor object. */
    Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
    { return Compare_y_at_x_left_2(*this); }

    class Compare_y_at_x_right_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>       Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Compare_y_at_x_right_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits)
      {}

      /*!
       * Compare the y value of two x-monotone curves immediately to the right
       * of their intersection point.
       * \param cv1 The first curve.
       * \param cv2 The second curve.
       * \param p The intersection point.
       * \pre The point p lies on both curves, and both of them must be also be
       *      defined(lexicographically) to its right.
       * \return The relative position of cv1 with respect to cv2 immdiately to
       *         the right of p: SMALLER, LARGER or EQUAL.
       */
      Comparison_result operator()(const X_monotone_curve_2& cv1,
                                   const X_monotone_curve_2& cv2,
                                   const Point_2& p) const
      {
        // Get the indices of the segments in cv1 and cv2 containing p and
        // defined to its right.
        std::size_t i1=m_poly_traits.locate_side(cv1, p, true);
        std::size_t i2=m_poly_traits.locate_side(cv2, p, true);

        CGAL_precondition(i1 != INVALID_INDEX);
        CGAL_precondition(i2 != INVALID_INDEX);

        // Compare cv1[i1] and cv2[i2] at p's right.
        return m_poly_traits.segment_traits_2()->
          compare_y_at_x_right_2_object()(cv1[i1], cv2[i2], p);
      }
    };

    /*! Get a Compare_y_at_x_right_2 functor object. */
    Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
    { return Compare_y_at_x_right_2(*this); }

    class Equal_2 {
    protected:

      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:

      /*! Constructor. */
      Equal_2(const Polyline_traits_2& poly_tr) : m_poly_traits(poly_tr) {}

      /*!
       * Check if the two points are the same.
       * \param p1 The first point.
       * \param p2 The second point.
       * \return (true) if the two point are the same;(false) otherwise.
       */
      bool operator()(const Point_2& p1, const Point_2& p2) const
      { return m_poly_traits.segment_traits_2()->equal_2_object()(p1, p2); }

      /*!
       * Check if the two x-monotone curves are the same(have the same graph).
       * \param cv1 The first curve.
       * \param cv2 The second curve.
       * \return(true) if the two curves are the same;(false) otherwise.
       */
      bool operator()(const X_monotone_curve_2 & cv1,
                      const X_monotone_curve_2 & cv2) const
      {
        // Check the pairwise equality of the contained segments.
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        typename Segment_traits_2::Equal_2 equal =
          seg_traits->equal_2_object();
        typename Segment_traits_2::Compare_x_2 compare_x =
          seg_traits->compare_x_2_object();
        typename Segment_traits_2::Compare_y_at_x_2 compare_y_at_x =
          seg_traits->compare_y_at_x_2_object();
        typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
          seg_traits->construct_max_vertex_2_object();
        typename Segment_traits_2::Compare_endpoints_xy_2 comp_endpt =
          seg_traits->compare_endpoints_xy_2_object();
        Is_vertical_2 is_vertical = m_poly_traits.is_vertical_2_object();
        Construct_min_vertex_2 xpoly_min_v =
          m_poly_traits.construct_min_vertex_2_object();
        Construct_max_vertex_2 xpoly_max_v =
          m_poly_traits.construct_max_vertex_2_object();

        // The first and last points of the segments should be equal.
        bool res = equal(xpoly_min_v(cv1), xpoly_min_v(cv2));
        if (!res) return false;
        res = equal(xpoly_max_v(cv1), xpoly_max_v(cv2));
        if (!res) return false;

        // If the first and last points are equal and the curves are vertical,
        // it means that it is equal.
        bool ver1 = is_vertical(cv1);
        bool ver2 = is_vertical(cv2);
        // both curves are vertical and therefore equal.
        if (ver1 && ver2) return true;
        // one is vertical and the other is not - hence not equal.
        if (ver1 || ver2) return false;

        // If we arrived here it means that the first and last point of the
        // curve are equal.
        std::size_t i = 0;
        std::size_t j = 0;
        std::size_t n1 = cv1.number_of_segments();
        std::size_t n2 = cv2.number_of_segments();
        Comparison_result is_cv1_left_to_right = comp_endpt(cv1[0]);
        Comparison_result is_cv2_left_to_right = comp_endpt(cv2[0]);

        while ((i < (n1-1)) || (j < (n2-1))) {
          Point_2 point1, point2;
          std::size_t cv1_seg_ind, cv2_seg_ind;
          if (SMALLER == is_cv1_left_to_right) {
            cv1_seg_ind = i;
            point1 = max_vertex(cv1[cv1_seg_ind]);
          }
          else {
            cv1_seg_ind = n1 - 1 - i;
            point1 = max_vertex(cv1[cv1_seg_ind]);
          }
          if (SMALLER == is_cv2_left_to_right) {
            cv2_seg_ind = j;
            point2 = max_vertex(cv2[cv2_seg_ind]);
          }
          else {
            cv2_seg_ind = n2 - 1 - j;
            point2 = max_vertex(cv2[cv2_seg_ind]);
          }

          bool res = equal(point1, point2);
          // Easy case - the two points are equal
          if (res) {
            ++i;
            ++j;
          }
          else {
            Comparison_result res_x = compare_x(point1,point2);
            // Check if the different point is a collinear point situated on
            // the line between its two neighbors.
            if (SMALLER == res_x) {
              Comparison_result res_y_at_x =
                compare_y_at_x(point1, cv2[cv2_seg_ind]);
              if (EQUAL == res_y_at_x) ++i;
              else return false;
            }
            else if (LARGER == res_x) {
              Comparison_result res_y_at_x =
                compare_y_at_x(point2,cv1[cv1_seg_ind]);
              if (EQUAL == res_y_at_x) ++j;
              else return false;
            }
            else {
              return false;
            }
          }
        }
        return true;
      }
    };

    /*! Get an Equal_2 functor object. */
    Equal_2 equal_2_object() const
    { return Equal_2(*this); }

    class Compare_endpoints_xy_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2> Polyline_traits_2;
      /*! The traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Compare_endpoints_xy_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits) {}

      /*!
       * Compare the endpoints of an \(x\)-monotone curve lexicographically.
       * (assuming the curve has a designated source and target points).
       * \param cv The curve.
       * \return SMALLER if the curve is oriented left-to-right;
       *         LARGER if the curve is oriented right-to-left.
       */
      Comparison_result operator()(const X_monotone_curve_2& xcv) const
      {
        return (m_poly_traits.segment_traits_2()->
                compare_endpoints_xy_2_object()(xcv[0]) == SMALLER) ?
          (SMALLER) : (LARGER);
      }
    };

    Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
    { return Compare_endpoints_xy_2(*this); }

    class Construct_opposite_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2> Polyline_traits_2;
      /*! The traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor */
      Construct_opposite_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits) {}

      /*!
       * Construct the reversed \(x\)-monotone polyline of the input.
       * Note that the functor constructs the opposites of _all_ segments
       * constituting xcv.
       * \param xcv the \(x\)-monotone polyline to be reveres
       * \pre xcv contains at least one segment
       * \return An \(x\)-monotone polyline with the same graph as the input xcv
       *         only with a reverse orientation.
       */
      X_monotone_curve_2 operator()(const X_monotone_curve_2& xcv) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        typename Segment_traits_2::Construct_opposite_2 const_op =
          seg_traits->construct_opposite_2_object();
        std::vector<X_monotone_segment_2> rev_segs;
        typename X_monotone_curve_2::Segment_const_iterator it;
        for (it = xcv.begin_segments(); it != xcv.end_segments(); ++it)
          rev_segs.push_back(const_op(*it));
        return X_monotone_curve_2(rev_segs.rbegin(),rev_segs.rend());
      }
    };

    Construct_opposite_2 construct_opposite_2_object() const
    { return Construct_opposite_2(*this); }

    ///@}

    /// \name Construction functors(based on the segment traits).
    //@{

    class Make_x_monotone_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /*! The traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Make_x_monotone_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits) {}

      /*!
       * Cut the given curve into x-monotone sub-curves and insert them into the
       * given output iterator.
       *
       * \pre if `cv` is not empty then it must be continuous and well-oriented.
       * \param cv The curve.
       * \param oi The output iterator, whose value-type is Object. The output
       *           object is a wrapper of a X_monotone_curve_2.
       * \return The past-the-end iterator.
       */
      template<class OutputIterator>
      OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const
      {
        std::cout << std::endl;

        typedef typename Curve_2::Segment_const_iterator const_seg_iterator;

        // If the polyline is empty, return.
        if (cv.number_of_segments() == 0) return oi;

        Construct_x_monotone_curve_2 ctr_x_curve =
          m_poly_traits.construct_x_monotone_curve_2_object();

        typename Segment_traits_2::Make_x_monotone_2 make_seg_x_monotone =
          m_poly_traits.segment_traits_2()->make_x_monotone_2_object();

        typename Segment_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
          m_poly_traits.segment_traits_2()->compare_endpoints_xy_2_object();

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
        typename Segment_traits_2::Construct_opposite_2 ctr_seg_opposite =
          m_poly_traits.segment_traits_2()->construct_opposite_2_object();
#endif

        // Convert the input polyline to a sequence of CGAL objects, such
        // that each Object wraps an x-monotone segment.
        std::vector<Object> x_seg_objects;
        const_seg_iterator it_segs;
        for (it_segs = cv.begin_segments(); it_segs != cv.end_segments();
             ++it_segs)
          make_seg_x_monotone(*it_segs, std::back_inserter(x_seg_objects));
        typename std::vector<Object>::iterator it = x_seg_objects.begin();
        X_monotone_segment_2 x_seg;
#if defined (CGAL_NO_ASSERTIONS)
        CGAL::assign(x_seg, *it);
#else
        bool assign_res = CGAL::assign(x_seg, *it);
        CGAL_assertion(assign_res);
#endif

        // If the polyline consists of a single x-monotone segment, return.
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
           typename Segment_traits_2::Construct_min_vertex_2 min_seg_v =
             m_poly_traits.segment_traits_2()->construct_min_vertex_2_object();
           typename Segment_traits_2::Construct_max_vertex_2 max_seg_v =
             m_poly_traits.segment_traits_2()->construct_max_vertex_2_object();
           typename Segment_traits_2::Equal_2 equal =
             m_poly_traits.segment_traits_2()->equal_2_object();
           Point_2 last_target = (cmp_seg_endpts(x_seg) == SMALLER) ?
             max_seg_v(x_seg) : min_seg_v(x_seg);
           Point_2 next_src;
           );

        // The polyline consists of at least 2 x-monotone segments:
        Push_back_2 push_back = m_poly_traits.push_back_2_object();
        typename Segment_traits_2::Is_vertical_2 is_seg_vertical =
          m_poly_traits.segment_traits_2()->is_vertical_2_object();

        bool is_start_vertical = is_seg_vertical(x_seg);
        Comparison_result start_dir = cmp_seg_endpts(x_seg);

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
        Push_front_2 push_front = m_poly_traits.push_front_2_object();
        if (cmp_seg_endpts(x_seg) == LARGER)
          x_seg = ctr_seg_opposite(x_seg);
#endif
        X_monotone_curve_2 x_polyline = ctr_x_curve(x_seg);

        for (++it; it != x_seg_objects.end(); ++it) {
          X_monotone_segment_2 x_seg;
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
            *oi++ = make_object(x_polyline);
            is_start_vertical = is_seg_vertical(x_seg);
            start_dir = cmp_seg_endpts(x_seg);
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
            if (cmp_seg_endpts(x_seg) == LARGER)
              x_seg = ctr_seg_opposite(x_seg);
#endif
            x_polyline = ctr_x_curve(x_seg);
          }
          else {
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
            if (cmp_seg_endpts(x_seg) == LARGER) {
              x_seg = ctr_seg_opposite(x_seg);
              push_front(x_polyline, x_seg);
            }
            else
              push_back(x_polyline, x_seg);
#else
            push_back(x_polyline, x_seg);
#endif
          }
        }
        if (x_polyline.number_of_segments() != 0)
          *oi++ = make_object(x_polyline);
        x_seg_objects.clear();
        return oi;
      }
    };

    /*! Get a Make_x_monotone_2 functor object. */
    Make_x_monotone_2 make_x_monotone_2_object() const
    { return Make_x_monotone_2(*this); }

    /* Functor to augment a polyline by either adding a vertex or a segment
     * at the back.
     * TODO: Test all the operator()'s. (Don't forget vertical cases!)
     */
    class Push_back_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /*! The traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Push_back_2(const Polyline_traits_2& traits) : m_poly_traits(traits) {}

      /* Append a point `p` to an existing polyline `cv` at the back. */
      void operator()(Curve_2& cv, const Point_2& p) const
      {
        typedef typename Curve_2::Segments_size_type size_type;
        size_type num_seg = cv.number_of_segments();
        CGAL_precondition(num_seg > 0);
        int last_seg = num_seg-1;

        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        typename Segment_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
          seg_traits->compare_endpoints_xy_2_object();

        /*
         * Since it is desired to maintain `cv` well-oriented, we have
         * to append the segment [cv[last_seg].target(),p]. The
         * following test determines which end of the last segment is
         * the target, i.e. the actual end of `cv`.
         */
        if (cmp_seg_endpts(cv[last_seg]) == SMALLER) {
          typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
            seg_traits->construct_max_vertex_2_object();
          cv.push_back(Segment_2(get_max_v(cv[last_seg]), p));
        }
        else {
          typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
            seg_traits->construct_min_vertex_2_object();
          cv.push_back(Segment_2(get_min_v(cv[last_seg]), p));
        }
      }


      /* Append a segment `seg` to an existing polyline `cv`. If `cv` is
         empty, `seg` will be its first segment. */
      void operator()(Curve_2& cv, const Segment_2& seg) const
      { cv.push_back(seg); }

      /* Append a point `p` to an existing polyline `xcv` at the back. */
      void operator()(X_monotone_curve_2& xcv, const Point_2& p) const
      {
        typedef typename X_monotone_curve_2::Segments_size_type size_type;
        size_type num_seg = xcv.number_of_segments();
        CGAL_precondition(num_seg > 0);

        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        CGAL_precondition_code
          (
           typename Segment_traits_2::Compare_x_2 comp_x =
             seg_traits->compare_x_2_object();
           typename Segment_traits_2::Compare_xy_2 comp_xy =
             seg_traits->compare_xy_2_object();
           Is_vertical_2 is_vertical = m_poly_traits.is_vertical_2_object();
           );

        if (seg_traits->compare_endpoints_xy_2_object()(xcv[0]) == SMALLER) {
          // xcv is oriented left-to-right
          typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
            seg_traits->construct_max_vertex_2_object();
          CGAL_precondition
            (
             (!is_vertical(xcv) &&
              (comp_x(get_max_v(xcv[num_seg-1]), p) == SMALLER)) ||
             (is_vertical(xcv) &&
              (comp_x(get_max_v(xcv[num_seg-1]), p) == EQUAL) &&
              (comp_xy(get_max_v(xcv[num_seg-1]), p) == SMALLER))
             );
          xcv.push_back(X_monotone_segment_2(get_max_v(xcv[num_seg-1]), p));
        }
        else {
          // xcv is oriented right-to-left
          typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
            seg_traits->construct_min_vertex_2_object();
          CGAL_precondition
            (
             (!is_vertical(xcv) &&
              (comp_x(get_min_v(xcv[num_seg-1]), p) == LARGER)) ||
             (is_vertical(xcv) &&
              (comp_x(get_min_v(xcv[num_seg-1]), p) == EQUAL) &&
              (comp_xy(get_min_v(xcv[num_seg-1]), p) == LARGER))
             );
          xcv.push_back(X_monotone_segment_2(get_min_v(xcv[num_seg-1]), p));
        }
      }

      /* Append a segment `seg` to an existing polyline `xcv` at the back. */
      void operator()(X_monotone_curve_2& xcv,
                      const X_monotone_segment_2& seg) const
      {
        CGAL_precondition_code
          (
           typedef typename X_monotone_curve_2::Segments_size_type size_type;
           size_type num_seg = xcv.number_of_segments();
           const Segment_traits_2* seg_traits =
             m_poly_traits.segment_traits_2();
           typename Segment_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
             seg_traits->compare_endpoints_xy_2_object();
           Comparison_result dir = cmp_seg_endpts(seg);
           typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
             seg_traits->construct_max_vertex_2_object();
           typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
             seg_traits->construct_min_vertex_2_object();
           typename Segment_traits_2::Equal_2 equal =
             seg_traits->equal_2_object();
           typename Segment_traits_2::Is_vertical_2 is_vertical =
             seg_traits->is_vertical_2_object();
           );

        CGAL_precondition_msg((num_seg == 0) ||
                              ((is_vertical(xcv[0]) && is_vertical(seg)) ||
                               (!is_vertical(xcv[0]) && !is_vertical(seg))),
                              "xcv is vertical and seg is not or vice versa!");

        CGAL_precondition_msg((num_seg == 0) || (cmp_seg_endpts(xcv[0]) == dir),
                              "xcv and seg do not have the same orientation!");

        CGAL_precondition_msg((num_seg == 0) ||
                              !equal(get_min_v(seg), get_max_v(seg)),
                              "Seg degenerates to a point!");

        CGAL_precondition_msg((num_seg == 0) ||
                              (((dir != SMALLER) ||
                                equal(get_max_v(xcv[num_seg-1]),
                                      get_min_v(seg)))),
                              "Seg does not extend to the right!");

        CGAL_precondition_msg((num_seg == 0) ||
                              (((dir != LARGER) ||
                                equal(get_min_v(xcv[num_seg-1]),
                                      get_max_v(seg)))),
                              "Seg does not extend to the left!");

        xcv.push_back(seg);
      }
    };

    /*! Get a Push_back_2 functor object. */
    Push_back_2 push_back_2_object() const
    { return Push_back_2(*this); }

    /* Functor to augment a polyline by either adding a vertex or a segment
     * at the front.
     * TODO: Test all the operator()'s. (Don't forget vertical cases!)
     */
    class Push_front_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /*! The traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Push_front_2(const Polyline_traits_2& traits) : m_poly_traits(traits) {}

      /* Append a point `p` to an existing polyline `cv` at the front. */
      void operator()(Curve_2& cv, const Point_2& p) const
      {
        CGAL_precondition_code
          (
           typedef typename Curve_2::Segments_size_type size_type;
           size_type num_seg = cv.number_of_segments();
           );
        CGAL_precondition(num_seg > 0);

        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        typename Segment_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
          seg_traits->compare_endpoints_xy_2_object();

        if (cmp_seg_endpts(cv[0]) == SMALLER) {
          typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
            seg_traits->construct_min_vertex_2_object();
          cv.push_front(Segment_2(p, get_min_v(cv[0])));
        }
        else {
          typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
            seg_traits->construct_max_vertex_2_object();
          cv.push_front(Segment_2(p, get_max_v(cv[0])));
        }
      }


      /* Append a segment `seg` to an existing polyline `cv` at the front. */
      void operator()(Curve_2& cv, const Segment_2& seg) const
      { cv.push_front(seg); }

      /* Append a point `p` to an existing polyline `xcv` at the front. */
      void operator()(const X_monotone_curve_2& xcv, Point_2& p) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        CGAL_precondition_code
          (
           typedef typename X_monotone_curve_2::Segments_size_type size_type;
           size_type num_seg = xcv.number_of_segments();
           typename Segment_traits_2::Compare_x_2 comp_x =
             seg_traits->compare_x_2_object();
           typename Segment_traits_2::Compare_xy_2 comp_xy =
             seg_traits->compare_xy_2_object();
           Is_vertical_2 is_vertical = m_poly_traits.is_vertical_2_object();
           );
        CGAL_precondition(num_seg > 0);

        if (seg_traits->compare_endpoints_xy_2_object()(xcv[0]) == SMALLER) {
          // xcv is oriented left-to-right
          typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
            seg_traits->construct_max_vertex_2_object();
          CGAL_precondition
            (
             (!is_vertical(xcv) &&
              (comp_x(get_max_v(xcv[0]), p) == LARGER)) ||
             (is_vertical(xcv) &&
              (comp_x(get_max_v(xcv[0]), p) == EQUAL) &&
              (comp_xy(get_max_v(xcv[0]), p) == LARGER))
             );
          xcv.push_front(X_monotone_segment_2(p, get_max_v(xcv[0])));
        }
        else {
          // xcv is oriented right-to-left
          typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
            seg_traits->construct_min_vertex_2_object();
          CGAL_precondition
            (
             (!is_vertical(xcv) &&
              (comp_x(get_min_v(xcv[0]), p) == SMALLER)) ||
             (is_vertical(xcv) &&
              (comp_x(get_min_v(xcv[0]), p) == EQUAL) &&
              (comp_xy(get_min_v(xcv[0]), p) == SMALLER))
             );
          xcv.push_front(X_monotone_segment_2(p, get_min_v(xcv[0])));
        }
      }

      /* Append a segment `seg` to an existing polyline `xcv` at the front. */
      void operator()(X_monotone_curve_2& xcv,
                      const X_monotone_segment_2& seg) const
      {
        CGAL_precondition_code
          (
           typedef typename X_monotone_curve_2::Segments_size_type size_type;
           size_type num_seg = xcv.number_of_segments();
           const Segment_traits_2* seg_traits =
             m_poly_traits.segment_traits_2();
           typename Segment_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
             seg_traits->compare_endpoints_xy_2_object();
           Comparison_result dir = cmp_seg_endpts(seg);
           typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
             seg_traits->construct_max_vertex_2_object();
           typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
             seg_traits->construct_min_vertex_2_object();
           typename Segment_traits_2::Compare_xy_2 comp_xy =
             seg_traits->compare_xy_2_object();
           typename Segment_traits_2::Equal_2 equal =
             seg_traits->equal_2_object();
           typename Segment_traits_2::Is_vertical_2 is_vertical =
             seg_traits->is_vertical_2_object();
           );

        CGAL_precondition_msg((num_seg == 0) ||
                              ((is_vertical(xcv[0]) && is_vertical(seg)) ||
                               (!is_vertical(xcv[0]) && !is_vertical(seg))),
                              "xcv is vertical and seg is not or vice versa!");

        CGAL_precondition_msg((num_seg == 0) || (cmp_seg_endpts(xcv[0]) == dir),
                              "xcv and seg do not have the same orientation!");

        CGAL_precondition_msg((num_seg == 0) ||
                              !equal(get_min_v(seg), get_max_v(seg)),
                              "Seg degenerates to a point!");

        CGAL_precondition_msg((num_seg == 0) ||
                              (((dir != SMALLER) ||
                                equal(get_min_v(xcv[0]), get_max_v(seg)))),
                              "Seg does not extend to the left!");

        CGAL_precondition_msg((num_seg == 0) ||
                              (((dir != LARGER) ||
                                equal(get_max_v(xcv[0]), get_min_v(seg)))),
                              "Seg does not extend to the right!");

        xcv.push_front(seg);
      }
    };

    /*! Get a Push_front_2 functor object. */
    Push_front_2 push_front_2_object() const
    { return Push_front_2(*this); }

    class Split_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>       Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Split_2(const Polyline_traits_2& traits) : m_poly_traits(traits) {}

    public:
      /*!
       * Split a given x-monotone curve at a given point into two sub-curves.
       * \param cv The curve to split
       * \param p The split point.
       * \param c1 Output: The left resulting subcurve(p is its right endpoint).
       * \param c2 Output: The right resulting subcurve(p is its left endpoint).
       * \pre p lies on cv but is not one of its end-points.
       */
      void operator()(const X_monotone_curve_2& xcv, const Point_2& p,
                      X_monotone_curve_2& xcv1, X_monotone_curve_2& xcv2) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
          seg_traits->construct_min_vertex_2_object();
        typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
          seg_traits->construct_max_vertex_2_object();
        typename Segment_traits_2::Equal_2 equal =
          seg_traits->equal_2_object();
        typename Segment_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
          seg_traits->compare_endpoints_xy_2_object();

        // Make sure the split point is not one of the curve endpoints.
        CGAL_precondition((!equal(m_poly_traits.
                                  construct_min_vertex_2_object()(xcv), p)));
        CGAL_precondition((!equal(m_poly_traits.
                                  construct_max_vertex_2_object()(xcv), p)));

        CGAL_precondition_msg(xcv.number_of_segments() > 0,
                              "Cannot split a polyline of length zero.");

        Comparison_result dir = cmp_seg_endpts(xcv[0]);

        // Locate the segment on the polyline xcv that contains p.
        std::size_t i = m_poly_traits.locate(xcv, p);

        CGAL_precondition(i != INVALID_INDEX);

        // Clear the output curves.
        xcv1.clear();
        xcv2.clear();

        // Push all segments labeled(0, 1, ... , i-1) into xcv1.
        for (std::size_t j = 0; j < i; ++j)
          xcv1.push_back(xcv[j]);

        if (dir == SMALLER){
          // Check whether the split point is xcv[i]'s source or target.
          if (equal(max_vertex(xcv[i]), p)) {
            // The entire i'th segment belongs to xcv1:
            xcv1.push_back(xcv[i]);
          } else if (equal(min_vertex(xcv[i]), p)) {
            // The entire i'th segments belongs to xcv2:
            xcv2.push_back(xcv[i]);
          } else {
            // The i'th segment should be split: The left part(seg1)
            // goes to xcv1, and the right part(seg2) goes to xcv2.
            X_monotone_segment_2 seg1, seg2;
            m_poly_traits.segment_traits_2()->
              split_2_object()(xcv[i], p, seg1, seg2);

            xcv1.push_back(seg1);
            xcv2.push_back(seg2);
          }
        }
        else{
          if (equal(min_vertex(xcv[i]), p)) {
            xcv1.push_back(xcv[i]);
          } else if (equal(max_vertex(xcv[i]), p)) {
            xcv2.push_back(xcv[i]);
          } else {
            X_monotone_segment_2 seg1, seg2;
            m_poly_traits.segment_traits_2()->
              split_2_object()(xcv[i], p, seg1, seg2);

            if (cmp_seg_endpts(seg2) == LARGER){
              xcv1.push_back(seg2);
            }
            else{
              // seg2 has to be reversed
              seg2 = m_poly_traits.segment_traits_2()->
                construct_opposite_2_object()(seg2);
              xcv1.push_back(seg2);
            }

            if (cmp_seg_endpts(seg1) == LARGER){
              xcv2.push_back(seg1);
            } else {
              // seg2 has to be reversed
              seg1 = m_poly_traits.segment_traits_2()->
                construct_opposite_2_object()(seg1);
              xcv1.push_back(seg1);
            }
          }
        }

        // Push all segments labeled(i+1, i+2, ... , n-1) into xcv1.
        std::size_t n = xcv.number_of_segments();

        for (std::size_t j = i+1; j < n; ++j)
          xcv2.push_back(xcv[j]);

        if (dir != SMALLER){
          X_monotone_curve_2 tmp;
          tmp = xcv1;
          xcv1 = xcv2;
          xcv2 = tmp;
        }
      }
    };

    /*! Get a Split_2 functor object. */
    Split_2 split_2_object() const
    { return Split_2(*this); }

    class Intersect_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>       Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Intersect_2(const Polyline_traits_2& traits) : m_poly_traits(traits) {}

      /*!
       * Find the intersections of the two given curves and insert them into the
       * given output iterator. As two segments may itersect only once, only a
       * single intersection will be contained in the iterator.
       * Note: If the intersection yields an overlap then it will be oriented
       *       from left-to-right.
       * \param cv1 The first curve.
       * \param cv2 The second curve.
       * \param oi The output iterator.
       * \return The past-the-end iterator.
       */
      template <typename OutputIterator>
      OutputIterator operator()(const X_monotone_curve_2& cv1,
                                const X_monotone_curve_2& cv2,
                                OutputIterator oi) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        Compare_y_at_x_2 cmp_y_at_x = m_poly_traits.compare_y_at_x_2_object();
        typename Segment_traits_2::Equal_2 equal = seg_traits->equal_2_object();
        typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
          seg_traits->construct_min_vertex_2_object();
        typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
          seg_traits->construct_max_vertex_2_object();
        typename Segment_traits_2::Intersect_2 intersect =
          seg_traits->intersect_2_object();
        typename Segment_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
          seg_traits->compare_endpoints_xy_2_object();

        Comparison_result dir1 = cmp_seg_endpts(cv1[0]);
        Comparison_result dir2 = cmp_seg_endpts(cv2[0]);

        const std::size_t n1 = cv1.number_of_segments();
        const std::size_t n2 = cv2.number_of_segments();

        std::size_t i1 = (dir1 == SMALLER) ? 0 : n1-1;
        std::size_t i2 = (dir2 == SMALLER) ? 0 : n2-1;

        X_monotone_curve_2 ocv;           // Used to represent overlaps.

        Compare_xy_2 compare_xy = m_poly_traits.compare_xy_2_object();
        Comparison_result left_res =
          compare_xy(cv1[i1], ARR_MIN_END, cv2[i2], ARR_MIN_END);

        if (left_res == SMALLER) {
          // cv1's left endpoint is to the left of cv2's left endpoint:
          // Locate the index i1 of the segment in cv1 which contains cv2's
          // left endpoint.
          i1 = m_poly_traits.locate_impl(cv1, cv2[i2], ARR_MIN_END,
                                         Are_all_sides_oblivious_tag());
          if (i1 == INVALID_INDEX) return oi;

          if (equal(max_vertex(cv1[i1]), min_vertex(cv2[i2]))) {
            if (((dir1 == SMALLER) && (i1 == n1-1)) ||
                ((dir1 == LARGER) && (i1 == 0))){
              // cv1's right endpoint equals cv2's left endpoint
              // Thus we can return this single(!) intersection point
              std::pair<Point_2, Multiplicity>  p(max_vertex(cv1[i1]), 0);
              *oi++ = make_object(p);
              return oi;
            }
            dir1 == SMALLER ? ++i1 : (i1 != 0) ? --i1 : std::size_t(INVALID_INDEX);
            left_res = EQUAL;
          }
        }
        else if (left_res == LARGER) {
          // cv1's left endpoint is to the right of cv2's left endpoint:
          // Locate the index i2 of the segment in cv2 which contains cv1's
          // left endpoint.
          i2 = m_poly_traits.locate_impl(cv2, cv1[i1], ARR_MIN_END,
                                         Are_all_sides_oblivious_tag());
          if (i2 == INVALID_INDEX) return oi;

          if (equal(max_vertex(cv2[i2]), min_vertex(cv1[i1]))) {
            if (((dir2 == SMALLER) && (i2 == n2-1)) ||
                ((dir2 == LARGER) && (i2 == 0))){
              // cv2's right endpoint equals cv1's left endpoint
              // Thus we can return this single(!) intersection point
              std::pair<Point_2, Multiplicity>  p(max_vertex(cv2[i2]), 0);
              *oi++ = make_object(p);
              return oi;
            }

            dir2 == SMALLER ? ++i2 : (i2 != 0) ? --i2 : std::size_t(INVALID_INDEX);
            left_res = EQUAL;
          }
        }

        // Check if the the left endpoint lies on the other polyline.
        bool left_coincides = (left_res == EQUAL);
        bool left_overlap = false;

        if (left_res == SMALLER)
          left_coincides = (cmp_y_at_x(cv2[i2], ARR_MIN_END, cv1[i1]) == EQUAL);
        else if (left_res == LARGER)
          left_coincides = (cmp_y_at_x(cv1[i1], ARR_MIN_END, cv2[i2]) == EQUAL);

        // The main loop: Go simultaneously over both polylines.
        Comparison_result right_res = left_res;
        bool right_coincides = left_coincides;
        bool right_overlap = false;

        while (((dir1==SMALLER) && (dir2 == SMALLER) && (i1 < n1) &&(i2 < n2))||
               ((dir1!=SMALLER) && (dir2 == SMALLER) &&
                (i1 != INVALID_INDEX) && (i2 < n2)) ||
               ((dir1==SMALLER) && (dir2 != SMALLER) && (i1 < n1) &&
                (i2 != INVALID_INDEX)) ||
               ((dir1!=SMALLER) && (dir2 != SMALLER) &&
                (i1 != INVALID_INDEX) && (i2 != INVALID_INDEX)))
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
              // Non of the endpoints of the current segment of one polyline
              // coincides with the curent segment of the other polyline:
              // Output the intersection if exists.
              oi = intersect(cv1[i1], cv2[i2], oi);
            }
            else if (right_coincides && left_coincides) {
              // An overlap exists between the current segments of the
              // polylines: Output the overlapping segment.
              right_overlap = true;
              if (left_res == SMALLER) {
                if (right_res == SMALLER) {
                  X_monotone_segment_2 seg(min_vertex(cv2[i2]),
                                           max_vertex(cv1[i1]));
                  ocv.push_back(seg);
                }
                else {
                  X_monotone_segment_2 seg(min_vertex(cv2[i2]),
                                           max_vertex(cv2[i2]));
                  ocv.push_back(seg);
                }
              }
              else {
                if (right_res == SMALLER) {
                  X_monotone_segment_2 seg(min_vertex(cv1[i1]),
                                           max_vertex(cv1[i1]));
                  ocv.push_back(seg);
                }
                else {
                  X_monotone_segment_2 seg(min_vertex(cv1[i1]),
                                           max_vertex(cv2[i2]));
                  ocv.push_back(seg);
                }
              }
            }
            else if (left_coincides && !right_coincides) {
              // The left point of the current segment of one polyline
              // coincides with the current segment of the other polyline.
              if (left_overlap) {
                // An overlap occured at the previous iteration:
                // Output the overlapping polyline.
                CGAL_assertion(ocv.number_of_segments() > 0);
                *oi++ = make_object(ocv);
                ocv.clear();
              }
              else {
                // The left point of the current segment of one
                // polyline coincides with the current segment of the
                // other polyline, and no overlap occured at the
                // previous iteration: Output the intersection
                // point. The derivative of at least one of the
                // polylines is not defined at this point, so we give
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
              if (dir2 == SMALLER)
                ++i2;
              else {
                if (i2 == 0)
                  i2 = INVALID_INDEX;
                else
                  --i2;
              }
            }
            if (right_res != LARGER) {
              if (dir1 == SMALLER)
                ++i1;
              else {
                if (i1 == 0)
                  i1 = INVALID_INDEX;
                else --i1;
              }
            }
            left_res = (right_res == SMALLER) ? LARGER :
              (right_res == LARGER) ? SMALLER : EQUAL;

            left_coincides = right_coincides;
            left_overlap = right_overlap;
          } // END of while loop

        // Output the remaining overlapping polyline, if necessary.
        if (ocv.number_of_segments() > 0) {
          *oi++ = make_object(ocv);
        }
        else if (right_coincides) {
          typedef std::pair<Point_2,Multiplicity> return_point;
          return_point ip;
          if (right_res == SMALLER) {
            ip = (dir1 == SMALLER) ?
              return_point(max_vertex(cv1[i1-1]), 0) :
              (i1 != INVALID_INDEX) ?
              return_point(max_vertex(cv1[i1+1]), 0) :
              return_point(max_vertex(cv1[0]), 0);
            *oi++ = make_object(ip);
          }
          else if (right_res == LARGER) {
            ip = (dir2 == SMALLER) ?
              return_point(max_vertex(cv2[i2-1]), 0) :
              (i2 != INVALID_INDEX) ?
              return_point(max_vertex(cv2[i2+1]), 0) :
              return_point(max_vertex(cv2[0]), 0);
            *oi++ = make_object(ip);
          }
          else if (((i1 > 0) && (dir1 == SMALLER)) ||
                   ((i1 < n1) && (dir1 != SMALLER)) ||
                   ((i1 == INVALID_INDEX) && (dir1 != SMALLER)))
          {
            ip = (dir1 == SMALLER) ?
              return_point(max_vertex(cv1[i1-1]), 0) :
              (i1 != INVALID_INDEX) ?
              return_point(max_vertex(cv1[i1+1]), 0) :
              return_point(max_vertex(cv1[0]), 0);
            *oi++ = make_object(ip);
          }
          else {
            CGAL_assertion_msg((dir2 == SMALLER && i2 > 0) ||
                               (dir2 != SMALLER && i2 < n2) ||
                               (dir2 != SMALLER &&
                                (i1 == INVALID_INDEX || i2 ==INVALID_INDEX)),
                               "Wrong index for xcv2 in Intersect_2 of "
                               "polylines.");
            ip = (dir2 == SMALLER) ?
              return_point(max_vertex(cv2[i2-1]), 0) :
              (i2 != INVALID_INDEX) ?
              return_point(max_vertex(cv2[i2+1]), 0) :
              return_point(max_vertex(cv2[0]), 0);
            *oi++ = make_object(ip);
          }
        }
        return oi;
      }
    };

    /*! Get an Intersect_2 functor object. */
    Intersect_2 intersect_2_object() const
    { return Intersect_2(*this); }

    class Are_mergeable_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>       Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Are_mergeable_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits) {}

      /*!
       * Check whether it is possible to merge two given x-monotone curves.
       * \param cv1 The first curve.
       * \param cv2 The second curve.
       * \return(true) if the two curves are mergeable, that is, they share a
       * common endpoint and the same orientation;(false) otherwise.
       */
      bool operator()(const X_monotone_curve_2& cv1,
                      const X_monotone_curve_2& cv2) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        Construct_min_vertex_2 min_vertex =
          m_poly_traits.construct_min_vertex_2_object();
        Construct_max_vertex_2 max_vertex =
          m_poly_traits.construct_max_vertex_2_object();
        typename Segment_traits_2::Equal_2 equal =
          seg_traits->equal_2_object();
        typename Segment_traits_2::Is_vertical_2 is_seg_vertical =
          seg_traits->is_vertical_2_object();

        Comparison_result dir1 =
          m_poly_traits.compare_endpoints_xy_2_object()(cv1);
        Comparison_result dir2 =
          m_poly_traits.compare_endpoints_xy_2_object()(cv2);

        if (dir1 != dir2)
          return false;

        bool ver1 = is_seg_vertical(cv1[0]);
        bool ver2 = is_seg_vertical(cv2[0]);

        return (
                (
                 (// Both are directed from left-to-right
                  (dir1 == SMALLER) &&
                  ((equal(max_vertex(cv1),min_vertex(cv2))) ||
                   (equal(max_vertex(cv2),min_vertex(cv1))))
                  ) ||
                 (// Both are directed from right-to-left
                  (dir1 == LARGER) &&
                  ((equal(min_vertex(cv1),max_vertex(cv2))) ||
                   (equal(max_vertex(cv1),min_vertex(cv2))))
                  )
                 ) &&
                (// Either both should be vertical or both should
                 // be NOT vertical.
                 (ver1 && ver2) || (!ver1 && !ver2)));
      }
    };

    /*! Get an Are_mergeable_2 functor object. */
    Are_mergeable_2 are_mergeable_2_object() const
    { return Are_mergeable_2(*this); }

    /*! \class Merge_2
     * A functor that merges two x-monotone curves into one.
     */
    /* Roadmap: Allow merging of overlapping polylines. This means also
     *          changing the segment traits class.
     */
    class Merge_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Geometry_traits;
      /*! The traits (in case it has state) */
      const Geometry_traits& m_poly_traits;

    public:
      /*! Constructor
       * \param traits the traits (in case it has state)
       */
      Merge_2(const Geometry_traits& traits) : m_poly_traits(traits) {}

      /*!
       * Merge two given x-monotone curves into a single curve(segment).
       * \param cv1 The first curve.
       * \param cv2 The second curve.
       * \param c Output: The merged curve.
       * \pre The two curves are mergeable.
       */
      void operator()(const X_monotone_curve_2 & cv1,
                      const X_monotone_curve_2 & cv2,
                      X_monotone_curve_2 & c) const
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
        if (
            // Either both are left-to-right and cv2 is to the right of cv1
            ((cmp_seg_endpts(cv1)==SMALLER) &&
             (equal(get_max_v(cv1),get_min_v(cv2)))) ||
            // or both are right-to-left and cv2 is to the left of cv1
            ((cmp_seg_endpts(cv1)==LARGER) &&
             (equal(get_min_v(cv1), get_max_v(cv2))))
            )
        {
          const std::size_t n1 = cv1.number_of_segments();
          const std::size_t n2 = cv2.number_of_segments();
          std::size_t       i;

          // cv2 extends cv1 to the right:
          for (i = 0; i < n1 - 1; ++i)
            c.push_back(cv1[i]);

          // Try to merge the to contiguous line segments:
          if (m_poly_traits.segment_traits_2()->
              are_mergeable_2_object()(cv1[n1 - 1], cv2[0])) {
            X_monotone_segment_2 seg;
            m_poly_traits.segment_traits_2()->
              merge_2_object()(cv1[n1 - 1], cv2[0], seg);
            c.push_back(seg);
          }
          else {
            c.push_back(cv1[n1 - 1]);
            c.push_back(cv2[0]);
          }

          for (i = 1; i < n2; ++i)
            c.push_back(cv2[i]);
        }
        else
          return this->operator()(cv2,cv1,c);
      }
    };

    /*! Get a Merge_2 functor object. */
    Merge_2 merge_2_object() const { return Merge_2(*this); }
    ///@}

    /// \name Functor definitions for the landmarks point-location strategy.
    //@{

#if 0
    // The following block assumes that the segment traits template parameter
    // is a model of the ArrangementLandmarkTraits concept; in other words, it
    // defines the nested types Approximate_number_type and Approximate_2 and
    // the member function approximate_2_object(). It cannot be used as is if
    // the segment traits does not model the ArrangementLandmarkTraits concept.
    // The functor Construct_x_monotone_curve_2 is provided regardless of the
    // segment traits.

    typedef typename Segment_traits_2::Approximate_number_type
      Approximate_number_type;
    typedef typename Segment_traits_2::Approximate_2    Approximate_2;

    /*! Get an Approximate_2 functor object. */
    Approximate_2 approximate_2_object() const
    { return segment_traits_2()->approximate_2_object(); }
#else
    // The following block defines the nested types Approximate_number_type and
    // Approximate_2 and the member function approximate_2_object() based on the
    // corresponding types and function definitions of the segment traits. If
    // the segment traits does not provide these definitions, they are defined
    // as dummies. Essentially, the polyline traits becomes a practical model of
    // the ArrangementLandmarkTraits concept only if the segment traits is a
    // model of this concept.
    //
    // The following implementation is inspired by
    // http://stackoverflow.com/a/11816999/1915421

    template <typename T>
    struct Void {
      typedef void type;
    };

    template <typename T, typename _ = void>
    struct has_approximate_2 {
      // Generic implementation
      typedef void                      Approximate_number_type;
      typedef void                      Approximate_2;
    };

    template <typename T>
    struct has_approximate_2<T, typename Void<typename T::Approximate_2>::type>
    {
      // Specialization for types holding a nested type T::Approximate_2
      typedef typename T::Approximate_number_type
        Approximate_number_type;
      typedef typename T::Approximate_2  Approximate_2;
    };

    typedef typename has_approximate_2<Segment_traits_2>::Approximate_number_type
      Approximate_number_type;
    typedef typename has_approximate_2<Segment_traits_2>::Approximate_2
      Approximate_2;

    /*! Get an Approximate_2 functor object. */
    Approximate_2 approximate_2_object_impl(boost::false_type) const
    { return segment_traits_2()->approximate_2_object(); }

    Approximate_2 approximate_2_object_impl(boost::true_type) const { }

    Approximate_2 approximate_2_object() const
    {
      typedef typename boost::is_same<void, Approximate_2>::type      Is_void;
      return approximate_2_object_impl(Is_void());
    }
#endif
    //@}

    class Construct_curve_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>       Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Construct_curve_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits) {}

      /* Returns an polyline connecting the two given endpoints. */
      Curve_2 operator()(const Point_2& p, const Point_2& q) const
      {
        CGAL_precondition_msg (!m_poly_traits.
                               segment_traits_2()->equal_2_object()(p,q),
                               "Cannot construct a degenerated segment");
        return Curve_2(Segment_2(p,q));
      }

      /* Returns a polyline consists of one given segment. */
      Curve_2 operator()(const Segment_2& seg) const
      {
        return Curve_2(seg);
      }

      /* Construct a well-oriented polyline from a range of either
         `SegmentTraits::Point_2` or `SegmentTraits::Segment_2`. */
      template <typename ForwardIterator>
      Curve_2 operator()(ForwardIterator begin, ForwardIterator end) const
      {
        typedef typename std::iterator_traits<ForwardIterator>::value_type VT;
        typedef typename boost::is_same<VT, Point_2>::type Is_point;
        // Dispatch the range to the appropriate implementation.
        return constructor_impl(begin, end, Is_point());
      }

      /*! Construction of a polyline from a range of points.
       * \pre The range contains at least two points
       * \pre Consecutive points are disjoint.
       * \return Well-oriented polyline connecting the given
       *         points. The order of the vertices is determined by
       *         their order in the range.  Furthermore, the
       *         orientation of the polyline is induced by their
       *         order.
       */
      template <typename ForwardIterator>
      Curve_2 constructor_impl(ForwardIterator begin, ForwardIterator end,
                               boost::true_type) const
      {
        // Container of the segments to be created.
        std::vector<Segment_2> segs;

        // The range must contain at least two points.
        CGAL_precondition_msg(std::distance(begin,end)>1,
                              "Range of points must contain at least 2 points");
        CGAL_precondition_code
          (
           typename Segment_traits_2::Equal_2 equal =
           m_poly_traits.segment_traits_2()->equal_2_object();
           );
        ForwardIterator curr = begin;
        ForwardIterator next = curr;
        ++next;
        while (next != end) {
          CGAL_precondition_msg(!equal(*curr,*next),
                                "Cannot construct a degenerated segment");
          segs.push_back(Segment_2(*curr,*next));
          ++next;
          ++curr;
        }

        return Curve_2(segs.begin(), segs.end());
      }

      /*! Construction implementation from a range of segments.
       *
       *  Note that the segments in the range are NOT necessarily x-monotone,
       *  thus it is impossible to test (even in precondition) whether the input
       *  forms a continuous and well oriented polyline.
       *  \pre Range should contain at least one segment.
       */
      template <typename ForwardIterator>
      Curve_2 constructor_impl(ForwardIterator begin, ForwardIterator end,
                               boost::false_type) const
      {
        // Range has to contain at least one segment
        CGAL_precondition(begin != end);

        return Curve_2(begin, end);
      }
    };

    /*! Get a Construct_curve_2 functor object. */
    Construct_curve_2 construct_curve_2_object() const
    { return Construct_curve_2(*this); }

    class Construct_x_monotone_curve_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>       Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Construct_x_monotone_curve_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits) {}

      /*! Returns an x-monotone polyline connecting the two given endpoints.
       * \param p The first point.
       * \param q The second point.
       * \pre p and q must not be the same.
       * \return A segment connecting p and q.
       */
      X_monotone_curve_2 operator()(const Point_2& p, const Point_2& q) const
      {
        CGAL_precondition_code
          (typename Segment_traits_2::Equal_2 equal =
             m_poly_traits.segment_traits_2()->equal_2_object(););
        CGAL_precondition_msg
          (!equal(p,q),
           "Cannot construct a degenerated segment as a polyline");
        X_monotone_segment_2 seg(p, q);

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
        if (m_poly_traits.segment_traits_2()->compare_xy_2_object()(p,q) ==
          LARGER)
          seg = m_poly_traits.segment_traits_2()->
            construct_opposite_2_object()(seg);
#endif

        return X_monotone_curve_2(seg);
      }

      /*! Returns an x-monotone polyline consists of one given segment.
       * \param seg input segment.
       * \pre seg is not degenerated.
       * \return An x-monotone polyline with one segment, namely seg.
       */
      X_monotone_curve_2 operator()(const X_monotone_segment_2& seg) const
      {
        CGAL_precondition_code
          (
           /*
            * Test that the segment is not degenerated. We do this test
            * independently from the SegmentTraits in use, as we do not allow
            * a polyline with degenerated segments.
            */
           const Segment_traits_2* seg_traits =
           m_poly_traits.segment_traits_2();
           typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
           seg_traits->construct_min_vertex_2_object();
           typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
           seg_traits->construct_max_vertex_2_object();
           typename Segment_traits_2::Equal_2 equal =
           seg_traits->equal_2_object();

           CGAL_precondition_msg(!equal(get_min_v(seg),get_max_v(seg)),
                                 "Cannot construct a degenerated segment");
           );

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
        if (m_poly_traits.segment_traits_2()->
            compare_endpoints_xy_2_object()(seg) == LARGER)
          return X_monotone_segment_2(m_poly_traits.segment_traits_2()->
                                      construct_opposite_2_object()(seg));
#endif

        return X_monotone_curve_2(seg);
      }

      /*!
       * Construct an x-monotone polyline which is well-oriented from a range of
       * elements.
       * \pre Range should from a continuous well-oriented x-monotone polyline.
       */
      template <typename ForwardIterator>
      X_monotone_curve_2 operator()(ForwardIterator begin,
                                    ForwardIterator end) const
      {
        typedef typename std::iterator_traits<ForwardIterator>::value_type VT;
        typedef typename boost::is_same<VT,Point_2>::type Is_point;
        // Dispatch the range to the appropriate implementation.
        return constructor_impl(begin, end, Is_point());
      }

      /*!
       * Construct an x-monotone polyline from a range of points. The
       * polyline may be oriented left-to-right or right-to-left
       * depending on the lexicographical order of the points in the
       * input.
       * \pre Range contains at least two points.
       * \pre No two consecutive points are the same.
       * \pre The points form an continuous well-oriented x-monotone polyline.
       * \post By the construction the returned polyline is well-oriented.
       */
      template <typename ForwardIterator>
      X_monotone_curve_2 constructor_impl(ForwardIterator begin,
                                          ForwardIterator end,
                                          boost::true_type) const
      {
        // Vector of the segments to be constructed from the range of points
        std::vector<X_monotone_segment_2> segs;
        // Make sure the range of points contains at least two points.
        ForwardIterator ps = begin;
        CGAL_precondition(ps != end);
        ForwardIterator pt = ps;
        ++pt;
        CGAL_precondition_msg((pt != end),
                              "Range of points must contain at least 2 points");

        CGAL_precondition_code
          (
           const Segment_traits_2* seg_traits =
             m_poly_traits.segment_traits_2();
           // Initialize two comparison functors
           typename Segment_traits_2::Compare_x_2 compare_x =
             seg_traits->compare_x_2_object();
           typename Segment_traits_2::Compare_xy_2 compare_xy =
             seg_traits->compare_xy_2_object();
           // Make sure there is no changed of directions.
           // Saves the comp_x between the first two points
           const Comparison_result cmp_x_res = compare_x(*ps, *pt);
           // Save the comp_xy between the first two points
           const Comparison_result cmp_xy_res = compare_xy(*ps, *pt);
           );

        // Assure that the first two points are not the same.
        // Note that this also assures that non of the consecutive
        // points are equal in the whole range.
        CGAL_precondition(cmp_xy_res != EQUAL);

        while (pt != end) {
          CGAL_precondition(compare_xy(*ps, *pt) == cmp_xy_res);
          CGAL_precondition(compare_x(*ps, *pt) == cmp_x_res);

          segs.push_back(X_monotone_segment_2(*ps,*pt));
          ++ps; ++pt;
        }

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
        if (m_poly_traits.segment_traits_2()->
          compare_endpoints_xy_2_object()(*segs.begin()) == LARGER)
          {
            X_monotone_curve_2 xcv(segs.begin(), segs.end());
            return m_poly_traits.construct_opposite_2_object()(xcv);
          }
#endif

        return X_monotone_curve_2(segs.begin(), segs.end());
      }

      /*! Returns an x-monotone polyline from a range of segments.
       * \param begin An iterator pointing to the first segment in the range.
       * \param end An iterator pointing to the past-the-end segment
       * in the range.
       * \pre The range contains at least one segment.
       * \pre Segments correspond to a well-oriented polyline. That
       *      is, the target of the i-th segment is an source of the
       *      (i+1)th segment.
       * \pre The sequence of segments in the range forms a weak x-monotone
       *      polyline.
       * \pre The container should support bidirectional iteration.
       * \return A continuous, well-oriented x-monotone polyline which
       *         is directed either left-to-right or right-to-left
       *         depending on the segments in the input.
       */
      template <typename ForwardIterator>
      X_monotone_curve_2 constructor_impl(ForwardIterator begin,
                                          ForwardIterator end,
                                          boost::false_type) const
      {
        CGAL_precondition_msg
          (
           begin != end,
           "Input range of segments has to contain at least"
           "one segment"
           );

        CGAL_precondition_code
          (
           const Segment_traits_2* seg_traits =
             m_poly_traits.segment_traits_2();
           typename Segment_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
             seg_traits->compare_endpoints_xy_2_object();
           typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
             seg_traits->construct_min_vertex_2_object();
           typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
             seg_traits->construct_max_vertex_2_object();
           typename Segment_traits_2::Equal_2 equal =
             seg_traits->equal_2_object();

           ForwardIterator curr = begin;
           ForwardIterator next = begin;
           ++next;
           );


        CGAL_precondition_msg
          (
           (next != end) || !equal(get_max_v(*curr),get_min_v(*curr)),
           "Cannot construct a polyline with degenerated segment"
           );

        CGAL_precondition_code
          (
           // Range contains at least two segments

           Comparison_result init_dir = cmp_seg_endpts(*curr);
           while (next != end){
             CGAL_precondition_msg
               (!equal(get_min_v(*next),get_max_v(*next)),
                "Cannot construct a polyline with degenerated segment"
                );
             CGAL_precondition_msg
               (
                init_dir == cmp_seg_endpts(*next),
                "Segments must form x-monotone polyline"
                );
             if (init_dir == SMALLER){
               CGAL_precondition_msg
                 (
                  equal(get_max_v(*curr),get_min_v(*next)),
                  "Segments should concatenate in source->target manner"
                  );
             }
             else{
               CGAL_precondition_msg
                 (
                  equal(get_min_v(*curr),get_max_v(*next)),
                  "Segments should concatenate in source->target manner"
                  );
             }
             ++curr;
             ++next;
           }
           );

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
        if (m_poly_traits.segment_traits_2()->
          compare_endpoints_xy_2_object()(*begin) == LARGER)
          {
            X_monotone_curve_2 xcv(begin, end);
            return m_poly_traits.construct_opposite_2_object()(xcv);
          }
#endif

        return X_monotone_curve_2(begin, end);
      }
    };

    /*! Get a Construct_x_monotone_curve_2 functor object. */
    Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const
    { return Construct_x_monotone_curve_2(*this); }

    /*! A function object that obtains the parameter space of a geometric
     * entity along the x-axis
     */
    class Parameter_space_in_x_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Parameter_space_in_x_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits)
      {}

      /*! Obtains the parameter space at the end of a curve along the x-axis .
       * Note that if the curve-end coincides with a pole, then unless the curve
       * coincides with the identification curve, the curve-end is considered to
       * be approaching the boundary, but not on the boundary.
       * If the curve coincides with the identification curve, it is assumed to
       * be smaller than any other object.
       * \param xcv the curve
       * \param ce the curve-end indicator:
       *     ARR_MIN_END - the minimal end of xc or
       *     ARR_MAX_END - the maximal end of xc
       * \return the parameter space at the ce end of the curve xcv.
       *   ARR_LEFT_BOUNDARY  - the curve approaches the identification curve
       *                        from the right at the curve left end.
       *   ARR_INTERIOR       - the curve does not approache the identification
       *                        curve.
       *   ARR_RIGHT_BOUNDARY - the curve approaches the identification curve
       *                        from the left at the curve right end.
       * \pre xcv does not coincide with the vertical identification curve.
       */
      Arr_parameter_space operator()(const X_monotone_curve_2& xcv,
                                     Arr_curve_end ce) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        Comparison_result direction =
          seg_traits->compare_endpoints_xy_2_object()(xcv[0]);
        const X_monotone_segment_2& xs =
          (((direction == SMALLER) && (ce == ARR_MAX_END)) ||
           ((direction == LARGER) && (ce == ARR_MIN_END))) ?
          xcv[0] : xcv[xcv.number_of_segments()-1];
        return seg_traits->parameter_space_in_x_2_object()(xs, ce);
      }

      /*! Obtains the parameter space at a point along the x-axis.
       * \param p the point.
       * \return the parameter space at p.
       * \pre p does not lie on the vertical identification curve.
       */
      Arr_parameter_space operator()(const Point_2 p) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        return seg_traits->parameter_space_in_x_2_object()(p);
      }
    };

    /*! Obtain a Parameter_space_in_x_2 function object */
    Parameter_space_in_x_2 parameter_space_in_x_2_object() const
    { return Parameter_space_in_x_2(*this); }

    /*! A function object that obtains the parameter space of a geometric
     * entity along the y-axis
     */
    class Parameter_space_in_y_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Parameter_space_in_y_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits)
      {}

      /*! Obtains the parameter space at the end of an curve along the y-axis .
       * Note that if the curve-end coincides with a pole, then unless the curve
       * coincides with the identification curve, the curve-end is considered to
       * be approaching the boundary, but not on the boundary.
       * If the curve coincides with the identification curve, it is assumed to
       * be smaller than any other object.
       * \param xcv the curve
       * \param ce the curve-end indicator:
       *     ARR_MIN_END - the minimal end of xcv or
       *     ARR_MAX_END - the maximal end of xcv
       * \return the parameter space at the ce end of the curve xcv.
       *   ARR_BOTTOM_BOUNDARY  - the curve approaches the south pole at the
       *                          curve left end.
       *   ARR_INTERIOR         - the curve does not approache a contraction
       *                          point.
       *   ARR_TOP_BOUNDARY     - the curve approaches the north pole at the
       *                          curve right end.
       * There are no horizontal identification curves!
       */
      Arr_parameter_space operator()(const X_monotone_curve_2& xcv,
                                     Arr_curve_end ce) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        Comparison_result direction =
          seg_traits->compare_endpoints_xy_2_object()(xcv[0]);
        const X_monotone_segment_2& xs =
          (((direction == SMALLER) && (ce == ARR_MAX_END)) ||
           ((direction == LARGER) && (ce == ARR_MIN_END))) ?
          xcv[0] : xcv[xcv.number_of_segments()-1];
        return seg_traits->parameter_space_in_y_2_object()(xs, ce);
      }

      /*! Obtains the parameter space at a point along the y-axis.
       * \param p the point.
       * \return the parameter space at p.
       * \pre p does not lie on the horizontal identification curve.
       * There are no horizontal identification curves!
       */
      Arr_parameter_space operator()(const Point_2 p) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        return seg_traits->parameter_space_in_y_2_object()(p);
      }
    };

    /*! Obtain a Parameter_space_in_y_2 function object */
    Parameter_space_in_y_2 parameter_space_in_y_2_object() const
    { return Parameter_space_in_y_2(*this); }

    /*! A functor that compares the x-coordinate of curve-ends and points on the
     * boundary of the parameter space.
     */
    class Compare_x_on_boundary_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Compare_x_on_boundary_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits)
      {}

      /*! Compare the x-limit of a point with the x-coordinate of an
       * x-curve-end on the boundary.
       * \param point the point.
       * \param xcv the x-curve, the endpoint of which is compared.
       * \param ce the x-curve-end indicator:
       *            ARR_MIN_END - the minimal end of xcv or
       *            ARR_MAX_END - the maximal end of xcv.
       * \return the comparison result:
       *         SMALLER - x(p) < x(xcv, ce);
       *         EQUAL   - x(p) = x(xcv, ce);
       *         LARGER  - x(p) > x(xcv, ce).
       * \pre p lies in the interior of the parameter space.
       * \pre the ce end of the x-curve xcv lies on the top boundary.
       * \pre xcv does not coincide with the vertical identification curve.
       */
      Comparison_result operator()(const Point_2& point,
                                   const X_monotone_curve_2& xcv,
                                   Arr_curve_end ce) const
      {
        const Segment_traits_2* seg_traits =
          m_poly_traits.segment_traits_2();
        Comparison_result direction =
          seg_traits->compare_endpoints_xy_2_object()(xcv[0]);
        const X_monotone_segment_2& xs =
          (((direction == SMALLER) && (ce == ARR_MAX_END)) ||
           ((direction == LARGER) && (ce == ARR_MIN_END))) ?
          xcv[0] : xcv[xcv.number_of_segments()-1];
        return seg_traits->compare_x_on_boundary_2_object()(point, xs, ce);
      }

      /*! Compare the x-coordinates of 2 curve-ends near the boundary of the
       * parameter space.
       * \param xcv1 the first curve.
       * \param ce1 the first curve-end indicator:
       *            ARR_MIN_END - the minimal end of xcv1 or
       *            ARR_MAX_END - the maximal end of xcv1.
       * \param xcv2 the second curve.
       * \param ce2 the second  curve-end indicator:
       *            ARR_MIN_END - the minimal end of xcv2 or
       *            ARR_MAX_END - the maximal end of xcv2.
       * \return the second comparison result:
       *         SMALLER - x(xcv1, ce1) < x(xcv2, ce2);
       *         EQUAL   - x(xcv1, ce1) = x(xcv2, ce2);
       *         LARGER  - x(xcv1, ce1) > x(xcv2, ce2).
       * \pre the ce1 end of the curve xcv1 lies on a pole (implying ce1 is
       *      vertical).
       * \pre the ce2 end of the curve xcv2 lies on a pole (implying ce2 is
       *      vertical).
       * \pre xcv1 does not coincide with the vertical identification curve.
       * \pre xcv2 does not coincide with the vertical identification curve.
       */
      Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                   Arr_curve_end ce1,
                                   const X_monotone_curve_2& xcv2,
                                   Arr_curve_end ce2) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        Comparison_result direction1 =
          seg_traits->compare_endpoints_xy_2_object()(xcv1[0]);
        const X_monotone_segment_2& xs1 =
          (((direction1 == SMALLER) && (ce1 == ARR_MAX_END)) ||
           ((direction1 == LARGER) && (ce1 == ARR_MIN_END))) ?
          xcv1[0] : xcv1[xcv1.number_of_segments()-1];
        Comparison_result direction2 =
          seg_traits->compare_endpoints_xy_2_object()(xcv2[0]);
        const X_monotone_segment_2& xs2 =
          (((direction2 == SMALLER) && (ce2 == ARR_MAX_END)) ||
           ((direction2 == LARGER) && (ce2 == ARR_MIN_END))) ?
          xcv2[0] : xcv2[xcv2.number_of_segments()-1];
        return seg_traits->compare_x_on_boundary_2_object()(xs1, ce1, xs2, ce2);
      }
    };

    /*! Obtain a Compare_x_on_boundary_2 function object. */
    Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const
    { return Compare_x_on_boundary_2(*this); }

    /*! A functor that compares the y-coordinate of two given points
     * that lie on the vertical identification curve.
     */
    class Compare_y_on_boundary_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Compare_y_on_boundary_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits)
      {}

      /*! Compare the y-coordinate of two given points that lie on the vertical
       * identification curve.
       * \param p1 the first point.
       * \param p2 the second point.
       * \return SMALLER - p1 is lexicographically smaller than p2;
       *         EQUAL   - p1 and p2 coincides;
       *         LARGER  - p1 is lexicographically larger than p2;
       * \pre p1 lies on the vertical identification curve.
       * \pre p2 lies on the vertical identification curve.
       */
      Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        return seg_traits->compare_y_on_boundary_2_object()(p1, p2);
      }
    };

    /*! Obtain a Compare_y_on_boundary_2 function object */
    Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const
    { return Compare_y_on_boundary_2(*this); }

    /*! A functor that compares the y-coordinates of curve-ends near the
     * boundary of the parameter space.
     */
    class Compare_y_near_boundary_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Compare_y_near_boundary_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits)
      {}

      /*! Compare the y-coordinates of 2 curves at their ends near the boundary
       * of the parameter space.
       * \param xcv1 the first curve.
       * \param xcv2 the second curve.
       * \param ce the curve-end indicator:
       *     ARR_MIN_END - the minimal end or
       *     ARR_MAX_END - the maximal end
       * \return the second comparison result.
       * \pre the ce ends of the curves xcv1 and xcv2 lie either on the left
       *      boundary or on the right boundary of the parameter space (implying
       *      that they cannot be vertical).
       * There is no horizontal identification curve!
       */
      Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                   const X_monotone_curve_2& xcv2,
                                   Arr_curve_end ce) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        Comparison_result direction1 =
          seg_traits->compare_endpoints_xy_2_object()(xcv1[0]);
        const X_monotone_segment_2& xs1 =
          (((direction1 == SMALLER) && (ce == ARR_MAX_END)) ||
           ((direction1 == LARGER) && (ce == ARR_MIN_END))) ?
          xcv1[0] : xcv1[xcv1.number_of_segments()-1];
        Comparison_result direction2 =
          seg_traits->compare_endpoints_xy_2_object()(xcv2[0]);
        const X_monotone_segment_2& xs2 =
          (((direction2 == SMALLER) && (ce == ARR_MAX_END)) ||
           ((direction2 == LARGER) && (ce == ARR_MIN_END))) ?
          xcv2[0] : xcv2[xcv2.number_of_segments()-1];
        return seg_traits->compare_y_near_boundary_2_object()(xs1, xs2, ce);
      }
    };

    /*! Obtain a Compare_y_near_boundary_2 function object */
    Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
    { return Compare_y_near_boundary_2(*this); }

    /*! A functor that indicates whether a geometric object lies on the
     * vertical identification arc.
     */
    class Is_on_y_identification_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Is_on_y_identification_2(const Polyline_traits_2& traits) :
        m_poly_traits(traits)
      {}

      /*! Determine whether a point lies in the vertical boundary.
       * \param p the point.
       * \return a Boolean indicating whether p lies in the vertical boundary.
       */
      bool operator()(const Point_2& p) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        return seg_traits->is_on_y_identification_2_object()(p);
      }

      /*! Determine whether an x-monotone curve lies in the vertical boundary.
       * \param xcv the x-monotone curve.
       * \return a Boolean indicating whether xcv lies in the vertical boundary.
       */
      bool operator()(const X_monotone_curve_2& xcv) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        typename X_monotone_curve_2::Segment_const_iterator it;
        for (it = xcv.begin_segments(); it != xcv.end_segments(); ++it)
          if (! seg_traits->is_on_y_identification_2_object()(*it))
            return false;
        return true;
      }
    };

    /*! Obtain a Is_on_y_identification_2 function object */
    Is_on_y_identification_2 is_on_y_identification_2_object() const
    { return Is_on_y_identification_2(*this); }

  private:
    /*
     * Roadmap: locate() should return an iterator to the located segment
     */
    /*!
     * Return the index of the segment in the polyline that contains the
     * point q in its x-range. The function performs a binary search, so if the
     * point q is in the x-range of the polyline with n segments, the segment
     * containing it can be located in O(log n) operations.
     * \param cv The polyline curve.
     * \param q The point.
     * \return An index i such that q is in the x-range of cv[i].
     *         If q is not in the x-range of cv, returns INVALID_INDEX.
     */
    template <typename Compare>
    std::size_t locate_gen(const X_monotone_curve_2& cv, Compare compare) const
    {
      // The direction of cv. SMALLER means left-to-right and
      // otherwise right-to-left
      Comparison_result direction = segment_traits_2()->
        compare_endpoints_xy_2_object()(cv[0]);
      std::size_t from, to;
      if (direction == SMALLER) {
        from = 0;
        to = cv.number_of_segments() - 1;
      }
      else {
        from = cv.number_of_segments() - 1;
        to = 0;
      }

      // Test if q is one of cv's end points
      Comparison_result res_from = compare(cv[from], ARR_MIN_END);
      if (res_from == EQUAL) return from;

      Comparison_result res_to = compare(cv[to], ARR_MAX_END);
      if (res_to == EQUAL) return to;

      // Check whether the point is either lexicographically to the left of
      // the curve or lexicographically to the right of the curve.
      if (res_to == res_from)
        // If the x-monotone polyline is vertical, return the index of the
        // segment that is closest to the point. Otherwise, the point is not
        // in the x-range of the polyline.
        return (is_vertical_2_object()(cv)) ?
          ((res_to == SMALLER) ? from : to) : std::size_t(INVALID_INDEX);

      // Perform a binary search to locate the segment that contains q in its
      // range:
      while (((direction == SMALLER) && (to > from)) ||
             ((direction == LARGER) && (to < from)))
      {
        std::size_t mid = (from + to) / 2;
        if (((direction == SMALLER) && (mid > from)) ||
            ((direction == LARGER) && (mid < from)))
        {
          Comparison_result res_mid = compare(cv[mid], ARR_MIN_END);
          if (res_mid == EQUAL) {
            // Ensure that the returned segment contains the query point
            // on its right end (if possible)
            if ((direction == SMALLER) && (mid > 0)) --mid;
             else if ((direction == LARGER) &&
                      ((mid + 1) < cv.number_of_segments()))
               ++mid;
            return mid;
          }
          if (res_mid == res_from) from = mid;
          else to = (direction == SMALLER) ? mid - 1 : mid + 1;
        }
        else {
          CGAL_assertion(((direction == SMALLER) && (mid < to)) ||
                         ((direction == LARGER) && (mid > to)));
          Comparison_result res_mid = compare(cv[mid], ARR_MAX_END);
          if (res_mid == EQUAL) return mid;
          if (res_mid == res_to) to = mid;
          else from = (direction == SMALLER) ? mid + 1 : mid - 1;
        }
      }
      // In case (from == to), and we know that the polyline contains the q:
      CGAL_assertion(from == to);
      return from;
    }

    // A utility class that compare a curve-end with a point.
    template <typename Comparer>
    class Compare_points {
    private:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;
      /*! The polyline traits (in case it has state) */
      const Polyline_traits_2& m_poly_traits;

      const Point_2& m_point;

      Comparer m_compare;

    public:
      // Constructor
      Compare_points(const Polyline_traits_2& traits, Comparer compare,
                     const Point_2& p) :
        m_poly_traits(traits),
        m_point(p),
        m_compare(compare)
      {}

      // Compare the given curve-end with the stored point.
      Comparison_result operator()(const X_monotone_segment_2& xs,
                                   Arr_curve_end ce)
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        const Point_2& p = (ce == ARR_MAX_END) ?
          seg_traits->construct_max_vertex_2_object()(xs) :
          seg_traits->construct_min_vertex_2_object()(xs);
        return m_compare(p, m_point);
      }
    };

    // A utility class that compare two curve-ends.
    template <typename Comparer>
    class Compare_curve_ends {
    private:
      const X_monotone_segment_2& m_x_monotone_segment;

      Arr_curve_end m_curve_end;

      Comparer m_compare;

    public:
      // Constructor
      Compare_curve_ends(Comparer compare,
                         const X_monotone_segment_2& xs, Arr_curve_end ce) :
        m_x_monotone_segment(xs),
        m_curve_end(ce),
        m_compare(compare)
      {}

      // Compare the given curve-end with the stored point.
      Comparison_result operator()(const X_monotone_segment_2& xs,
                                   Arr_curve_end ce)
      { return m_compare(xs, ce, m_x_monotone_segment, m_curve_end); }
    };

    //
    std::size_t locate_impl(const X_monotone_curve_2& xcv,
                             const X_monotone_segment_2& xs,
                             Arr_curve_end ce,
                             Arr_not_all_sides_oblivious_tag) const
    {
      const Segment_traits_2* seg_traits = segment_traits_2();
      if (seg_traits->is_vertical_2_object()(xcv[0])) {
        // Verify that q has the same x-coord as xcv (which is vertical)
        Compare_x_2 compare_x = compare_x_2_object();
        Comparison_result res = compare_x(xcv[0], ARR_MIN_END, xs, ce);
        if (res != EQUAL) return INVALID_INDEX;

        Compare_curve_ends<Compare_xy_2> compare(compare_xy_2_object(), xs, ce);
        return locate_gen(xcv, compare);
      }

      Compare_curve_ends<Compare_x_2> compare(compare_x_2_object(), xs, ce);
      return locate_gen(xcv, compare);
    }

    //
    std::size_t locate_impl(const X_monotone_curve_2& xcv,
                             const X_monotone_segment_2& xs,
                             Arr_curve_end ce,
                             Arr_all_sides_oblivious_tag) const
    {
      const Segment_traits_2* seg_traits = segment_traits_2();
      const Point_2& p = (ce == ARR_MAX_END) ?
        seg_traits->construct_max_vertex_2_object()(xs) :
        seg_traits->construct_min_vertex_2_object()(xs);
      return locate(xcv, p);
    }

    //
    std::size_t locate(const X_monotone_curve_2& xcv, const Point_2& q) const
    {
      const Segment_traits_2* seg_traits = segment_traits_2();
      if (seg_traits->is_vertical_2_object()(xcv[0])) {
        // Verify that q has the same x-coord as cv (which is vertical)
        typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
          seg_traits->construct_min_vertex_2_object();
        typename Segment_traits_2::Compare_x_2 compare_x =
          seg_traits->compare_x_2_object();
        Comparison_result res = compare_x(min_vertex(xcv[0]), q);
        if (res != EQUAL) return INVALID_INDEX;

        Compare_points<Compare_xy_2> compare(seg_traits, compare_xy_2_object(), q);
        return locate_gen(xcv, compare);
      }

      Compare_points<Compare_x_2> compare(seg_traits, compare_x_2_object(), q);
      return locate_gen(xcv, compare);
    }

    /*!
     * Find the index of the segment in the polyline that is defined to the
     * left(or to the right) of the point q.
     * \param cv The polyline curve.
     * \param q The point.
     * \param to_right(true) if we wish to locate a segment to the right of q,
     *               (false) if we wish to locate a segment to its right.
     * \return An index i such that segments[i] is defined to the left(or to the
     *         right) of q, or INVALID_INDEX if no such segment exists.
     */
    std::size_t locate_side(const X_monotone_curve_2& cv,
                             const Point_2& q, const bool& to_right) const
    {
      // First locate a segment segments[i] that contains q in its x-range.
      std::size_t i = locate(cv, q);
      if (i == INVALID_INDEX) return INVALID_INDEX;

      typename Segment_traits_2::Equal_2 equal =
        segment_traits_2()->equal_2_object();
      typename Segment_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
        segment_traits_2()->compare_endpoints_xy_2_object();
      typename Segment_traits_2::Compare_x_2 comp_x =
        segment_traits_2()->compare_x_2_object();
      typename Segment_traits_2::Is_vertical_2 is_vert =
        segment_traits_2()->is_vertical_2_object();
      typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
        segment_traits_2()->construct_max_vertex_2_object();
      typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
        segment_traits_2()->construct_min_vertex_2_object();

      Comparison_result direction = cmp_seg_endpts(cv[i]);

      if ((!is_vert(cv[0]) && (comp_x(get_min_v(cv[i]), q) == EQUAL)) ||
          (is_vert(cv[0]) && equal(get_min_v(cv[i]), q))){
        // q is the left endpoint of the i'th segment:
        if (to_right)
          return i;
        else {
          // to_left
          if (direction == SMALLER)
            if (i == 0)
              return INVALID_INDEX;
            else
              return i - 1;
          else {
            if (i == cv.number_of_segments()-1)
              return INVALID_INDEX;
            else
              return i+1;
          }
        }
      }

      if ((!is_vert(cv[0]) && (comp_x(get_max_v(cv[i]), q) == EQUAL)) ||
          (is_vert(cv[0]) && equal(get_max_v(cv[i]), q)))
      {
        // q is the right endpoint of the i'th segment:
        if (!to_right)
          return i;
        else {
          if (direction == SMALLER){
            if (i == (cv.number_of_segments() - 1))
              return INVALID_INDEX;
            else
              return i + 1;
          }
          else {
            if (i == 0)
              return INVALID_INDEX;
            else
              return i-1;
          }
        }
      }

      // In case q is in cv[i]'s interior:
      return i;
    }
  };
} //namespace CGAL

#endif
