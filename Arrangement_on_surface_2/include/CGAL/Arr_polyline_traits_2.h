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

/*
 * TODO: Complete the documentation of the changes derived from the cleaning
 *       In particular, doxygen only the things that have to be exposed
 *       to the user.
 */

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
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_geometry_traits/Polyline_2.h>

namespace CGAL {

  template <typename SegmentTraits_2 = Arr_segment_traits_2<> >
  class Arr_polyline_traits_2 {
  public:
    typedef SegmentTraits_2                          Segment_traits_2;

    // Tag definitions:
    typedef Tag_true                                 Has_left_category;
    typedef Tag_true                                 Has_merge_category;
    typedef Tag_false                                Has_do_intersect_category;

    typedef Arr_oblivious_side_tag                   Left_side_category;
    typedef Arr_oblivious_side_tag                   Bottom_side_category;
    typedef Arr_oblivious_side_tag                   Top_side_category;
    typedef Arr_oblivious_side_tag                   Right_side_category;

  private:
    typedef Arr_polyline_traits_2<Segment_traits_2>  Self;

    // Data members:
    const Segment_traits_2*   m_seg_traits;    // The base segment-traits class.
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
    typedef typename Segment_traits_2::Compare_x_2        Compare_x_2;

    /*! Get a Compare_x_2 functor object. */
    Compare_x_2 compare_x_2_object() const
    { return segment_traits_2()->compare_x_2_object(); }

    /*! Compare two points lexigoraphically: by x, then by y. */
    typedef typename Segment_traits_2::Compare_xy_2       Compare_xy_2;

    /*! Get a Compare_xy_2 functor object. */
    Compare_xy_2 compare_xy_2_object() const
    { return segment_traits_2()->compare_xy_2_object(); }

    ///@}

    /// \name Basic predicate functors(based on the segment traits).
    //@{

    class Number_of_points_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Geometry_traits_2;
      /* The polyline traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:
      /* Constructor. */
      Number_of_points_2(const Geometry_traits_2& traits) :
        m_poly_traits(traits)
      {}

      const int operator()(const Curve_2& cv) const
      {
        int num_seg = cv.number_of_segments();
        return (num_seg == 0) ? 0 : num_seg + 1;
      }
    };

    Number_of_points_2 number_of_points_2_object() const
    { return Number_of_points_2(this); }


    class Construct_min_vertex_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Geometry_traits_2;
      /*! The polyline traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:
      /* Constructor. */
      Construct_min_vertex_2(const Geometry_traits_2& traits) :
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
      typedef Arr_polyline_traits_2<Segment_traits_2>     Geometry_traits_2;
      /*! The polyline traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Construct_max_vertex_2(const Geometry_traits_2& traits) :
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
      typedef Arr_polyline_traits_2<Segment_traits_2>     Geometry_traits_2;
      /*! The polyline traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Is_vertical_2(const Geometry_traits_2& traits) : m_poly_traits(traits) {}

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
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>       Geometry_traits_2;
      /*! The polyline traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Compare_y_at_x_2(const Geometry_traits_2& traits) :
        m_poly_traits(traits) {}

      /*!
       * Return the location of the given point with respect to the input curve.
       * \param cv The polyline curve.
       * \param p The point.
       * \pre p is in the x-range of cv.
       * \return SMALLER if y(p) < cv(x(p)), i.e. the point is below the curve;
       *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
       *         EQUAL if p lies on the curve.
       */
      Comparison_result operator()(const Point_2& p,
                                   const X_monotone_curve_2& cv) const
      {
        // Get the index of the segment in cv containing p.
        unsigned int i = // What is a smart way to call _locate()?
          m_poly_traits.locate(cv, p);
        CGAL_precondition(i != INVALID_INDEX);

        // Compare the segment cv[i] and p.
        return m_poly_traits.segment_traits_2()->
          compare_y_at_x_2_object()(p, cv[i]);
      }
    };

    /*! Get a Compare_y_at_x_2 functor object. */
    Compare_y_at_x_2 compare_y_at_x_2_object() const
    { return Compare_y_at_x_2(*this); }

    class Compare_y_at_x_left_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>       Geometry_traits_2;
      /*! The polyline traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Compare_y_at_x_left_2(const Geometry_traits_2& traits) :
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
        unsigned int i1=m_poly_traits.locate_side(cv1, p, false);
        unsigned int i2=m_poly_traits.locate_side(cv2, p, false);

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
      typedef Arr_polyline_traits_2<Segment_traits_2>       Geometry_traits_2;
      /*! The polyline traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Compare_y_at_x_right_2(const Geometry_traits_2& traits) :
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
        unsigned int i1=m_poly_traits.locate_side(cv1, p, true);
        unsigned int i2=m_poly_traits.locate_side(cv2, p, true);

        CGAL_precondition(i1 != INVALID_INDEX);
        CGAL_precondition(i2 != INVALID_INDEX);

        // Compare cv1[i1] and cv2[i2] at p's right.
        return m_poly_traits.segment_traits_2()->
          compare_y_at_x_right_2_object()(cv1[i1],cv2[i2], p);
      }
    };

    /*! Get a Compare_y_at_x_right_2 functor object. */
    Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
    { return Compare_y_at_x_right_2(*this); }

    class Equal_2 {
    protected:

      typedef Arr_polyline_traits_2<Segment_traits_2>     Geometry_traits_2;
      /*! The polyline traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:

      /*! Constructor. */
      Equal_2(const Geometry_traits_2& poly_tr) : m_poly_traits(poly_tr) {}

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
        typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
          seg_traits->construct_min_vertex_2_object();
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
        bool res = equal(xpoly_min_v(cv1),xpoly_min_v(cv2));
        if (!res) return false;
        res = equal(xpoly_max_v(cv1),xpoly_max_v(cv2));
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
        Point_2 point1, point2;
        Comparison_result res_x;
        Comparison_result res_y_at_x;
        unsigned int i = 0, j = 0;
        unsigned int n1 = cv1.number_of_segments();
        unsigned int n2 = cv2.number_of_segments();
        Comparison_result is_cv1_left_to_right = comp_endpt(cv1[0]);
        Comparison_result is_cv2_left_to_right = comp_endpt(cv2[0]);

        while ((i < n1-1) || (j < n2-1)) {
          int cv1_seg_ind,cv2_seg_ind;
          if (is_cv1_left_to_right == SMALLER){
            cv1_seg_ind = i;
            point1 = max_vertex(cv1[cv1_seg_ind]);
          }
          else{
            cv1_seg_ind=n1-1-i;
            point1 = max_vertex(cv1[cv1_seg_ind]);
          }
          if (is_cv2_left_to_right == SMALLER){
            cv2_seg_ind=j;
            point2 = max_vertex(cv2[cv2_seg_ind]);
          }
          else{
            cv2_seg_ind=n2-1-j;
            point2 = max_vertex(cv2[cv2_seg_ind]);
          }

          res = equal(point1, point2);
          // Easy case - the two points are equal
          if (res) {
            ++i;
            ++j;
          }
          else {
            res_x = compare_x(point1,point2);
            // Check if the different point is a collinear point situated on
            // the line between its two neighbors.
            if (res_x == SMALLER) {
              res_y_at_x = compare_y_at_x(point1,cv2[cv2_seg_ind]);
              if (res_y_at_x == EQUAL)
                ++i;
              else
                return false;
            }
            else if(res_x == LARGER) {
              res_y_at_x = compare_y_at_x(point2,cv1[cv1_seg_ind]);
              if (res_y_at_x == EQUAL)
                ++j;
              else
                return false;
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
      typedef Arr_polyline_traits_2<Segment_traits_2> Geometry_traits_2;
      /*! The traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Compare_endpoints_xy_2(const Geometry_traits_2& traits) :
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
      typedef Arr_polyline_traits_2<Segment_traits_2> Geometry_traits_2;
      /*! The traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:
      /*! Constructor */
      Construct_opposite_2(const Geometry_traits_2& traits) :
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
      typedef Arr_polyline_traits_2<Segment_traits_2>     Geometry_traits_2;
      /*! The traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Make_x_monotone_2(const Geometry_traits_2& traits) :
        m_poly_traits(traits) {}

      /*!
       * Cut the given curve into x-monotone subcurves and insert them into the
       * given output iterator. As segments are always x_monotone, only one
       * object will be contained in the iterator.
       * \param cv The curve.
       * \param oi The output iterator, whose value-type is Object. The output
       *           object is a wrapper of either an X_monotone_curve_2, or - in
       *           case the input segment is degenerate - a Point_2 object.
       * \return The past-the-end iterator.
       */
      template<class OutputIterator>
      OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const
      {
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
        bool check = CGAL::assign(x_seg, *it);
        CGAL_assertion(check);

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

        // The polyline consists of at least 2 x-monotone segments:
        Push_back_2 push_back = m_poly_traits.push_back_2_object();
        
        typename Segment_traits_2::Is_vertical_2 is_seg_vertical =
          m_poly_traits.segment_traits_2()->is_vertical_2_object();

        bool is_start_vertical = is_seg_vertical(x_seg);
        Comparison_result start_dir = cmp_seg_endpts(x_seg);

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
        if (cmp_seg_endpts(x_seg) == LARGER)
          x_seg = ctr_seg_opposite(x_seg);
#endif
        X_monotone_curve_2 x_polyline = ctr_x_curve(x_seg);

        for (++it; it != x_seg_objects.end(); ++it) {
          X_monotone_segment_2 x_seg;
          bool check = CGAL::assign(x_seg, *it);
          CGAL_assertion(check);
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
            if (cmp_seg_endpts(x_seg) == LARGER)
              x_seg = ctr_seg_opposite(x_seg);
#endif
            push_back(x_polyline, x_seg);
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

    /* see documentation in
     * ../../doc/Arrangement_on_surface_2/CGAL/Arr_polyline_traits_2.h
     * TODO: Test all the operator()'s. (Don't forget vertical cases!)
     *       Should we add tests of this functor to the test suite?
     */
    class Push_back_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>     Geometry_traits_2;
      /*! The traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Push_back_2(const Geometry_traits_2& traits) : m_poly_traits(traits) {}

      /* see documentation in
         ../../doc/Arrangement_on_surface_2/CGAL/Arr_polyline_traits_2.h */
      void operator()(Curve_2& cv, const Point_2& p) const
      {
        int num_seg = cv.number_of_segments();
        CGAL_precondition(num_seg >0);
        int last_seg = num_seg-1;

        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        typename Segment_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
          seg_traits->compare_endpoints_xy_2_object();

        /*
         * Since we assume that the segments of cv are well oriented,
         * pushing a single point to an existing polyline means that we
         * have to append the segment [cv[last_seg].target(),p]. The
         * following test determines which end of the last segment is
         * the target.
         */
        if (cmp_seg_endpts(cv[last_seg]) == SMALLER) {
          typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
            seg_traits->construct_max_vertex_2_object();
          cv.push_back(X_monotone_segment_2(get_max_v(cv[last_seg]),p));
        }
        else {
          typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
            seg_traits->construct_min_vertex_2_object();
          cv.push_back(X_monotone_segment_2(get_min_v(cv[last_seg]),p));
        }
      }


      /* see documentation in
         ../../doc/Arrangement_on_surface_2/CGAL/Arr_polyline_traits_2.h */
      void operator()(Curve_2& cv, const Segment_2& seg) const
      {
        int num_seg = cv.number_of_segments();

        // cv is empty
        if (num_seg == 0) {
          cv.push_back(seg);
          return;
        }

        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
          seg_traits->construct_min_vertex_2_object();
        typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
          seg_traits->construct_max_vertex_2_object();
        typename Segment_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
          seg_traits->compare_endpoints_xy_2_object();
        typename Segment_traits_2::Equal_2 equal =
          seg_traits->equal_2_object();

        Point_2 last_v;

        Comparison_result xcv_dir = cmp_seg_endpts(cv[0]);
        Comparison_result seg_dir = cmp_seg_endpts(seg);

        CGAL_precondition_msg(xcv_dir == seg_dir,
                              "Appended segment must have the same "
                              "orientation as the x-monotone polyline.");

        if (xcv_dir == SMALLER) {
          last_v = get_max_v(cv[num_seg-1]);
          CGAL_precondition(equal(last_v,get_min_v(seg)));
          cv.push_back(seg);
        }
        else {
          last_v = get_min_v(cv[num_seg-1]);
          CGAL_precondition(equal(last_v,get_max_v(seg)));
          cv.push_back(seg);
        }
      }

      /* see documentation in
         ../../doc/Arrangement_on_surface_2/CGAL/Arr_polyline_traits_2.h */
      void operator()(const X_monotone_curve_2& xcv, Point_2& p) const
      {
        int num_seg = xcv.number_of_segments();

        CGAL_precondition(num_seg > 0);

        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
          seg_traits->construct_max_vertex_2_object();

        if (seg_traits->compare_endpoints_xy_2_object()(xcv[0]) == SMALLER) {
          // xcv is oriented left-to-right
          CGAL_precondition_code(
            typename Segment_traits_2::Compare_x_2 comp_x =
            seg_traits->compare_x_2_object();
            CGAL_precondition(comp_x(get_max_v(xcv[num_seg-1]), p) == LARGER);
                                 );
          xcv.push_back(X_monotone_segment_2(get_max_v(xcv[num_seg-1]),p));

        }
        else {
          // xcv is oriented right-to-left
          CGAL_precondition_code(
            typename Segment_traits_2::Compare_x_2 comp_x =
            seg_traits->compare_x_2_object();
            CGAL_precondition(comp_x(get_min_v(xcv[num_seg-1]), p) == SMALLER);
                                 );
          xcv.push_back(X_monotone_segment_2(get_min_v(xcv[num_seg-1]), p));
        }
      }

      /* see documentation in
         ../../doc/Arrangement_on_surface_2/CGAL/Arr_polyline_traits_2.h */
      void operator()(X_monotone_curve_2& xcv,
                      const X_monotone_segment_2& seg) const
      {
        int num_seg = xcv.number_of_segments();

        if (num_seg == 0) {
          xcv.push_back(seg);
          return;
        }

        CGAL_precondition_code(
          const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
          typename Segment_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
          seg_traits->compare_endpoints_xy_2_object();
                               );

        CGAL_precondition_msg(cmp_seg_endpts(xcv[0]) == cmp_seg_endpts(seg),
                              "xcv and seg should have the same orientation");

        CGAL_precondition_code(
          typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
            seg_traits->construct_max_vertex_2_object();
          typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
            seg_traits->construct_min_vertex_2_object();
          typename Segment_traits_2::Compare_xy_2 comp_xy =
            seg_traits->compare_xy_2_object();
          typename Segment_traits_2::Equal_2 equal =
            seg_traits->equal_2_object();
                               );
          CGAL_precondition_msg
            (((cmp_seg_endpts(seg) == SMALLER) &&
              equal(get_max_v(xcv[num_seg-1]), get_min_v(seg)) &&
              (comp_xy(get_min_v(seg), get_max_v(seg)) != EQUAL)) ||
             ((cmp_seg_endpts(seg) == LARGER) &&
              equal(get_min_v(xcv[num_seg-1]), get_max_v(seg)) &&
              (comp_xy(get_min_v(seg), get_max_v(seg)) != EQUAL)),
             "Seg must extend either to the left or to the "
             "right of xcv depending on the orientation");
          
        xcv.push_back(seg);
      }
    };

    /*! Get a Push_Back_2 functor object. */
    Push_back_2 push_back_2_object() const
    { return Push_back_2(*this); }

    class Split_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>       Geometry_traits_2;
      /*! The polyline traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Split_2(const Geometry_traits_2& traits) : m_poly_traits(traits) {}

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
        unsigned int i = m_poly_traits.locate(xcv, p);

        CGAL_precondition(i != INVALID_INDEX);

        // Clear the output curves.
        xcv1.clear();
        xcv2.clear();

        // Push all segments labeled(0, 1, ... , i-1) into xcv1.
        for (int j = 0; j < i; ++j)
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
        unsigned int n = xcv.number_of_segments();

        for (int j = i+1; j < n; ++j)
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
      typedef Arr_polyline_traits_2<Segment_traits_2>       Geometry_traits_2;
      /*! The polyline traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Intersect_2(const Geometry_traits_2& traits) : m_poly_traits(traits) {}

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
      template<class OutputIterator>
      OutputIterator operator()(const X_monotone_curve_2& cv1,
                                const X_monotone_curve_2& cv2,
                                OutputIterator oi) const
      {
        const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
        typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
          seg_traits->construct_min_vertex_2_object();
        typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
          seg_traits->construct_max_vertex_2_object();
        typename Segment_traits_2::Equal_2                equal =
          seg_traits->equal_2_object();
        typename Segment_traits_2::Compare_xy_2           compare_xy =
          seg_traits->compare_xy_2_object();
        typename Segment_traits_2::Intersect_2            intersect =
          seg_traits->intersect_2_object();
        typename Segment_traits_2::Compare_y_at_x_2       compare_y_at_x =
          seg_traits->compare_y_at_x_2_object();
        typename Segment_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
          seg_traits->compare_endpoints_xy_2_object();

        Comparison_result dir1 = cmp_seg_endpts(cv1[0]);
        Comparison_result dir2 = cmp_seg_endpts(cv2[0]);

        const unsigned int n1 = cv1.number_of_segments();
        const unsigned int n2 = cv2.number_of_segments();

        unsigned int i1 = (dir1 == SMALLER) ? 0 : n1-1;
        unsigned int i2 = (dir2 == SMALLER) ? 0 : n2-1;

        X_monotone_curve_2 ocv;           // Used to represent overlaps.

        Comparison_result left_res =
          compare_xy(min_vertex(cv1[i1]), min_vertex(cv2[i2]));

        if (left_res == SMALLER) {
          // cv1's left endpoint is to the left of cv2's left endpoint:
          // Locate the index i1 of the segment in cv1 which contains cv2's
          // left endpoint.
          i1 = m_poly_traits.locate(cv1, min_vertex(cv2[i2]));
          if (i1 == INVALID_INDEX) return oi;

          if (equal(max_vertex(cv1[i1]), min_vertex(cv2[i2]))) {
            if (((dir1 == SMALLER) && (i1 == n1-1)) ||
                ((dir1 == LARGER) && (i1 == 0))){
              // cv1's right endpoint equals cv2's left endpoint
              // Thus we can return this single(!) intersection point
              std::pair<Point_2,Multiplicity>  p(max_vertex(cv1[i1]), 0);
              *oi++ = make_object(p);
              return oi;
            }
            dir1 == SMALLER ? ++i1 : (i1 != 0) ? --i1 : INVALID_INDEX;
            left_res = EQUAL;
          }
        }
        else if (left_res == LARGER) {
          // cv1's left endpoint is to the right of cv2's left endpoint:
          // Locate the index i2 of the segment in cv2 which contains cv1's
          // left endpoint.
          i2 = m_poly_traits.locate(cv2, min_vertex(cv1[i1]));
          if (i2 == INVALID_INDEX) return oi;

          if (equal(max_vertex(cv2[i2]), min_vertex(cv1[i1]))) {
            if (((dir2 == SMALLER) && (i2 == n2-1)) ||
                ((dir2 == LARGER) && (i2 == 0))){
              // cv2's right endpoint equals cv1's left endpoint
              // Thus we can return this single(!) intersection point
              std::pair<Point_2,Multiplicity>  p(max_vertex(cv2[i2]), 0);
              *oi++ = make_object(p);
              return oi;
            }

            dir2 == SMALLER ? ++i2 : (i2 != 0) ? --i2 : INVALID_INDEX;
            left_res = EQUAL;
          }
        }

        // Check if the the left endpoint lies on the other polyline.
        bool left_coincides = (left_res == EQUAL);
        bool left_overlap = false;

        if (left_res == SMALLER)
          left_coincides =
            (compare_y_at_x(min_vertex(cv2[i2]), cv1[i1]) == EQUAL);
        else if (left_res == LARGER)
          left_coincides =
            (compare_y_at_x(min_vertex(cv1[i1]), cv2[i2]) == EQUAL);

        // The main loop: Go simultaneously over both polylines.
        Comparison_result right_res = left_res;
        bool right_coincides = left_coincides;
        bool right_overlap = false;

        while (((dir1==SMALLER) && (dir2 == SMALLER) && (i1 < n1) &&(i2 < n2))||
               ((dir1!=SMALLER) && (dir2 == SMALLER) && (i1 >= 0) &&
                (i1 != INVALID_INDEX) && (i2 < n2)) ||
               ((dir1==SMALLER) && (dir2 != SMALLER) && (i1 < n1) &&
                (i2 >= 0) && (i2 != INVALID_INDEX)) ||
               ((dir1!=SMALLER) && (dir2 != SMALLER) && (i1 >= 0) &&(i2 >= 0)&&
                (i1 != INVALID_INDEX) && (i2 != INVALID_INDEX)))
          {
            right_res = compare_xy(max_vertex(cv1[i1]), max_vertex(cv2[i2]));

            right_coincides = (right_res == EQUAL);
            if (right_res == SMALLER)
              right_coincides = (compare_y_at_x(max_vertex(cv1[i1]),
                                                cv2[i2]) == EQUAL);
            else if (right_res == LARGER)
              right_coincides = (compare_y_at_x(max_vertex(cv2[i2]),
                                                cv1[i1]) == EQUAL);

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
                  std::pair<Point_2, Multiplicity>  p(min_vertex(cv2[i2]), 0);
                  *oi++ = make_object(p);
                }
                else {
                  std::pair<Point_2, Multiplicity>  p(min_vertex(cv1[i1]), 0);
                  *oi++ = make_object(p);
                }
              }
            }

            // Proceed forward.
            if (right_res != SMALLER)
              if (dir2 == SMALLER)
                ++i2;
              else if (i2 == 0)
                i2 = INVALID_INDEX;
              else
                --i2;
            if (right_res != LARGER)
              if (dir1 == SMALLER)
                ++i1;
              else if (i1 == 0)
                i1 = INVALID_INDEX;
              else --i1;
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
                   (i1 == INVALID_INDEX) && (dir1 != SMALLER))
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
      typedef Arr_polyline_traits_2<Segment_traits_2>       Geometry_traits_2;
      /*! The polyline traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Are_mergeable_2(const Geometry_traits_2& traits) :
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

        const unsigned int n1 = cv1.number_of_segments();
        const unsigned int n2 = cv2.number_of_segments();

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
     * A functor that merges two x-monotone arcs into one.
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
          const unsigned int n1 = cv1.number_of_segments();
          const unsigned int n2 = cv2.number_of_segments();
          unsigned int       i;

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
    typedef typename Segment_traits_2::Approximate_number_type
    Approximate_number_type;
    typedef typename Segment_traits_2::Approximate_2    Approximate_2;


    /*! Get an Approximate_2 functor object. */
    Approximate_2 approximate_2_object() const
    { return segment_traits_2()->approximate_2_object(); }

    class Construct_curve_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>       Geometry_traits_2;
      /*! The polyline traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Construct_curve_2(const Geometry_traits_2& traits) :
        m_poly_traits(traits) {}

      /* see documentation in
         ../../doc/Arrangement_on_surface_2/CGAL/Arr_polyline_traits_2.h */
      Curve_2 operator()(const Point_2& p, const Point_2& q) const
      {
        CGAL_precondition_msg (!m_poly_traits.
                               segment_traits_2()->equal_2_object()(p,q),
                               "Cannot construct a degenerated segment");
        return Curve_2(X_monotone_segment_2(p,q));
      }

      /* see documentation in
         ../../doc/Arrangement_on_surface_2/CGAL/Arr_polyline_traits_2.h */
      Curve_2 operator()(const Segment_2& seg) const
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
        return Curve_2(seg);
      }

      /* see documentation in
         ../../doc/Arrangement_on_surface_2/CGAL/Arr_polyline_traits_2.h */
      template <typename ForwardIterator>
      Curve_2 operator()(ForwardIterator begin, ForwardIterator end) const
      {
        typedef typename std::iterator_traits<ForwardIterator>::value_type VT;
        typedef typename boost::is_same<VT,Point_2>::type Is_point;
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
       *  \pre Range should contain at least one segment.
       *  \pre The segments form a continuous polyline.
       *  \pre The polyline is well oriented, that is the target of the i-th
       *       segment should be the source of the (i+1)-th segment.
       *  \return Well-oriented polyline.
       */
      template <typename ForwardIterator>
      Curve_2 constructor_impl(ForwardIterator begin, ForwardIterator end,
                               boost::false_type) const
      {
        // Range has to contain at least one segment
        CGAL_precondition(begin != end);

        ForwardIterator curr = begin;
        ForwardIterator next = curr;

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
           );

        if (++next == end) {
          CGAL_precondition_msg(!equal(get_min_v(*curr),get_max_v(*curr)),
                                "Cannot construct degenerated segment");
          // Construct a polyline with one segment.
          return Curve_2 (begin,end);
        }

        while (next != end) {
          CGAL_precondition_msg(!equal(get_min_v(*curr),get_max_v(*curr)),
                                "Cannot construct degenerated segment");

          // Verify that the segments' ends match and well-oriented
          CGAL_precondition_code(
            Point_2 curr_target;
            Point_2 next_source;

            if (cmp_seg_endpts(*curr) == SMALLER)
              curr_target = get_max_v(*curr);
            else
              curr_target = get_min_v(*curr);

            if (cmp_seg_endpts(*next) == SMALLER)
              next_source = get_min_v(*next);
            else
              next_source = get_max_v(*next);

            CGAL_precondition_msg(
              equal(curr_target,next_source),
              "Input must form a continuous, well-oriented polyline");
                                 );
          ++next;
          ++curr;
        }
        // Verify that the last segment is not degenerated
        CGAL_precondition_msg(!equal(get_min_v(*curr), get_max_v(*curr)),
                              "Cannot construct degenerated segment");

        return Curve_2(begin, end);
      }
    };

    /*! Get a Construct_curve_2 functor object. */
    Construct_curve_2 construct_curve_2_object() const
    { return Construct_curve_2(*this); }

    class Construct_x_monotone_curve_2 {
    protected:
      typedef Arr_polyline_traits_2<Segment_traits_2>       Geometry_traits_2;
      /*! The polyline traits (in case it has state) */
      const Geometry_traits_2& m_poly_traits;

    public:
      /*! Constructor. */
      Construct_x_monotone_curve_2(const Geometry_traits_2& traits) :
        m_poly_traits(traits) {}

      /*! Returns an x-monotone polyline connecting the two given endpoints.
       * \param p The first point.
       * \param q The second point.
       * \pre p and q must not be the same.
       * \return A segment connecting p and q.
       */
      X_monotone_curve_2 operator()(const Point_2& p, const Point_2& q) const
      {
        CGAL_precondition_code(
                               typename Segment_traits_2::Equal_2 equal =
                               m_poly_traits.segment_traits_2()->
                               equal_2_object();
                               );
        CGAL_precondition_msg(
          !equal(p,q),
          "Cannot construct a degenerated segment as a polyline");
        X_monotone_segment_2 seg(p, q);
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

        CGAL_precondition_code(
          const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
          // Initialize two comparison functors
          typename Segment_traits_2::Compare_x_2 compare_x =
          seg_traits->compare_x_2_object();
          typename Segment_traits_2::Compare_xy_2 compare_xy =
          seg_traits->compare_xy_2_object();
          // Make sure there is no changed of directions.
          // Save the comp_x between the first two points
          const Comparison_result cmp_x_res = compare_x(*ps, *pt);
          // Save the comp_xy between the first two points
          const Comparison_result cmp_xy_res = compare_xy(*ps, *pt);
                               );

        // Assure that the first two points are not the same.
        // Note that this also assures that no to consecutive points are equal
        // in the whole range.
        CGAL_precondition(cmp_xy_res != EQUAL);

        while (pt != end) {
          CGAL_precondition(compare_xy(*ps, *pt) == cmp_xy_res);
          CGAL_precondition(compare_x(*ps, *pt) == cmp_x_res);

          segs.push_back(X_monotone_segment_2(*ps,*pt));
          ++ps; ++pt;
        }

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
        // std::cout << "Construct x-poly from the following range:\n";
        // for (ForwardIterator it = begin; it !=end; ++it)
        //   std::cout << *it << std::endl;

        CGAL_precondition_msg(begin != end,
                              "Input range of segments has to contain at least"
                              "one segment");

        CGAL_precondition_code(
          const Segment_traits_2* seg_traits = m_poly_traits.segment_traits_2();
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

          if (next == end){
            CGAL_precondition_msg(!equal(get_max_v(*curr),get_min_v(*curr)),
                                  "Cannot construct a polyline "
                                  "with degenerated segment");
            return X_monotone_curve_2(begin,end);
          }

          // Range contains at least two segments

          Comparison_result init_dir = cmp_seg_endpts(*curr);
          while (next != end){
            CGAL_precondition_msg(!equal(get_min_v(*next),get_max_v(*next)),
                                  "Cannot construct a polyline "
                                  "with degenerated segment");
            CGAL_precondition_msg(init_dir == cmp_seg_endpts(*next),
                                  "Segments must form x-monotone polyline");
            if (init_dir == SMALLER){
              CGAL_precondition_msg(equal(get_max_v(*curr),get_min_v(*next)),
                                    "Segments should concatenate in "
                                    "source->target manner");
            }
            else{
              CGAL_precondition_msg(equal(get_min_v(*curr),get_max_v(*next)),
                                    "Segments should concatenate in "
                                    "source->target manner");
            }
            ++curr;
            ++next;
          }
                               );

        return X_monotone_curve_2(begin, end);
      }
    };

    /*! Get a Construct_x_monotone_curve_2 functor object. */
    Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const
    { return Construct_x_monotone_curve_2(*this); }
    //@}

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
    unsigned int locate(const X_monotone_curve_2& cv,
                        const Point_2& q) const
    {
      typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
        segment_traits_2()->construct_min_vertex_2_object();
      typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
        segment_traits_2()->construct_max_vertex_2_object();
      typename Segment_traits_2::Compare_x_2 compare_x =
        segment_traits_2()->compare_x_2_object();


      // The direction of cv. SMALLER means left-to-right and
      // otherwise right-to-left
      Comparison_result direction = segment_traits_2()->
        compare_endpoints_xy_2_object()(cv[0]);
      unsigned int from, to;
      if (direction == SMALLER) {
        from = 0;
        to = cv.number_of_segments() - 1;
      }
      else {
        from = cv.number_of_segments() - 1;
        to = 0;
      }

      // Test the case that cv is vertical
      if (segment_traits_2()->is_vertical_2_object()(cv[0])) {
        typename Segment_traits_2::Compare_xy_2 compare_xy =
          segment_traits_2()->compare_xy_2_object();

        // Test if q is one of cv's end points
        Comparison_result res_from;
        res_from = compare_xy(min_vertex(cv[from]), q);
        if (res_from == EQUAL) return from;

        Comparison_result res_to;
        res_to = compare_xy(max_vertex(cv[to]), q);
        if (res_to == EQUAL) return to;

        // Verify that q has the same x-coord as cv (which is vertical)
        Comparison_result res = compare_x(max_vertex(cv[to]), q);
        if ((res != EQUAL) || (res_to == res_from)) return INVALID_INDEX;

        // Perform a binary search to locate the segment that contains q in its
        // xy-range:
        while (((direction == SMALLER) && (to > from)) ||
               ((direction == LARGER) && (to < from))) {
          unsigned int mid = (from + to) / 2;
          if (((direction == SMALLER) && (mid > from)) ||
              ((direction == LARGER) && (mid < from))) {
            Comparison_result res_mid = compare_xy(min_vertex(cv[mid]), q);
            if (res_mid == EQUAL) return mid;
            if (res_mid == res_from) from = mid;
            else if (direction == SMALLER)
              to = mid - 1;
            else
              to = mid + 1;
          }
          else {
            CGAL_assertion(((direction == SMALLER) && (mid < to)) ||
                           ((direction == LARGER) && (mid > to)));
            Comparison_result res_mid = compare_xy(max_vertex(cv[mid]), q);
            if (res_mid == EQUAL) return mid;
            if (res_mid == res_to) to = mid;
            else if (direction == SMALLER)
              from = mid + 1;
            else
              from = mid - 1;
          }
        }
        // In case (from == to), and we know that the polyline contains the q:
        CGAL_assertion(from == to);
        return from;
      }

      // First check whether the polyline curve really contains q in
      // its x-range.
      Comparison_result res_from = compare_x(min_vertex(cv[from]), q);
      if (res_from == EQUAL) return from;

      Comparison_result res_to = compare_x(max_vertex(cv[to]), q);
      if (res_to == EQUAL) return to;

      // q is not in the x-range of cv:
      if (res_from == res_to) return INVALID_INDEX;

      // Perform a binary search and locate the segment that contains q in its
      // x-range:
      while (((direction == SMALLER) && (to > from)) ||
             ((direction == LARGER)  && (to < from))) {
        unsigned int mid = (from + to) / 2;

        if (((direction == SMALLER) && (mid > from)) ||
            ((direction == LARGER)  && (mid < from))) {
          Comparison_result res_mid = compare_x(min_vertex(cv[mid]), q);
          if (res_mid == EQUAL) {
            // Ensure that the returned segment contains the query point
            // on its right end (if possible)
            if ((direction == SMALLER) && (mid > 0))
              --mid;
            else if ((direction == LARGER) &&
                     ((mid + 1) < cv.number_of_segments()))
              ++mid;
            return mid;
          }
          if (res_mid == res_from) from = mid;
          else if (direction == SMALLER)
            to = mid - 1;
          else
            to = mid + 1;
        }
        else {
          CGAL_assertion(((direction == SMALLER) && (mid < to)) ||
                         ((direction == LARGER) && (mid > to)));
          Comparison_result res_mid = compare_x(max_vertex(cv[mid]), q);
          if (res_mid == EQUAL) return mid;
          if (res_mid == res_to) to = mid;
          else if (direction == SMALLER)
            from = mid + 1;
          else
            from = mid - 1;
        }
      }

      // In case(from == to), and we know that the polyline contains the q:
      CGAL_assertion(from == to);
      return from;
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
    unsigned int locate_side(const X_monotone_curve_2& cv,
                             const Point_2& q, const bool& to_right) const
    {
      // First locate a segment segments[i] that contains q in its x-range.
      unsigned int i = locate(cv, q);
      if (i == INVALID_INDEX)
        return INVALID_INDEX;

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
