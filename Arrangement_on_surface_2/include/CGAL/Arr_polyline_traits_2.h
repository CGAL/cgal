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
 * TODO: What to do with all the calls to push_back(seg)? Probably I should
 *       use the functor Push_Back_2 from this traits class or make sure that
 *       they maintain the well-oriented invariant.
 * TODO: Model the concept ArrangementDirectionalXMonotoneTraits_2.
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
#include <CGAL/Arr_geometry_traits/Polyline_2.h>

namespace CGAL {

template <class SegmentTraits_2>
class Arr_polyline_traits_2 {
public:
  typedef SegmentTraits_2                          Segment_traits_2;

  // Tag definitions:
  typedef Tag_true                                   Has_left_category;
  typedef Tag_true                                   Has_merge_category;
  /*
   * Roadmap:
   *       Add a Do_intersect_2 functor which stops after the first
   *       intersection found.  Either copy/paste the code of the
   *       Intersect_2 only without using output iterator and stopping
   *       after the first intersection.  Or, use a dummy output
   *       iterator which throws an exception as soon as it is
   *       changed. Then from Do_intersect_2 call the Intersect_2 and
   *       provide this dummy output iterator as the third argument.
   *       Note that this tag is false also for Arr_segment_traits_2.h
   */
  typedef Tag_false                                  Has_do_intersect_category;

  typedef Arr_oblivious_side_tag                     Left_side_category;
  typedef Arr_oblivious_side_tag                     Bottom_side_category;
  typedef Arr_oblivious_side_tag                     Top_side_category;
  typedef Arr_oblivious_side_tag                     Right_side_category;

private:
  typedef Arr_polyline_traits_2<Segment_traits_2>    Self;

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
  typedef typename Segment_traits_2::Curve_2            Segment_2;

  /*!
   * A polyline represents a general continuous piecewise-linear curve, without
   * degenerated segments.
   */
  typedef POLYLINE::Polyline_2<Segment_traits_2>            Curve_2;
  /*!
   * An x monotone polyline represents a continuous piecewise-linear curve which
   * is either strongly x-monotone or vertical. Again, the polyline is without
   * degenerated segments.
   */
  typedef POLYLINE::X_monotone_polyline_2<Segment_traits_2> X_monotone_curve_2;

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
    const Geometry_traits_2* m_poly_traits;

  public:
    /* Constructor. */
    Number_of_points_2(const Geometry_traits_2* traits) :
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
    const Geometry_traits_2* m_poly_traits;

  public:
    /* Constructor. */
    Construct_min_vertex_2(const Geometry_traits_2* traits) :
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
      return m_poly_traits->
        segment_traits_2()->construct_min_vertex_2_object()(cv[0]);
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(this); }

  class Construct_max_vertex_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>     Geometry_traits_2;
    /*! The polyline traits (in case it has state) */
    const Geometry_traits_2* m_poly_traits;

  public:
    /*! Constructor. */
    Construct_max_vertex_2(const Geometry_traits_2* traits) :
      m_poly_traits(traits)
    {}

    /*!
     * Get the right endpoint of the x-monotone curve(segment).
     * \param cv The polylinecurve.
     * \return The right endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2& cv) const
    {
      CGAL_assertion(cv.number_of_segments() > 0);
      return m_poly_traits->segment_traits_2()->
        construct_max_vertex_2_object()(cv[cv.number_of_segments() - 1]);
    }
  };

  /*! Get a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(this); }

  class Is_vertical_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>     Geometry_traits_2;
    /*! The polyline traits (in case it has state) */
    const Geometry_traits_2* m_poly_traits;

  public:
    /*! Constructor. */
    Is_vertical_2(const Geometry_traits_2* traits) : m_poly_traits(traits) {}

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
      return (m_poly_traits->segment_traits_2()->is_vertical_2_object()(cv[0]));
    }
  };

  /*! Get an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object() const
  { return Is_vertical_2(this); }

  class Compare_y_at_x_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>       Geometry_traits_2;
    /*! The polyline traits (in case it has state) */
    const Geometry_traits_2* m_poly_traits;

  public:
    /*! Constructor. */
    Compare_y_at_x_2(const Geometry_traits_2* traits) : m_poly_traits(traits) {}

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
        m_poly_traits->locate(cv, p);
      CGAL_precondition(i != INVALID_INDEX);

      // Compare the segment cv[i] and p.
      return m_poly_traits->segment_traits_2()->
        compare_y_at_x_2_object()(p, cv[i]);
    }
  };

  /*! Get a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(this); }

  class Compare_y_at_x_left_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>       Geometry_traits_2;
    /*! The polyline traits (in case it has state) */
    const Geometry_traits_2* m_poly_traits;

  public:
    /*! Constructor. */
    Compare_y_at_x_left_2(const Geometry_traits_2* traits) :
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
      unsigned int i1=m_poly_traits->locate_side(cv1, p, false);
      unsigned int i2=m_poly_traits->locate_side(cv2, p, false);

      CGAL_precondition(i1 != INVALID_INDEX);
      CGAL_precondition(i2 != INVALID_INDEX);

      // Compare cv1[i1] and cv2[i2] at p's left.
      return m_poly_traits->segment_traits_2()->
        compare_y_at_x_left_2_object()(cv1[i1], cv2[i2], p);
    }
  };

  /*! Get a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(this); }

  class Compare_y_at_x_right_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>       Geometry_traits_2;
    /*! The polyline traits (in case it has state) */
    const Geometry_traits_2* m_poly_traits;

  public:
    /*! Constructor. */
    Compare_y_at_x_right_2(const Geometry_traits_2* traits) :
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
      unsigned int i1=m_poly_traits->locate_side(cv1, p, true);
      unsigned int i2=m_poly_traits->locate_side(cv2, p, true);

      CGAL_precondition(i1 != INVALID_INDEX);
      CGAL_precondition(i2 != INVALID_INDEX);

      // Compare cv1[i1] and cv2[i2] at p's right.
      return m_poly_traits->segment_traits_2()->
        compare_y_at_x_right_2_object()(cv1[i1],cv2[i2], p);
    }
  };

  /*! Get a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(this); }

  class Equal_2 {
  protected:

    typedef Arr_polyline_traits_2<Segment_traits_2>     Geometry_traits_2;
    /*! The polyline traits (in case it has state) */
    const Geometry_traits_2* m_poly_traits;

  public:

    /*! Constructor. */
    Equal_2(const Geometry_traits_2* poly_tr) : m_poly_traits(poly_tr) {}

    /*!
     * Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same;(false) otherwise.
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const
    { return m_poly_traits->segment_traits_2()->equal_2_object()(p1, p2); }

    /*!
     * Check if the two x-monotone curves are the same(have the same graph).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return(true) if the two curves are the same;(false) otherwise.
     */
    bool operator()(const X_monotone_curve_2 & cv1,
                    const X_monotone_curve_2 & cv2) const
    {
      unsigned int n1 = cv1.number_of_segments();
      unsigned int n2 = cv2.number_of_segments();

      // Check the pairwise equality of the contained segments.
      const Segment_traits_2* seg_traits = m_poly_traits->segment_traits_2();
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
      Is_vertical_2 is_vertical = m_poly_traits->is_vertical_2_object();
      Point_2 point1,point2;
      Comparison_result res_x;
      Comparison_result res_y_at_x;
      unsigned int i = 0, j = 0;

      // the first and last points of the segments should be equal.
      bool res = equal(min_vertex(cv1[0]),min_vertex(cv2[0]));
      if (!res)
        return false;
      res = equal(max_vertex(cv1[n1-1]),max_vertex(cv2[n2-1]));
      if (!res)
        return false;

      // If the first and last points are equal and the curves are vertical,
      // it means that it is equal.
      bool ver1 = is_vertical(cv1);
      bool ver2 = is_vertical(cv2);
      // both curves are vertical and therefore equal.
      if (ver1 && ver2)
        return true;
      // one is vertical and the other is not - hence not equal.
      if ((ver1 && !ver2) || (ver2 && !ver1))
        return false;

      // If we arrived here it means that the first and last point of the
      // curve are equal.
      while ((i < n1-1) || (j < n2-1)) {
        point1 = max_vertex(cv1[i]);
        point2 = max_vertex(cv2[j]);

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
            res_y_at_x = compare_y_at_x(point1,cv2[j]);
            if (res_y_at_x == EQUAL)
              ++i;
            else
              return false;
          }
          else if(res_x == LARGER) {
            res_y_at_x = compare_y_at_x(point2,cv1[i]);
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
  { return Equal_2(this); }
  ///@}

  /// \name Construction functors(based on the segment traits).
  //@{

  class Make_x_monotone_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>     Geometry_traits_2;
    /*! The traits (in case it has state) */
    const Geometry_traits_2* m_poly_traits;

  public:
    /*! Constructor. */
    Make_x_monotone_2(const Geometry_traits_2* traits) :
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
      const_seg_iterator start_seg = cv.begin_segments();
      const_seg_iterator end_seg = cv.end_segments();

      // Empty polyline:
      if (start_seg == end_seg)
        return oi;

      const_seg_iterator it_next = start_seg;
      ++it_next;

      Construct_x_monotone_curve_2 construct_x_monotone_curve =
        m_poly_traits->construct_x_monotone_curve_2_object();

      if (it_next == end_seg) {
        // The polyline contains a single segment:
        *oi++ = make_object(construct_x_monotone_curve(*start_seg));
        return oi;
      }

      // Polyline contains at least 2 segments!

      const Segment_traits_2* seg_traits = m_poly_traits->segment_traits_2();
      // Functor defined in ArrangementDirectionalXMonotoneTraits_2
      typename Segment_traits_2::Compare_endpoints_xy_2 comp_endpts_xy =
        seg_traits->compare_endpoints_xy_2_object();
      typename Segment_traits_2::Is_vertical_2 is_vertical =
        seg_traits->is_vertical_2_object();

      const_seg_iterator it_start = start_seg;
      const_seg_iterator it_curr = start_seg;

      bool is_start_vertical = is_vertical(*it_start);
      Comparison_result start_dir = comp_endpts_xy(*it_start);

      for (/*it_next was advanced earlier*/; it_next != end_seg; ++it_next)
      {
        if (
            comp_endpts_xy(*it_next) != start_dir ||
            is_vertical(*it_next) != is_start_vertical
            )
          {
            // Construct an x-monotone curve from the sub-range which
            // was found
            *oi++ =
              make_object(construct_x_monotone_curve(it_start, it_next));
            it_start = it_next;
            is_start_vertical = is_vertical(*it_start);
            start_dir = comp_endpts_xy(*it_start);
          }
        it_curr = it_next;
      }
      *oi++ = make_object(construct_x_monotone_curve(it_start, it_next));
      return oi;
    }
  };

  /*! Get a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(this); }

  /*! Functor to enable pushing back of either points or segments to an
   *  existing polyline.
   */
  class Push_back_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>     Geometry_traits_2;
    /*! The traits (in case it has state) */
    const Geometry_traits_2* m_poly_traits;

  public:
    /*! Constructor. */
    Push_back_2(const Geometry_traits_2* traits) : m_poly_traits(traits) {}

    /*!
     * Append a point p to an existing polyline cv.
     * \param cv a polyline. Note, cv is not (necessarily) x-monotone.
     * \param p a point to be appended to cv
     * \pre cv contains at least one segment
     */
    void operator()(Curve_2& cv, const Point_2& p) const
    {
      int num_seg = cv.number_of_segments();
      CGAL_precondition(num_seg >0);
      int last_seg = num_seg-1;

      const Segment_traits_2* seg_traits = m_poly_traits->segment_traits_2();
      typename Segment_traits_2::Compare_endpoints_xy_2 comp_endpts=
        seg_traits->compare_endpoints_xy_2_object();

      /*
       * Since we assume that the segments of cv are well oriented,
       * pushing a single point to an existing polyline means that we
       * have to append the segment [cv[last_seg].target(),p]. The
       * following test determines which end of the last segment is
       * the target.
       */
      if (comp_endpts(cv[last_seg]) == SMALLER)
        {
          typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
            seg_traits->construct_max_vertex_2_object();
          cv.push_back(Segment_2(get_max_v(cv[last_seg]),p));
        }
      else
        {
          typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
            seg_traits->construct_min_vertex_2_object();
          cv.push_back(Segment_2(get_min_v(cv[last_seg]),p));
        }
    }


    /*!
     * Append a segment seg to an existing polyline cv. If cv is empty, seg will
     * be its first segment.
     * \param cv a polyline. Note, cv is not (necessarily) x-monotone.
     * \param seg a segment to be appended to cv
     * \pre One of the ends of seg should be equal to the last vertex of cv
     * \post The resulting cv is well oriented.
     */
    void operator()(Curve_2& cv, const Segment_2& seg) const
    {
      int num_seg = cv.number_of_segments();

      const Segment_traits_2* seg_traits = m_poly_traits->segment_traits_2();
      typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
        seg_traits->construct_min_vertex_2_object();
      typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
        seg_traits->construct_max_vertex_2_object();
      typename Segment_traits_2::Compare_endpoints_xy_2 comp_endpts =
        seg_traits->compare_endpoints_xy_2_object();
      typename Segment_traits_2::Equal_2 equal =
        seg_traits->equal_2_object();

      // cv is empty
      if (num_seg == 0)
        {
          cv.push_back(seg);
          return;
        }

      Point_2 last_v;

      /*
       * ROADMAP: In order to ensure that the resulting cv is well
       *          oriented the following test must(?) be made. If, in
       *          the future we allow ill-oriented polylines, then at
       *          this point seg can simply be appended.
       */
      if (comp_endpts(cv[num_seg-1]) == SMALLER)
        last_v = get_max_v(cv[num_seg-1]);
      else
        last_v = get_min_v(cv[num_seg-1]);

      CGAL_precondition(equal(last_v,get_min_v(seg))||
                        equal(last_v,get_max_v(seg)));

      if (equal(last_v,get_min_v(seg)))
        {
          if (comp_endpts(seg) == SMALLER)
            cv.push_back(seg);
          else
            cv.push_back(seg_traits->construct_opposite_2_object()(seg));
        }
      else
        {
          if (comp_endpts(seg) == SMALLER)
            cv.push_back(seg_traits->construct_opposite_2_object()(seg));
          else
            cv.push_back(seg);
        }
    }

    /*!
     * Append a point p to an existing x-monotone polyline cv.
     * \param cv the existing x-monotone polyline
     * \param p the point to be pushed back.
     * \pre cv contains at least one segment
     * \pre p is to the right of cv
     */
    void operator()(const X_monotone_curve_2& xcv, Point_2& p) const
    {
      int num_seg = xcv.number_of_segments();

      CGAL_precondition(num_seg > 0);

      const Segment_traits_2* seg_traits = m_poly_traits->segment_traits_2();
      typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
        seg_traits->construct_max_vertex_2_object();

      CGAL_precondition_code(
        typename Segment_traits_2::Compare_x_2 comp_x =
          seg_traits->compare_x_2_object();
        CGAL_precondition(comp_x(get_max_v(xcv[num_seg-1]),p)==LARGER);
                             );
      xcv.push_back(Segment_2(get_max_v(xcv[num_seg-1]),p));
    }

    /*!
     * Append a segment seg to an existing x-monotone polyline cv.
     * \param cv existing x-monotone polyline
     * \param seg the segment to be added
     * \pre cv is not an isolated point (in case it contains only one segment)
     * \pre seg is (strongly) right to cv, that it extends cv in
     *      a strong x-monotone manner.
     */
    void operator()(const X_monotone_curve_2& cv, Segment_2& seg) const
    {
      int num_seg = cv.number_of_segments();

      if (num_seg == 0)
        {
          cv.push_back(seg);
          return;
        }

      const Segment_traits_2* seg_traits = m_poly_traits->segment_traits_2();
      typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
        seg_traits->construct_max_vertex_2_object();
      typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
        seg_traits->construct_min_vertex_2_object();
      typename Segment_traits_2::Compare_x_2 comp_x =
        seg_traits->compare_x_2_object();
      typename Segment_traits_2::Equal_2 equal =
        seg_traits->equal_2_object();

      CGAL_precondition_code(
        if (num_seg == 1);
        CGAL_precondition(!equal(get_min_v(cv[0]),get_max_v(cv[0])));
                             );

      CGAL_precondition(equal(get_max_v(cv[num_seg-1]),get_min_v(seg)));
      CGAL_precondition(comp_x(get_min_v(seg),get_max_v(seg))==LARGER);

      cv.push_back(seg);
    }
  };

  /*! Get a Push_Back_2 functor object. */
  Push_back_2 push_back_2_object() const
  { return Push_back_2(this); }

  class Split_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>       Geometry_traits_2;
    /*! The polyline traits (in case it has state) */
    const Geometry_traits_2* m_poly_traits;

  public:
    /*! Constructor. */
    Split_2(const Geometry_traits_2* traits) : m_poly_traits(traits) {}

  public:
    /*!
     * Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve(p is its right endpoint).
     * \param c2 Output: The right resulting subcurve(p is its left endpoint).
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator()(const X_monotone_curve_2& cv, const Point_2& p,
                    X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      const Segment_traits_2* seg_traits = m_poly_traits->segment_traits_2();
      typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
        seg_traits->construct_min_vertex_2_object();
      typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
        seg_traits->construct_max_vertex_2_object();
      typename Segment_traits_2::Equal_2 equal =
        seg_traits->equal_2_object();

      // Make sure the split point is not one of the curve endpoints.
      CGAL_precondition(!equal(min_vertex(cv[0]), p));
      CGAL_precondition(!equal(max_vertex(cv[cv.number_of_segments() - 1]), p));

      // Locate the segment on the polyline cv that contains p.
      unsigned int i = m_poly_traits->locate(cv, p);
      CGAL_precondition(i != INVALID_INDEX);

      // Clear the output curves.
      c1.clear();
      c2.clear();

      // Push all segments labeled(0, 1, ... , i-1) into c1.
      for (int j = 0; j < i; ++j)
        c1.push_back(cv[j]);

      // Check whether the split point is cv[i]'s source or target.
      if (equal(max_vertex(cv[i]), p)) {
        // The entire i'th segment belongs to c1:
        c1.push_back(cv[i]);
      } else if (equal(min_vertex(cv[i]), p)) {
        // The entire i'th segments belongs to c2:
        c2.push_back(cv[i]);
      } else {
        // The i'th segment should be split: The left part(seg1) goes to cv1,
        // and the right part(seg2) goes to cv2.
        Segment_2   seg1, seg2;
        m_poly_traits->segment_traits_2()->
          split_2_object()(cv[i], p, seg1, seg2);

        c1.push_back(seg1);
        c2.push_back(seg2);
      }

      // Push all segments labeled(i+1, i+2, ... , n-1) into cv1.
      unsigned int n = cv.number_of_segments();

      for (int j = i+1; j < n; ++j)
        c2.push_back(cv[j]);
    }
  };

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object() const
  { return Split_2(this); }

  class Intersect_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>       Geometry_traits_2;
    /*! The polyline traits (in case it has state) */
    const Geometry_traits_2* m_poly_traits;

  public:
    /*! Constructor. */
    Intersect_2(const Geometry_traits_2* traits) : m_poly_traits(traits) {}

    /*!
     * Find the intersections of the two given curves and insert them into the
     * given output iterator. As two segments may itersect only once, only a
     * single intersection will be contained in the iterator.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template<class OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              OutputIterator oi)
    {
      const Segment_traits_2* seg_traits = m_poly_traits->segment_traits_2();
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

      const unsigned int n1 = cv1.number_of_segments();
      const unsigned int n2 = cv2.number_of_segments();
      unsigned int       i1 = 0;
      unsigned int       i2 = 0;
      X_monotone_curve_2 ocv;           // Used to represent overlaps.

      Comparison_result left_res =
        compare_xy(min_vertex(cv1[i1]), min_vertex(cv2[i2]));

      if (left_res == SMALLER) {
        // cv1's left endpoint is to the left of cv2's left endpoint:
        // Locate the index i1 of the segment in cv1 which contains cv2's
        // left endpoint.
        i1 = m_poly_traits->locate(cv1, min_vertex(cv2[i2]));
        if (i1 == INVALID_INDEX)
          return oi;

        if (equal(max_vertex(cv1[i1]), min_vertex(cv2[i2]))) {
          ++i1;
          left_res = EQUAL;
        }
      }  else if (left_res == LARGER) {
        // cv1's left endpoint is to the right of cv2's left endpoint:
        // Locate the index i2 of the segment in cv2 which contains cv1's
        // left endpoint.
        i2 = m_poly_traits->locate(cv2, min_vertex(cv1[i1]));

        if (i2 == INVALID_INDEX)
          return oi;

        if (equal(max_vertex(cv2[i2]), min_vertex(cv1[i1]))) {
          ++i2;
          left_res = EQUAL;
        }
      }

      // Check if the the left endpoint lies on the other polyline.
      bool left_coincides = (left_res == EQUAL);
      bool left_overlap = false;

      if (left_res == SMALLER)
        left_coincides =(compare_y_at_x(min_vertex(cv2[i2]), cv1[i1]) == EQUAL);
      else if (left_res == LARGER)
        left_coincides =(compare_y_at_x(min_vertex(cv1[i1]), cv2[i2]) == EQUAL);

      // The main loop: Go simultaneously over both polylines.
      Comparison_result right_res = left_res;
      bool right_coincides = left_coincides;
      bool right_overlap = false;

      while (i1 < n1 && i2 < n2) {
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
        } else if (right_coincides && left_coincides) {
          // An overlap exists between the current segments of the polylines:
          // Output the overlapping segment.
          right_overlap = true;
          if (left_res == SMALLER) {
            if (right_res == SMALLER) {
              Segment_2 seg(min_vertex(cv2[i2]), max_vertex(cv1[i1]));
              ocv.push_back(seg);
            } else {
              Segment_2 seg(min_vertex(cv2[i2]), max_vertex(cv2[i2]));
              ocv.push_back(seg);
            }
          } else {
            if (right_res == SMALLER) {
              Segment_2 seg(min_vertex(cv1[i1]), max_vertex(cv1[i1]));
              ocv.push_back(seg);
            } else {
              Segment_2 seg(min_vertex(cv1[i1]), max_vertex(cv2[i2]));
              ocv.push_back(seg);
            }
          }
        } else if (left_coincides && !right_coincides) {
          // The left point of the current segment of one polyline
          // coincides with the current segment of the other polyline.
          if (left_overlap) {
            // An overlap occured at the previous iteration:
            // Output the overlapping polyline.
            CGAL_assertion(ocv.number_of_segments() > 0);
            *oi++ = make_object(ocv);
            ocv.clear();
          } else {
            // The left point of the current segment of one polyline
            // coincides with the current segment of the other polyline, and
            // no overlap occured at the previous iteration:
            // Output the intersection point. The derivative of at least one of
            // the polylines is not defined at this point, so we give it
            // multiplicity 0.
            if (left_res == SMALLER) {
              std::pair<Point_2,Multiplicity>  p(min_vertex(cv2[i2]), 0);
              *oi++ = make_object(p);
            } else {
              std::pair<Point_2,Multiplicity>  p(min_vertex(cv1[i1]), 0);
              *oi++ = make_object(p);
            }
          }
        }

        // Proceed forward.
        if (right_res != SMALLER)
          ++i2;
        if (right_res != LARGER)
          ++i1;

        left_res = (right_res == SMALLER) ? LARGER :
          (right_res == LARGER) ? SMALLER : EQUAL;

        left_coincides = right_coincides;
        left_overlap = right_overlap;
      }

#if 0
      std::cout << "left coincides: " << left_coincides << std::endl;
      std::cout << "right coincides: " << right_coincides << std::endl;
      std::cout << "right res: " << right_res << std::endl;
      std::cout << "left res: " << left_res << std::endl;
#endif

      // Output the remaining overlapping polyline, if necessary.
      if (ocv.number_of_segments() > 0) {
        *oi++ = make_object(ocv);
      } else if (right_coincides) {
        if (right_res == SMALLER) {
          std::pair<Point_2,Multiplicity> ip(max_vertex(cv1[i1 - 1]), 0);
          *oi++ = make_object(ip);
        } else if (right_res == LARGER) {
          std::pair<Point_2,Multiplicity> ip(max_vertex(cv2[i2 - 1]), 0);
          *oi++ = make_object(ip);
        } else if (i1 > 0) {
          std::pair<Point_2,Multiplicity> ip(max_vertex(cv1[i1 - 1]), 0);
          *oi++ = make_object(ip);
        } else {
          CGAL_assertion(i2 > 0);
          std::pair<Point_2,Multiplicity> ip(max_vertex(cv2[i2 - 1]), 0);
          *oi++ = make_object(ip);
        }
      }
      return oi;
    }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object() const
  { return Intersect_2(this); }

  class Are_mergeable_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>       Geometry_traits_2;
    /*! The polyline traits (in case it has state) */
    const Geometry_traits_2* m_poly_traits;

  public:
    /*! Constructor. */
    Are_mergeable_2(const Geometry_traits_2* traits) : m_poly_traits(traits) {}

    /*!
     * Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return(true) if the two curves are mergeable, that is, they share a
     * common endpoint;(false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const
    {
      const Segment_traits_2* seg_traits = m_poly_traits->segment_traits_2();
      typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
        seg_traits->construct_min_vertex_2_object();
      typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
        seg_traits->construct_max_vertex_2_object();
      typename Segment_traits_2::Equal_2 equal =
        seg_traits->equal_2_object();
      typename Segment_traits_2::Is_vertical_2 is_vertical =
        seg_traits->is_vertical_2_object();

      const unsigned int n1 = cv1.number_of_segments();
      const unsigned int n2 = cv2.number_of_segments();

      bool ver1 = is_vertical(cv1[0]);
      bool ver2 = is_vertical(cv2[0]);

      return ((equal(max_vertex(cv1[n1 - 1]), min_vertex(cv2[0])) ||
               equal(max_vertex(cv2[n2 - 1]), min_vertex(cv1[0]))) &&
              ((ver1 && ver2) || (!ver1 && !ver2)));
    }
  };

  /*! Get an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(this); }

  /*! \class Merge_2
   * A functor that merges two x-monotone arcs into one.
   */
  class Merge_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>     Geometry_traits;
    /*! The traits (in case it has state) */
    const Geometry_traits* m_poly_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Merge_2(const Geometry_traits* traits) : m_poly_traits(traits) {}

  public:
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
      CGAL_precondition(m_poly_traits->are_mergeable_2_object()(cv1, cv2));

      Construct_min_vertex_2 get_min_v =
        m_poly_traits->construct_min_vertex_2_object();
      Construct_max_vertex_2 get_max_v =
        m_poly_traits->construct_max_vertex_2_object();
      Equal_2 equal = m_poly_traits->equal_2_object();

      c.clear();
      if (equal(get_max_v(cv1), get_min_v(cv2))) {

        const unsigned int n1 = cv1.number_of_segments();
        const unsigned int n2 = cv2.number_of_segments();
        unsigned int       i;

        // cv2 extends cv1 to the right:
        for (i = 0; i < n1 - 1; ++i)
          c.push_back(cv1[i]);

        // Try to merge the to contiguous line segments:
        if (m_poly_traits->segment_traits_2()->
            are_mergeable_2_object()(cv1[n1 - 1], cv2[0])) {
          Segment_2       seg;
          m_poly_traits->segment_traits_2()->
            merge_2_object()(cv1[n1 - 1], cv2[0], seg);
          c.push_back(seg);
        } else {
          c.push_back(cv1[n1 - 1]);
          c.push_back(cv2[0]);
        }

        for (i = 1; i < n2; ++i)
          c.push_back(cv2[i]);
      } else
        return this->operator()(cv2,cv1,c);
    }
  };

  /*! Get a Merge_2 functor object. */
  Merge_2 merge_2_object() const { return Merge_2(this); }
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
    const Geometry_traits_2* m_poly_traits;

  public:
    /*! Constructor. */
    Construct_curve_2 (const Geometry_traits_2* traits) :
      m_poly_traits(traits) {}

    /*! Returns an polyline connecting the two given endpoints.
     * \param p The first point.
     * \param q The second point.
     * \pre p and q are distinct.
     * \return A segment connecting p and q.
     */
    Curve_2 operator()(const Point_2& p, const Point_2& q) const
    {
      CGAL_precondition_code
        (
         typename Segment_traits_2::Equal_2 equal =
         m_poly_traits->segment_traits_2()->equal_2_object();
         );
      CGAL_precondition_msg (!equal(p,q),
                             "Cannot construct a degenerated segment");
      return Curve_2(Segment_2(p,q));
    }

    /*! Returns an polyline consists of one given segment.
     * \param seg input segment
     * \return A polyline with one segment, namely seg.
     */
    Curve_2 operator()(const Segment_2& seg) const
    {
      CGAL_precondition_code
        (
         /*
          * Test that the segment is not degenerated. We do this test
          * independently from the SegmentTraits in use, as we do not allow
          * a polyline with degenerated segments.
          */
         const Segment_traits_2* seg_traits = m_poly_traits->segment_traits_2();
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

    /*! Construct a polyline from a range of objects.
     *  \param begin An iterator pointing to the first segment in the range.
     *  \param end An iterator pointing to the past-the-end segment in the range
     *  \return A polyline using the corresponding construction implementation.
     */
    template <typename ForwardIterator>
    Curve_2 operator()(ForwardIterator begin, ForwardIterator end) const
    {
      typedef typename std::iterator_traits<ForwardIterator>::value_type VT;
      typedef typename boost::is_same<VT,Point_2>::type Is_point;
      return constructor_impl(begin, end, Is_point());
    }

    /*! Construction of a polyline from a range of points.
     * \pre The range contains at least two points
     * \pre Consecutive points are disjoint.
     * \return Polyline connecting the points in the input using the same
     *         order in which they were given.
     */
    template <typename ForwardIterator>
    Curve_2 constructor_impl (ForwardIterator begin, ForwardIterator end,
                             boost::true_type) const
    {
      // Container of the segments to be created.
      std::vector<Segment_2> segs;

      // The range must contain at least two points.
      CGAL_precondition_msg (std::distance(begin,end)>1,
                             "Cannot construct a polyline from one point");
      CGAL_precondition_code
        (
         typename Segment_traits_2::Equal_2 equal =
         m_poly_traits->segment_traits_2()->equal_2_object();
         );
      ForwardIterator curr = begin;
      ForwardIterator next = curr;
      ++next;
      while (next!=end)
        {
          CGAL_precondition_msg(!equal(*curr,*next),
                                "Cannot construct a degenerated segment");
          segs.push_back(Segment_2(*curr,*next));
          ++next;
          ++curr;
        }

      return Curve_2(segs.begin(),segs.end());
    }

    /*! Construction implementation from a range of segments.
     *  \pre The segments form a polyline, that is one of the ends of the i-th
     *       segment should coincide with one of the ends of the (i+1)-th
     *       segment. NOTE: The segments are not assumed to have any
     *       orientation. The end points are extracted using the
     *       Construct_max_vertex_2 Construct_min_vertex_2 functors which are
     *       provided by the SegmentTraits class.
     */
    template <typename ForwardIterator>
    Curve_2 constructor_impl (ForwardIterator begin, ForwardIterator end,
                             boost::false_type) const
    {
      // Range has to contain at least one segment
      CGAL_precondition(begin != end);

      ForwardIterator curr = begin;
      ForwardIterator next = curr;

      CGAL_precondition_code
        (
         const Segment_traits_2* seg_traits = m_poly_traits->segment_traits_2();
         typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
         seg_traits->construct_min_vertex_2_object();
         typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
         seg_traits->construct_max_vertex_2_object();
         typename Segment_traits_2::Compare_xy_2 comp_xy =
         seg_traits->compare_xy_2_object();
         typename Segment_traits_2::Equal_2 equal =
         seg_traits->equal_2_object();
         );

      if (++next == end)
        {
          CGAL_precondition_msg (!equal(get_min_v(*curr),get_max_v(*curr)),
                                 "Cannot construct degenerated segment");
          // Construct a polyline with one segment.
          return Curve_2 (begin,end);
        }

      while (next != end)
        {
          CGAL_precondition_msg (!equal(get_min_v(*curr),get_max_v(*curr)),
                                 "Cannot construct degenerated segment");
          // Verify that the segments' ends match
          CGAL_precondition(
                       comp_xy (get_min_v(*curr),get_min_v(*next)) == EQUAL ||
                       comp_xy (get_min_v(*curr),get_max_v(*next)) == EQUAL ||
                       comp_xy (get_max_v(*curr),get_min_v(*next)) == EQUAL ||
                       comp_xy (get_max_v(*curr),get_max_v(*next)) == EQUAL );
          ++next;
          ++curr;
        }
      CGAL_precondition_msg (!equal(get_min_v(*curr),get_max_v(*curr)),
                             "Cannot construct degenerated segment");

      return Curve_2 (begin, end);
    }
  };

  /*! Get a Construct_curve_2 functor object. */
  Construct_curve_2 construct_curve_2_object() const
  { return Construct_curve_2(this); }

  class Construct_x_monotone_curve_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>       Geometry_traits_2;
    /*! The polyline traits (in case it has state) */
    const Geometry_traits_2* m_poly_traits;

  public:
    /*! Constructor. */
    Construct_x_monotone_curve_2(const Geometry_traits_2* traits) :
      m_poly_traits(traits)
    {}


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
                             m_poly_traits->segment_traits_2()->
                             equal_2_object();
                             );
      CGAL_precondition_msg(!equal(p,q),
                       "Cannot construct a degenerated segment as a polyline");
      typename Segment_traits_2::Compare_xy_2 comp_xy =
        m_poly_traits->segment_traits_2()->compare_xy_2_object();

      if (comp_xy(p,q) == SMALLER)
        return X_monotone_curve_2(Segment_2(p,q));
      else
        return X_monotone_curve_2(Segment_2(q,p));
    }

    /*! Returns an x-monotone polyline consists of one given segment.
     * \param seg input segment
     * \return An x-monotone polyline with one segment, namely seg.
     */
    X_monotone_curve_2 operator()(const Segment_2& seg) const
    {
      CGAL_precondition_code
        (
         /*
          * Test that the segment is not degenerated. We do this test
          * independently from the SegmentTraits in use, as we do not allow
          * a polyline with degenerated segments.
          */
         const Segment_traits_2* seg_traits = m_poly_traits->segment_traits_2();
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

    template <typename ForwardIterator>
    X_monotone_curve_2 operator()(ForwardIterator begin,
                                  ForwardIterator end) const
    {
      typedef typename std::iterator_traits<ForwardIterator>::value_type VT;
      typedef typename boost::is_same<VT,Point_2>::type Is_point;
      return constructor_impl(begin, end, Is_point());
    }

    /*
     * Construct a polyline from a range of points.
     * As the input is a range of points it is natural to orient the
     * resulting polyline from left to right.
     * \pre no two consecutive points are the same
     * \pre The points form an x-monotone polyline
     * \post The constructed polyline is x-monotone and
     *       oriented from left to right.
     */
    template <typename ForwardIterator>
    X_monotone_curve_2 constructor_impl(ForwardIterator begin,
                                        ForwardIterator end,
                                        boost::true_type) const
    {
      // Vector of the segments to be constructed from the range of points
      std::vector<Segment_2> segs;

      // Make sure the range of points contains at least two points.
      ForwardIterator ps = begin;
      CGAL_precondition (ps != end);
      ForwardIterator pt = ps;
      ++pt;
      CGAL_precondition (pt != end);

      const Segment_traits_2* seg_traits = m_poly_traits->segment_traits_2();

      // Initialize two comparison functors
      CGAL_precondition_code(typename Segment_traits_2::Compare_x_2 compare_x =
                             seg_traits->compare_x_2_object();
                             );
      typename Segment_traits_2::Compare_xy_2 compare_xy =
        seg_traits->compare_xy_2_object();

      // Make sure there is no change of directions as we traverse the polyline.
      // Save the comp_x between the first two points
      CGAL_precondition_code (const Comparison_result cmp_x_res =
                              compare_x(*ps, *pt););
      // Save the comp_xy between the first two points
      const Comparison_result cmp_xy_res = compare_xy(*ps, *pt);

      // Assure that the first two points are not the same.
      // Note that this also assures that no to consecutive points are equal
      // in the whole range.
      CGAL_precondition (cmp_xy_res != EQUAL);

      while (pt != end) {
        CGAL_precondition (compare_xy(*ps, *pt) == cmp_xy_res);
        CGAL_precondition (compare_x(*ps, *pt) == cmp_x_res);
        cmp_xy_res == LARGER ?
          segs.push_back(Segment_2(*pt,*ps)) :
          segs.push_back(Segment_2(*ps,*pt));
        ++ps; ++pt;
      }

      // Reverse the polyline so it always directed from left to right.
      if (cmp_xy_res == LARGER)
        {
          // The constructed polyline has to be reversed.
          return X_monotone_curve_2(segs.rbegin(),segs.rend());
        }
      else
        // No need to reverse. Returning the polyline
        return X_monotone_curve_2(segs.begin(), segs.end());
    }

    /*! Returns an x-monotone polyline from a range of segments.
     * \param begin An iterator pointing to the first segment in the range.
     * \param end An iterator pointing to the past-the-end segment in the range.
     * \pre The range contains at least one segment.
     * \pre One endpoint of the i-th segment is an endpoint of the
     *      (i+1)th segment.
     * \pre The sequence of segments in the range forms a weak x-monotone
     *      polyline.
     * \pre The container should support bidirectional iteration.
     * \post The resulting x-monotone polyline directed from left to
     *      right.
     * \return An x-monotone polyline directed from left to right.
     */
    template <typename ForwardIterator>
    X_monotone_curve_2 constructor_impl(ForwardIterator begin,
                                        ForwardIterator end,
                                        boost::false_type) const
    {
      CGAL_precondition_msg(begin != end,
                            "Input range of segments has to contain at least"
                            "one segment");

      // Functors that have to be used always
      const Segment_traits_2* seg_traits = m_poly_traits->segment_traits_2();
      typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
        seg_traits->construct_min_vertex_2_object();
      typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
        seg_traits->construct_max_vertex_2_object();
      typename Segment_traits_2::Equal_2 equal =
        seg_traits->equal_2_object();

      CGAL_precondition_code
        (
         // A functor which is used only when validity tests of the
         // input have to be ran.
         typename Segment_traits_2::Is_vertical_2 is_vertical =
           seg_traits->is_vertical_2_object();

         ForwardIterator curr = begin;
         // Ensure that the first segment does not degenerate to a point.
         CGAL_precondition_msg(!equal(get_min_v(*curr), get_max_v(*curr)),
                               "Cannot construct a degenerated segment");

         ForwardIterator next = curr;

         if (++next != end) {
           // Ensure that the second segment does not degenerate to a point.
           CGAL_precondition(!equal(get_min_v(*next), get_max_v(*next)));

           // Ensure that either both are vertical or both are not vertical.

           CGAL_precondition((is_vertical(*curr) && is_vertical(*next)) ||
                             (!is_vertical(*curr) && !is_vertical(*next)));

           // Ensure that the segment connect.
           CGAL_precondition(equal(get_max_v(*curr), get_min_v(*next)) ||
                             equal(get_min_v(*curr), get_max_v(*next)));

           // Record the initial direction.
           bool left_to_right = equal(get_max_v(*curr), get_min_v(*next));

           for (curr = next++; next != end; ++next) {
             // Ensure that the next segment does not degenerate to a point
             CGAL_precondition(!equal(get_min_v(*next), get_max_v(*next)));

             // Ensure that either both are vertical or both are not vertical.
             CGAL_precondition((is_vertical(*curr) && is_vertical(*next)) ||
                               (!is_vertical(*curr) && !is_vertical(*next)));

             // Ensure the direction and connectivity.
             CGAL_precondition((left_to_right &&
                                equal(get_max_v(*curr), get_min_v(*next))) ||
                               (!left_to_right &&
                                equal(get_max_v(*next), get_min_v(*curr))));
             ++curr;
           }
         }
         );

      bool rev = false;
      if (std::distance(begin, end) >= 2) {
        ForwardIterator second = begin;
        ++second;
        rev = equal(get_min_v(*begin), get_max_v(*second));
      }
      // The following statement assumes that the begin (and end) iterators
      // are biderctional.
      if (rev)
        {
          // Reverse the _order_ of the segments in the container.
          return X_monotone_curve_2(
                         std::reverse_iterator<ForwardIterator>(end),
                         std::reverse_iterator<ForwardIterator>(begin));
        }
      else
        return X_monotone_curve_2(begin, end);
    }
  };

  /*! Get a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const
  { return Construct_x_monotone_curve_2(this); }
  //@}

private:
  /*
   * Roadmap: Improve the implementation _locate()
   *  - _locate() should return an iterator to the located segment
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

    unsigned int from = 0;
    unsigned int to = cv.number_of_segments() - 1;

    if (segment_traits_2()->is_vertical_2_object()(cv[0])) {
      typename Segment_traits_2::Compare_xy_2 compare_xy =
        segment_traits_2()->compare_xy_2_object();

      // First check whether the polyline curve really contains q in its
      // xy-range:

      Comparison_result res_from = compare_xy(min_vertex(cv[from]), q);
      if (res_from == EQUAL) return from;

      Comparison_result res_to = compare_xy(max_vertex(cv[to]), q);
      if (res_to == EQUAL) return to;

      typename Segment_traits_2::Compare_x_2 compare_x =
        segment_traits_2()->compare_x_2_object();

      Comparison_result res = compare_x(max_vertex(cv[to]), q);
      if (res != EQUAL) return INVALID_INDEX;

      //// q is not in the x-range of cv:
      //if (res_from == res_to) return INVALID_INDEX;

      // Perform a binary search to locate the segment that contains q in its
      // xy-range:
      while (to > from) {
        unsigned int mid = (from + to) / 2;
        if (mid > from) {
          Comparison_result res_mid = compare_xy(min_vertex(cv[mid]), q);
          if (res_mid == EQUAL) return mid;
          if (res_mid == res_from) from = mid;
          else to = mid - 1;
        } else {
          CGAL_assertion(mid < to);
          Comparison_result res_mid = compare_xy(max_vertex(cv[mid]), q);
          if (res_mid == EQUAL) return mid;
          if (res_mid == res_to) to = mid;
          else from = mid + 1;
        }
      }
      // In case (from == to), and we know that the polyline contains the q:
      CGAL_assertion(from == to);
      return from;
    }

    typename Segment_traits_2::Compare_x_2 compare_x =
      segment_traits_2()->compare_x_2_object();

    // First check whether the polyline curve really contains q in its x-range.
    Comparison_result res_from = compare_x(min_vertex(cv[from]), q);
    if (res_from == EQUAL) return from;

    Comparison_result res_to = compare_x(max_vertex(cv[to]), q);
    if (res_to == EQUAL) return to;

    // q is not in the x-range of cv:
    if (res_from == res_to) return INVALID_INDEX;

    // Perform a binary search and locate the segment that contains q in its
    // x-range:
    while (to > from) {
      unsigned int mid = (from + to) / 2;

      if (mid > from) {
        Comparison_result res_mid = compare_x(min_vertex(cv[mid]), q);
        if (res_mid == EQUAL) return mid;
        if (res_mid == res_from) from = mid;
        else to = mid - 1;
      } else {
        CGAL_assertion(mid < to);
        Comparison_result res_mid = compare_x(max_vertex(cv[mid]), q);
        if (res_mid == EQUAL) return mid;
        if (res_mid == res_to) to = mid;
        else from = mid + 1;
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

    typename Segment_traits_2::Equal_2 equal = segment_traits_2()->
      equal_2_object();

    if (equal(segment_traits_2()->construct_min_vertex_2_object()(cv[i]), q)) {
      // q is the left endpoint of the i'th segment:
      if (to_right)
        return i;
      else if (i == 0)
        return INVALID_INDEX;
      else
        return i - 1;
    }

    if (equal(segment_traits_2()->construct_max_vertex_2_object()(cv[i]), q)) {
      // q is the right endpoint of the i'th segment:
      if (!to_right)
        return i;
      else if (i == (cv.number_of_segments() - 1))
        return INVALID_INDEX;
      else
        return i + 1;
    }

    // In case q is in cv[i]'s interior:
    return i;
  }
};


} //namespace CGAL

#endif
