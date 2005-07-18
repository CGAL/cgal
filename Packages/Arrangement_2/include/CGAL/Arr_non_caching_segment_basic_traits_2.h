// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Iddo Hanniel      <hanniel@math.tau.ac.il>
//                 Eyal Flato        <flato@post.tau.ac.il>
//                 Oren Nechushtan   <theoren@math.tau.ac.il>
//                 Eti Ezra          <estere@post.tau.ac.il>
//                 Shai Hirsch       <shaihi@post.tau.ac.il>
//                 Eugene Lipovetsky <eug@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
//                 Ron Wein          <wein@post.tau.ac.il>

#ifndef CGAL_ARR_NON_CACHING_SEGMENT_BASIC_TRAITS_2_H
#define CGAL_ARR_NON_CACHING_SEGMENT_BASIC_TRAITS_2_H

/*! \file The basic non-caching segment traits-class for the arrangement
 * package. This traits class handles x-monotone non-intersecting segments.
 * It is a model of the ArrangementBasicTraits_2 concept. The class is
 * templated by a kernel and inherits from it all the types and many of the
 * functors required by the concept it models.
 */

#include <CGAL/tags.h>
#include <CGAL/assertions.h>

CGAL_BEGIN_NAMESPACE

/*! A model of the ArrangementBasicTraits_2 concept that handles x-monotone
 * non-intersecting segments
 */
template <class T_Kernel>
class Arr_non_caching_segment_basic_traits_2 : public T_Kernel
{
public:
  typedef T_Kernel                              Kernel;

  // Categories:
  //#define HAS_LEFT_NOT
#if !defined(HAS_LEFT_NOT)
  typedef Tag_true                              Has_left_category;
#else
  typedef Tag_false                             Has_left_category;
#endif

  //#define HAS_REFLECT
#if !defined(HAS_REFLECT)
  typedef Tag_false                             Has_reflect_category;
#else
  typedef Tag_true                              Has_reflect_category;
#endif
    
  /*! Default Constructor */
  Arr_non_caching_segment_basic_traits_2() {}
  
  /// \name Types and functor inherited from the kernel
  //@{

  // Traits types:
  typedef typename Kernel::Point_2              Point_2;
  typedef typename Kernel::Segment_2            X_monotone_curve_2;

  /*! Compare the x-coordinates of two points */
  typedef typename Kernel::Compare_x_2                  Compare_x_2;

  /*! Compare two points lexigoraphically; by x, then by y */
  typedef typename Kernel::Compare_xy_2                 Compare_xy_2;

  // The following 2 should be added to the kernel and enabled accordingly
  /*! Obtain the left endpoint of a given segment */
  typedef typename Kernel::Construct_min_vertex_2       Construct_min_vertex_2;

  /*! Obtain the right endpoint of a given segment */
  typedef typename Kernel::Construct_max_vertex_2       Construct_max_vertex_2;
  
  /*! Check whether a given segment is vertical */
  typedef typename Kernel::Is_vertical_2                Is_vertical_2;
  
  /*! Return the location of a given point with respect to an input segment */
  typedef typename Kernel::Compare_y_at_x_2             Compare_y_at_x_2;
  
  /*! Check if two segments or if two points are identical */
  typedef typename Kernel::Equal_2                      Equal_2;

  ///@}

  /// \name Functor introduced here (based on the kernel)
  //@{

#if !defined(HAS_LEFT_NOT)

  /*! A functor for comparing two segments to the left of a point */
  class Compare_y_at_x_left_2 {
  public:
    /* Compare the y value of two segments immediately to the left of their
     * intersection point.
     * \param cv1 The first segment.
     * \param cv2 The second segment.
     * \param p The intersection point.
     * \pre The point p lies on both segments, and both of them must be also be
     *      defined (lexicographically) to its left.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the left of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator()(const X_monotone_curve_2 & cv1,
                                 const X_monotone_curve_2 & cv2,
                                 const Point_2 & p) const
    {
      Kernel kernel;

      // The two segments must be defined at q and also to its left.
      CGAL_precondition_code(
        Compare_y_at_x_2 compare_y_at_x = kernel.compare_y_at_x_2_object();
        );
      
      CGAL_precondition(compare_y_at_x(p, cv1) == EQUAL &&
                        compare_y_at_x(p, cv2) == EQUAL);

      CGAL_precondition_code(
        Compare_xy_2 compare_xy = kernel.compare_xy_2_object();
        typename Kernel::Construct_vertex_2 construct_vertex =
          kernel.construct_vertex_2_object();
        const Point_2 & source1 = construct_vertex(cv1, 0);
        const Point_2 & target1 = construct_vertex(cv1, 1);
        const Point_2 & left1 =
          (kernel.less_xy_2_object()(source1, target1)) ? source1 : target1;
        const Point_2 & source2 = construct_vertex(cv2, 0);
        const Point_2 & target2 = construct_vertex(cv2, 1);
        const Point_2 & left2 =
          (kernel.less_xy_2_object()(source2, target2)) ? source2 : target2;
        );
      
      CGAL_precondition(compare_xy(left1, p) == SMALLER &&
                        compare_xy(left2, p) == SMALLER);
    
      // Compare the slopes of the two segments to determine thir relative
      // position immediately to the left of q.
      // Notice that we swap the order of the curves in order to obtain the
      // correct result to the left of p.
      return kernel.compare_slope_2_object()(cv2, cv1);
    }
  };

  /*! Obtain a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  {
    return Compare_y_at_x_left_2();
  }
  
#else

  /*! Afunctor for reflecting a point or a segment through the origin */
  class Reflect_xy_2 {
  public:
    /*! Reflects the given point through the origin.
     * \param pt The point to be reflected.
     * \return The reflected point.
     */
    Point_2 operator()(const Point_2 & pt) const
    {
      Kernel kernel;
      Point_2 org(ORIGIN);      
      typename Kernel::Vector_2 v = kernel.construct_vector_2_object()(pt, org);
      Point_2 reflected_pt = org + v;
      return reflected_pt;
    }

    /*! Reflects the given segment through the origin.
     * \param cv The segment to be reflected.
     * \return The reflected segment.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2 & cv) const
    {
      Construct_min_vertex_2 max_vertex = construct_min_vertex_2_object();
      Construct_max_vertex_2 max_vertex = construct_max_vertex_2_object();

        
      X_monotone_curve_2 reflected_cv(operator()(max_vertex(cv)),
                                      operator()(min_vertex(cv)));
      return reflected_cv;
    }
  };

  /*! Obtain a Reflect_xy_2 functor object. */
  Reflect_xy_2 reflect_xy_2_object() const
  {
    return Reflect_xy_2();
  }
#endif
    
  /*! A functor for comparing two segments to the right of a point */
  class Compare_y_at_x_right_2 {
  public:
    /*! Compare the y value of two segments immediately to the right of their
     * intersection point.
     * \param cv1 The first segment.
     * \param cv2 The second segment.
     * \param p The intersection point.
     * \pre The point p lies on both segments, and both of them must be also be
     *      defined (lexicographically) to its right.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the right of p: SMALLER, LARGER or EQUAL.
     */

    Comparison_result operator()(const X_monotone_curve_2 & cv1,
                                 const X_monotone_curve_2 & cv2,
                                 const Point_2 & p) const
    {
      Kernel kernel;

      // The two segments must be defined at q and also to its right.
      CGAL_precondition_code(
        Compare_y_at_x_2 compare_y_at_x = kernel.compare_y_at_x_2_object();
        );

      CGAL_precondition(compare_y_at_x(p, cv1) == EQUAL &&
                        compare_y_at_x(p, cv2) == EQUAL);

      CGAL_precondition_code(
        Compare_xy_2 compare_xy = kernel.compare_xy_2_object();
        typename Kernel::Construct_vertex_2 construct_vertex =
          kernel.construct_vertex_2_object();
        const Point_2 & source1 = construct_vertex(cv1, 0);
        const Point_2 & target1 = construct_vertex(cv1, 1);
        const Point_2 & right1 =
          (kernel.less_xy_2_object()(source1, target1)) ? target1 : source1;
        const Point_2 & source2 = construct_vertex(cv2, 0);
        const Point_2 & target2 = construct_vertex(cv2, 1);
        const Point_2 & right2 =
          (kernel.less_xy_2_object()(source2, target2)) ? target2 : source2;
        );

      CGAL_precondition(compare_xy(right1, p) == LARGER &&
                        compare_xy(right2, p) == LARGER);

      // Compare the slopes of the two segments to determine thir relative
      // position immediately to the left of q.
       return kernel.compare_slope_2_object()(cv1, cv2);
    }
  };

  /*! Obtain a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  {
    return Compare_y_at_x_right_2();
  }
  
  ///@}
};

CGAL_END_NAMESPACE

#endif
