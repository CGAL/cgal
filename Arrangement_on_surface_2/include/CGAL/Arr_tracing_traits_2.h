// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL$
// $Id$
//
// Author(s)     : Efi Fogel    <efif@post.tau.ac.il>

#ifndef CGAL_ARR_TRACING_TRAITS_H
#define CGAL_ARR_TRACING_TRAITS_H

/*! \file
 * A tracing traits-class for the arrangement package.
 * This is a meta-traits class. It is parameterized with another traits class
 * and inherits from it. For each traits method it prints out its input
 * parameters and its output result
 */

#include <iostream>
#include <list>

#include <CGAL/basic.h>
#include <CGAL/Arr_enums.h>

CGAL_BEGIN_NAMESPACE

/*! \class 
 * A model of the ArrangementTraits_2 concept that counts the methods invoked.
 */
template <typename Base_traits>
class Arr_tracing_traits_2 : public Base_traits {
public:
  enum Operation_id {
    COMPARE_X_OP = 0,
    COMPARE_XY_OP,
    CONSTRUCT_MIN_VERTEX_OP,
    CONSTRUCT_MAX_VERTEX_OP,
    BOUNDARY_IN_X_OP,
    BOUNDARY_IN_Y_OP,
    IS_VERTICAL_OP,
    COMPARE_Y_AT_X_OP,
    EQUAL_POINT_OP,
    EQUAL_CURVE_OP,
    COMPARE_Y_AT_X_LEFT_OP,
    COMPARE_Y_AT_X_RIGHT_OP,
    MAKE_X_MONOTONE_OP,
    SPLIT_OP,
    INTERSECT_OP,
    ARE_MERGEABLE_OP,
    MERGE_OP,
    CONSTRUCT_OPPOSITE_OP,
    COMPARE_ENDPOINTS_XY_OP,
    NUMBER_OF_OPERATIONS
  };

private:
  typedef Base_traits                           Base;
  typedef Arr_tracing_traits_2<Base>            Self;

  /*! A set of bits that indicate whether operations should be traced */
  unsigned int m_flags;
  
  bool compare_x_op() const
  { return m_flags & (0x1 << COMPARE_X_OP); }

  bool compare_xy_op() const
  { return m_flags & (0x1 << COMPARE_XY_OP); }
  
  bool construct_min_vertex_op() const
  { return m_flags & (0x1 << CONSTRUCT_MIN_VERTEX_OP); }
  
  bool construct_max_vertex_op() const
  { return m_flags & (0x1 << CONSTRUCT_MAX_VERTEX_OP); }
  
  bool boundary_in_x_op() const
  { return m_flags & (0x1 << BOUNDARY_IN_X_OP); }
  
  bool boundary_in_y_op() const
  { return m_flags & (0x1 << BOUNDARY_IN_Y_OP); }
  
  bool is_vertical_op() const
  { return m_flags & (0x1 << IS_VERTICAL_OP); }
  
  bool compare_y_at_x_op() const
  { return m_flags & (0x1 << COMPARE_Y_AT_X_OP); }
  
  bool equal_point_op() const
  { return m_flags & (0x1 << EQUAL_POINT_OP); }
  
  bool equal_curve_op() const
  { return m_flags & (0x1 << EQUAL_CURVE_OP); }
  
  bool compare_y_at_x_left_op() const
  { return m_flags & (0x1 << COMPARE_Y_AT_X_LEFT_OP); }
  
  bool compare_y_at_x_right_op() const
  { return m_flags & (0x1 << COMPARE_Y_AT_X_RIGHT_OP); }
  
  bool make_x_monotone_op() const
  { return m_flags & (0x1 << MAKE_X_MONOTONE_OP); }
  
  bool split_op() const
  { return m_flags & (0x1 << SPLIT_OP); }
  
  bool intersect_op() const
  { return m_flags & (0x1 << INTERSECT_OP); }
  
  bool are_mergeable_op() const
  { return m_flags & (0x1 << ARE_MERGEABLE_OP); }
  
  bool merge_op() const
  { return m_flags & (0x1 << MERGE_OP); }
  
  bool construct_opposite_op() const
  { return m_flags & (0x1 << CONSTRUCT_OPPOSITE_OP); }
  
  bool compare_endpoints_xy_op() const
  { return m_flags & (0x1 << COMPARE_ENDPOINTS_XY_OP); }
  
public:
  /*! Default constructor */
  Arr_tracing_traits_2() :
    Base()
  {
    enable_all_traces();
  }

  /*! Enable the trace of a traits operation
   * \param id the operation identifier
   */
  void enable_trace(Operation_id id) { m_flags[id] |= 0x1 << id; }

  /*! Enable the trace of all traits operations
   * \param id the operation identifier
   */
  void enable_all_traces() { m_flags = 0xffffffff; }
  
  /*! Disable the trace of a traits operation
   * \param id the operation identifier
   */
  void disable_trace(Operation_id id) { m_flags[id] &= ~(0x1 << id); }

  /*! Disable the trace of all traits operations
   * \param id the operation identifier
   */
  void disable_all_traces() { m_flags = 0x0; }
  
  /// \name Types and functors inherited from the base
  //@{

  // Traits types:
  typedef typename Base::Has_left_category      Has_left_category;
  typedef typename Base::Has_boundary_category  Has_boundary_category;
  typedef typename Base::Has_merge_category     Has_merge_category;
  
  typedef typename Base::Point_2                Point_2;
  typedef typename Base::X_monotone_curve_2     X_monotone_curve_2;
  typedef typename Base::Curve_2                Curve_2;

  /*! Compare the x-coordinates of two points */
  class Compare_x_2 {
  private:
    typename Base::Compare_x_2 m_object;
    bool m_enabled;
    
  public:
    /*! Construct */
    Compare_x_2(const Base * base, bool enabled = true) :
      m_object(base->compare_x_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param p1 first point
     * \param p2 second point
     * \return the comparison result
     */
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      if (!m_enabled) return m_object(p1, p2);
      std::cout << "compare_x 1" << std::endl
                << "  p1: " << p1 << std::endl
                << "  p2: " << p2 << std::endl;
      Comparison_result cr = m_object(p1, p2);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }

    /*! Operate
     * \param p the first point 
     * \param xc the curve the end of which is to be compared
     * \param ind the curve-end index
     * \return the comparison result
     */
    Comparison_result operator()(const Point_2 & p,
                                 const X_monotone_curve_2 & xc, Curve_end ind)
      const
    {
      if (!m_enabled) return m_object(p, xc, ind);
      std::cout << "compare_x 2" << std::endl
                << "  p: " << p << std::endl
                << "  ind: " << ind << ", xc: " << xc << std::endl;
      Comparison_result cr = m_object(p, xc, ind);
      std::cout << "  result: " << std::endl;
      return cr;
    }

    /*! Operate
     * \param xc1 the first curve the end of which is to be compared
     * \param ind1 the index of the end of the first curve
     * \param xc2 the second curve the end of which is to be compared
     * \param ind2 the index of the end of the second curve
     * \return the comparison result
     */
    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 Curve_end ind1,
                                 const X_monotone_curve_2 & xc2,
                                 Curve_end ind2) const
    {
      if (!m_enabled) return m_object(xc1, ind1, xc2, ind2);
      std::cout << "compare_x 2" << std::endl
                << "  ind1: " << ind1 << ", xc1: " << xc1 << std::endl
                << "  ind2: " << ind2 << ", xc2: " << xc2 << std::endl;
      Comparison_result cr = m_object(xc1, ind1, xc2, ind2);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /*! Compare two points lexigoraphically; by x, then by y */
  class Compare_xy_2 {
  private:
    typename Base::Compare_xy_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_xy_2(const Base * base, bool enabled = true) :
      m_object(base->compare_xy_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param p1 the first point
     * \param p2 the second point
     * \return the comparison result
     */
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      if (!m_enabled) return m_object(p1, p2);
      std::cout << "compare_xy" << std::endl
                << "  p1: " << p1 << std::endl
                << "  p2: " << p2 << std::endl;
      Comparison_result cr = m_object(p1, p2);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /*! Determine whether an endpoint of an x-monotone curve lies on an
   * x-boundary.
   */
  class Boundary_in_x_2 {
  private:
    typename Base::Boundary_in_x_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Boundary_in_x_2(const Base * base, bool enabled = true) :
      m_object(base->boundary_in_x_2_object()), m_enabled(enabled)
    {}
    
    /*! Operate
     * \param xc the curve the end of which is tested
     * \param ind the curve-end index
     * \return the boundary type
     */
    Boundary_type operator()(const X_monotone_curve_2 & xc,
                             Curve_end ind) const
    {
      if (!m_enabled) return m_object(xc, ind);
      std::cout << "boundary_in_x" << std::endl
                << "  ind: " << ind << ", xc: " << xc << std::endl;
      Boundary_type bt = m_object(xc, ind);
      std::cout << "  result: " << bt << std::endl;
      return bt;
    }
  };

  /*! Determine whether an endpoint of an x-monotone curve lies on a
   * y-boundary.
   */
  class Boundary_in_y_2 {
  private:
    typename Base::Boundary_in_y_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Boundary_in_y_2(const Base * base, bool enabled = true) :
      m_object(base->boundary_in_y_2_object()), m_enabled(enabled) {}
    
    /*! Operate
     * \param xc the curve the end of which is tested
     * \param ind the curve-end index
     * \return the boundary type
     */
    Boundary_type operator()(const X_monotone_curve_2 & xc, Curve_end ind) const
    {
      if (!m_enabled) return m_object(xc, ind);
        std::cout << "boundary_in_y" << std::endl
                  << "  ind: " << ind << ", xc: " << xc << std::endl;
      Boundary_type bt = m_object(xc, ind);
      std::cout << "  result: " << bt << std::endl;
      return bt;
    }
  };
  
  /*! Obtain the left endpoint of a given x-monotone curve */
  class Construct_min_vertex_2 {
  private:
    typename Base::Construct_min_vertex_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Construct_min_vertex_2(const Base * base, bool enabled = true) :
      m_object(base->construct_min_vertex_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xc the curev the left endpoint of which is obtained
     * \return the left endpoint
     */
    const Point_2 operator()(const X_monotone_curve_2 & xc) const
    {
      if (!m_enabled) return m_object(xc);
      std::cout << "construct_min_vertex" << std::endl
                << "  xc: " << xc << std::endl;
      Point_2 p = m_object(xc);
      std::cout << "  result: " << p << std::endl;
      return p;
    }
  };

  /*! Obtain the right endpoint of a given x-monotone curve */
  class Construct_max_vertex_2 {
  private:
    typename Base::Construct_max_vertex_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Construct_max_vertex_2(const Base * base, bool enabled = true) :
      m_object(base->construct_max_vertex_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xc the curev the right endpoint of which is obtained
     * \return the right endpoint
     */
    const Point_2 operator()(const X_monotone_curve_2 & xc) const
    {
      if (!m_enabled) return m_object(xc);
      std::cout << "construct_max_vertex" << std::endl
                << "  xc: " << xc << std::endl;
      Point_2 p = m_object(xc);
      std::cout << "  result: " << p << std::endl;      
      return p;
    }
  };
  
  /*! Check whether a given x-monotone curve is vertical */
  class Is_vertical_2 {
  private:
    typename Base::Is_vertical_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Is_vertical_2(const Base * base, bool enabled = true) :
      m_object(base->is_vertical_2_object()), m_enabled(enabled) {}
    
    /*! Operate
     * \param xc the curve
     * \return a Boolean that indicates whether the curve is vertical or not
     */
    bool operator()(const X_monotone_curve_2 & xc) const
    {
      if (!m_enabled) return m_object(xc);
      std::cout << "is_vertical" << std::endl
                << "  xc: " << xc << std::endl;
      bool is_vertical = m_object(xc);
      std::cout << "  result: " << is_vertical << std::endl;
      return is_vertical;
    }
  };
  
  /*! Return the location of a given point with respect to a given x-monotone
   * curve
   */
  class Compare_y_at_x_2 {
  private:
    typename Base::Compare_y_at_x_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_y_at_x_2(const Base * base, bool enabled = true) :
      m_object(base->compare_y_at_x_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param p the point
     * \param xc the curve
     * \return the comparison result
     */
    Comparison_result operator()(const Point_2 & p,
                                 const X_monotone_curve_2 & xc) const
    {
      if (!m_enabled) return m_object(p, xc);
      std::cout << "compare_y_at_x 1" << std::endl
                << "  p: " << p << std::endl
                << "  xc: " << xc << std::endl;
      Comparison_result cr = m_object(p, xc);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }

    /*! Operate
     * \param xc1 the first curve the end point of which is tested
     * \param xc2 the second curve the end point of which is tested
     * \param ind the curve-end index
     * \return the comparison result
     */
    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 const X_monotone_curve_2 & xc2, 
                                 Curve_end ind) const
    {
      if (!m_enabled) return m_object(xc1, xc2, ind);
      std::cout << "compare_y_at_x 2" << std::endl
                << "  ind: " << ind << std::endl
                << "  xc1: " << xc1 << std::endl
                << "  xc2: " << xc2 << std::endl;
      Comparison_result cr = m_object(xc1, xc2, ind);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }

    /*! Operate
     * \param xc1 the first curve
     * \param ind1 the index of the end of the first curve
     * \param xc2 the second curve
     * \param ind2 the index of the end of the second curve
     * \return the comparison result
     */
    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 Curve_end ind1,
                                 const X_monotone_curve_2 & xc2, 
                                 Curve_end ind2) const
    {
      if (!m_enabled) return m_object(xc1, ind1, xc2, ind2);
      std::cout << "compare_y_at_x 3" << std::endl
                << "  ind1: " << ind1 << ", xc1: " << xc1 << std::endl
                << "  ind2: " << ind2 << ", xc2: " << xc2 << std::endl;
      Comparison_result cr = m_object(xc1, ind1, xc2, ind2);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };
  
  /*! Check if two x-monotone curves or if two points are identical */
  class Equal_2 {
  private:
    typename Base::Equal_2 m_object;
    bool m_enabled_point;
    bool m_enabled_curve;

  public:
    /*! Construct */
    Equal_2(const Base * base,
            bool enabled_point = true, bool enabled_curve = true) :
      m_object(base->equal_2_object()),
      m_enabled_point(enabled_point),
      m_enabled_curve(enabled_curve)
    {}

    /*! Operate
     * \param xc1 the first curve
     * \param xc2 the second curve
     * \return true if the x-monotone curves are equal and false otherwise
     */
    bool operator()(const X_monotone_curve_2 & xc1,
                    const X_monotone_curve_2 & xc2) const
    {
      if (!m_enabled_curve) return m_object(xc1, xc2);
      std::cout << "equal 1" << std::endl
                << "  xc1: " << xc1 << std::endl
                << "  xc1: " << xc1 << std::endl;
      bool equal = m_object(xc1, xc2);
      std::cout << "  result: " << equal << std::endl;
      return equal;
    }

    /*! Operate
     * \param p1 the first point
     * \param p2 the second point
     * \return true if the points are equal and false otherwise
     */
    bool operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      if (!m_enabled_point) return m_object(p1, p2);
      std::cout << "equal 2" << std::endl
                << "  p1: " << p1 << std::endl
                << "  p2: " << p2 << std::endl;
      bool equal = m_object(p1, p2);
      std::cout << "  result: " << equal << std::endl;
      return equal;
    }    
  };

  /*! Compare the y value of two x-monotone curves immediately to the left of
   * their intersection point
   */
  class Compare_y_at_x_left_2 {
  private:
    typename Base::Compare_y_at_x_left_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_y_at_x_left_2(const Base * base, bool enabled = true) :
      m_object(base->compare_y_at_x_left_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xc1 the first curve
     * \param xc2 the second curve 
     * \param p the reference point
     * \return the comparison result
     */
    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 const X_monotone_curve_2 & xc2,
                                 const Point_2 & p) const
    {
      if (!m_enabled) return m_object(xc1, xc2, p);
      std::cout << "compare_y_at_x_left" << std::endl
                << "  p: " << p << std::endl
                << "  xc1: " << xc1 << std::endl
                << "  xc2: " << xc2 << std::endl;
      Comparison_result cr = m_object(xc1, xc2, p);
      std::cout << "  result:" << cr << std::endl;
      return cr;
    }
  };

  /*! Compare the y value of two x-monotone curves immediately to the right of
   * their intersection point
   */
  class Compare_y_at_x_right_2 {
  private:
    typename Base::Compare_y_at_x_right_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_y_at_x_right_2(const Base * base, bool enabled = true) :
      m_object(base->compare_y_at_x_right_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xc1 the first curve
     * \param xc2 the second curve
     * \param p the reference point
     * \return the comparison result
     */
    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 const X_monotone_curve_2 & xc2,
                                 const Point_2 & p) const
    {
      if (!m_enabled) return m_object(xc1, xc2, p);
      std::cout << "compare_y_at_x_right" << std::endl
                << "  p: " << p << std::endl
                << "  xc1: " << xc1 << std::endl
                << "  xc2: " << xc2 << std::endl;
      Comparison_result cr = m_object(xc1, xc2, p);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };
  
  /*! Split a curve into x-monotone pieces */
  class Make_x_monotone_2 {
  private:
    typename Base::Make_x_monotone_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Make_x_monotone_2(Base * base, bool enabled = true) :
      m_object(base->make_x_monotone_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param cv the curve
     * \param oi an output iterator that contains the result. It's value
     * type is CGAL::Object, which wraps either an x-monotone curve or a point
     * \return the output iterator
     */
    template<typename OutputIterator>
    OutputIterator operator()(const Curve_2 & cv, OutputIterator oi)
    {
      if (!m_enabled) return m_object(cv, oi);
      std::cout << "make_x_monotone" << std::endl
                << "  cv: " << cv << std::endl;        
      std::list<CGAL::Object> container;
      m_object(cv, std::back_inserter(container));
      if (container.empty()) return oi;

      std::list<CGAL::Object>::iterator it;
      unsigned int i = 0;
      for (it = container.begin(); it != container.end(); ++it) {
        X_monotone_curve_2 xc;
        if (assign (xc, *it)) {
          std::cout << "  result[" << i++ << "]: xc: " << xc << std::endl;
          continue;
        }

        Point_2 p;
        if (assign (p, *it)) {
          std::cout << "  result[" << i++ << "]: p: " << p << std::endl;
          continue;
        }
      }
      
      for (it = container.begin(); it != container.end(); ++it) *oi++ = *it;
      container.clear();
      return oi;
    }
  };

  /*! Split an x-monotone curve into two */
  class Split_2 {
  private:
    typename Base::Split_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Split_2(Base * base, bool enabled = true) :
      m_object(base->split_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xc
     * \param p
     * \param xc1
     * \param xc2
     */
    void operator()(const X_monotone_curve_2 & xc, const Point_2 & p,
                    X_monotone_curve_2 & xc1, X_monotone_curve_2 & xc2)
    {
      std::cout << "split: " << std::endl
                << "  xc: " << xc << std::endl
                << "  p: " << p << std::endl;
      m_object(xc, p, xc1, xc2);
      std::cout << "  result xc1: " << xc1 << std::endl
                << "         xc2: " << xc2 << std::endl;
    }
  };

  /*! compute intersections */
  class Intersect_2 {
  private:
    typename Base::Intersect_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Intersect_2(Base * base, bool enabled = true) :
      m_object(base->intersect_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xc1 the first curve
     * \param xc2 the ssecond curve
     * \param oi an output iterator that contains the result. It's value
     * type is CGAL::Object, which wraps either an x-monotone overlapping
     * curve or pair that consists of an intersection point and its
     * multiplicity
     * \return the output iterator
     */
    template<typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2 & xc1,
                              const X_monotone_curve_2 & xc2,
                              OutputIterator oi)
    {
      if (!m_enabled) return m_object(xc1, xc2, oi);
      std::cout << "intersect" << std::endl
                << "  xc1: " << xc1 << std::endl
                << "  xc2: " << xc2 << std::endl;
      std::list<CGAL::Object> container;
      m_object(xc1, xc2, std::back_inserter(container));
      if (container.empty()) return oi;

      std::list<CGAL::Object>::iterator it;
      unsigned int i = 0;
      for (it = container.begin(); it != container.end(); ++it) {
        X_monotone_curve_2 xc;
        if (assign (xc, *it)) {
          std::cout << "  result[" << i++ << "]: xc: " << xc << std::endl;
          continue;
        }

        std::pair<Point_2,unsigned int> point_pair;
        if (assign (point_pair, *it)) {
          std::cout << "  result[" << i++ << "]: p: " << point_pair.first
                    << ", multiplicity: " << point_pair.second << std::endl;
          continue;
        }
      }
      
      for (it = container.begin(); it != container.end(); ++it) *oi++ = *it;
      container.clear();
      return oi;
    }
  };

  /*! Test whether two x-monotone curves are mergeable */
  class Are_mergeable_2 {
  private:
    typename Base::Are_mergeable_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Are_mergeable_2(const Base * base, bool enabled = true) :
      m_object(base->are_mergeable_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xc1 the first curve
     * \param xc2 the second curve
     * \return true if the the two curve are mergeable and false otherwise.
     * Two curves are mergeable if they have the same underlying theoretical
     * curve
     */
    bool operator()(const X_monotone_curve_2 & xc1,
                    const X_monotone_curve_2 & xc2) const
    {
      if (!m_enabled) return m_object(xc1, xc2);
      std::cout << "are_mergeable" << std::endl
                << "  xc1: " << xc1 << std::endl
                << "  xc2: " << xc2 << std::endl;
      bool are_mergeable = m_object(xc1, xc2);
      std::cout << "  result: " << are_mergeable << std::endl;
      return are_mergeable;
    }
  };

  /*! Merge two x-monotone curves into one */
  class Merge_2 {
  private:
    typename Base::Merge_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Merge_2(Base * base, bool enabled = true) :
      m_object(base->merge_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xc1 the first curve
     * \param xc2 the second curve
     * \param xc the merged curve
     */
    void operator()(const X_monotone_curve_2 & xc1,
                    const X_monotone_curve_2 & xc2,
                    X_monotone_curve_2 & xc)
    {
      std::cout << "merge" << std::endl
                << "  xc1: " << xc1 << std::endl
                << "  xc2: " << xc2 << std::endl;
      return m_object(xc1, xc2, xc);
      std::cout << "  result: " << xc << std::endl;
    }
  };

  /*! Construct an opposite x-monotone curve */
  class Construct_opposite_2 {
  private:
    typename Base::Construct_opposite_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Construct_opposite_2(const Base * base, bool enabled = true) :
      m_object(base->construct_opposite_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xc the curve
     * \return the opposite curve
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2 & xc)
    {
      if (!m_enabled) return m_object(xc);
      std::cout << "construct_opposite" << std::endl;
      << "  xc: " << xc << std::endl;
      X_monotone_curve_2 xc = m_object(xc);
      std::cout << "  result: " << xc << std::endl;
      return xc;
    }
  };

  /*! Compare the two endpoints of a given curve lexigoraphically */
  class Compare_endpoints_xy_2 {
  private:
    typename Base::Compare_endpoints_xy_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_endpoints_xy_2(const Base * base, bool enabled = true) :
      m_object(base->compare_endpoints_xy_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xc the curve
     * \return the comparison result
     */
    Comparison_result operator()(const X_monotone_curve_2 & xc)
    {
      if (!m_enabled) return m_object(xc);
      std::cout << "compare_endpoints_xy" << std::endl
                << "  xc: " << xc << std::endl;
      Comparison_result cr = m_object(xc);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  //@}

  /// \name Obtain the appropriate functor
  //@{
  Compare_x_2 compare_x_2_object() const
  { return Compare_x_2(this, compare_x_op()); }
  
  Compare_xy_2 compare_xy_2_object() const
  { return Compare_xy_2(this, compare_xy_op()); }

  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(this, construct_min_vertex_op()); }

  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(this, construct_max_vertex_op()); }

  Boundary_in_x_2 boundary_in_x_2_object() const
  { return Boundary_in_x_2(this, boundary_in_x_op()); }

  Boundary_in_y_2 boundary_in_y_2_object() const
  { return Boundary_in_y_2(this, boundary_in_y_op()); }
  
  Is_vertical_2 is_vertical_2_object() const
  { return Is_vertical_2(this, is_vertical_op()); }
  
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(this, compare_y_at_x_op()); }
  
  Equal_2 equal_2_object() const
  { return Equal_2(this, equal_point_op(), equal_curve_op()); }

  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(this, compare_y_at_x_left_op()); }

  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(this, compare_y_at_x_right_op()); }
  
  Make_x_monotone_2 make_x_monotone_2_object()
  { return Make_x_monotone_2(this, make_x_monotone_op()); }

  Split_2 split_2_object()
  { return Split_2(this, split_op()); }

  Intersect_2 intersect_2_object()
  { return Intersect_2(this, intersect_op()); }

  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(this, are_mergeable_op()); }

  Merge_2 merge_2_object()
  { return Merge_2(this, merge_op()); }

  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(this, construct_opposite_op()); }

  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(this, compare_endpoints_xy_op()); }

  //@}
};

template <class OutputStream>
OutputStream & operator<<(OutputStream & os, Comparison_result cr)
{
  os << ((cr == SMALLER) ? "SMALLER" : (cr == EQUAL) ? "EQUAL" : "LARGER");
  return os;
}

CGAL_END_NAMESPACE

#endif
