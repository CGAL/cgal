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

#ifndef CGAL_ARR_COUNTING_TRAITS_H
#define CGAL_ARR_COUNTING_TRAITS_H

/*! \file
 * A counting traits-class for the arrangement package.
 * This is a meta-traits class. It is parameterized with another traits class
 * and inherits from it. For each traits method it maintains a counter that
 * counts the number of invokations into the method.
 */

#include <iostream>
#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

/*! \class 
 * A model of the ArrangementTraits_2 concept that counts the methods invoked.
 */
template <class Base_traits>
class Arr_counting_traits_2 : public Base_traits {
public:
  enum Operation_id {
    COMPARE_X = 0,
    COMPARE_XY,
    CONSTRUCT_MIN_VERTEX, 
    CONSTRUCT_MAX_VERTEX,
    BOUNDARY_IN_X,
    BOUNDARY_IN_Y,
    IS_VERTICAL, 
    COMPARE_Y_AT_X, 
    EQUAL_POINTS, 
    EQUAL_CURVES, 
    COMPARE_Y_AT_X_LEFT, 
    COMPARE_Y_AT_X_RIGHT, 
    MAKE_X_MONOTONE, 
    SPLIT, 
    INTERSECT, 
    ARE_MERGEABLE, 
    MERGE, 
    CONSTRUCT_OPPOSITE, 
    COMPARE_ENDPOINTS_XY,
    NUMBER_OF_OPERATIONS
  };

  typedef Base_traits                           Base;
  typedef Arr_counting_traits_2<Base>           Self;

  /*! Default constructor */
  Arr_counting_traits_2() : Base()
  { increment(); bzero(m_counters, sizeof(m_counters)); }

  /*! Obtain the counter of the given operation */
  unsigned int get_count(Operation_id id) const { return m_counters[id]; }
  
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
    unsigned int & m_counter;
  public:
    Compare_x_2(const Base * base, unsigned int & counter) :
      m_object(base->compare_x_2_object()), m_counter(counter) {}

    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    { ++m_counter; return m_object(p1, p2); }

    Comparison_result operator()(const Point_2 & p,
                                 const X_monotone_curve_2 & xc, Curve_end ind)
      const
    { ++m_counter; return m_object(p, xc, ind); }

    Comparison_result operator()(const X_monotone_curve_2 & xc1, Curve_end ind1,
                                 const X_monotone_curve_2 & xc2, Curve_end ind2)
      const
    { ++m_counter; return m_object(xc1, ind1, xc2, ind2); }
  };

  /*! Compare two points lexigoraphically; by x, then by y */
  class Compare_xy_2 {
  private:
    typename Base::Compare_xy_2 m_object;
    mutable unsigned int & m_counter;
  public:
    Compare_xy_2(const Base * base, unsigned int & counter) :
      m_object(base->compare_xy_2_object()), m_counter(counter) {}
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    { ++m_counter; return m_object(p1, p2); }
  };

  /*! Determine whether an endpoint of an x-monotone curve lies on an
   * x-boundary.
   */
  class Boundary_in_x_2 {
  private:
    typename Base::Boundary_in_x_2 m_object;
    mutable unsigned int & m_counter;
  public:
    Boundary_in_x_2(const Base * base, unsigned int & counter) :
      m_object(base->boundary_in_x_2_object()), m_counter(counter) {}
    Boundary_type operator()(const X_monotone_curve_2 & xc, Curve_end ind) const
    { ++m_counter; return m_object(xc, ind); }
  };

  /*! Determine whether an endpoint of an x-monotone curve lies on an
   * y-boundary.
   */
  class Boundary_in_y_2 {
  private:
    typename Base::Boundary_in_y_2 m_object;
    mutable unsigned int & m_counter;
  public:
    Boundary_in_y_2(const Base * base, unsigned int & counter) :
      m_object(base->boundary_in_y_2_object()), m_counter(counter) {}
    Boundary_type operator()(const X_monotone_curve_2 & xc, Curve_end ind) const
    { ++m_counter; return m_object(xc, ind); }
  };
  
  /*! Obtain the left endpoint of a given x-monotone curve */
  class Construct_min_vertex_2 {
  private:
    typename Base::Construct_min_vertex_2 m_object;
    unsigned int & m_counter;
  public:
    Construct_min_vertex_2(const Base * base, unsigned int & counter) :
      m_object(base->construct_min_vertex_2_object()), m_counter(counter) {}
    const Point_2 operator()(const X_monotone_curve_2 & xc) const
    { ++m_counter; return m_object(xc); }
  };

  /*! Obtain the right endpoint of a given x-monotone curve */
  class Construct_max_vertex_2 {
  private:
    typename Base::Construct_max_vertex_2 m_object;
    unsigned int & m_counter;
  public:
    Construct_max_vertex_2(const Base * base, unsigned int & counter) :
      m_object(base->construct_max_vertex_2_object()), m_counter(counter) {}
    const Point_2 operator()(const X_monotone_curve_2 & xc) const
    { ++m_counter; return m_object(xc); }
  };
  
  /*! Check whether a given x-monotone curve is vertical */
  class Is_vertical_2 {
  private:
    typename Base::Is_vertical_2 m_object;
    unsigned int & m_counter;
  public:
    Is_vertical_2(const Base * base, unsigned int & counter) :
      m_object(base->is_vertical_2_object()), m_counter(counter) {}
    bool operator()(const X_monotone_curve_2 & xc) const
    { ++m_counter; return m_object(xc); }
  };
  
  /*! Return the location of a given point with respect to a given x-monotone
   * curve
   */
  class Compare_y_at_x_2 {
  private:
    typename Base::Compare_y_at_x_2 m_object;
    unsigned int & m_counter;
  public:
    Compare_y_at_x_2(const Base * base, unsigned int & counter) :
      m_object(base->compare_y_at_x_2_object()), m_counter(counter) {}
    Comparison_result operator()(const Point_2 & p,
                                 const X_monotone_curve_2 & xc) const
    { ++m_counter; return m_object(p, xc); }

    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 const X_monotone_curve_2 & xc2, 
                                 Curve_end ind) const
    { ++m_counter; return m_object(xc1, xc2, ind); }

    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 Curve_end ind1,
                                 const X_monotone_curve_2 & xc2, 
                                 Curve_end ind2) const
    { ++m_counter; return m_object(xc1, ind1, xc2, ind2); }
  };
  
  /*! Check if two x-monotone curves or if two points are identical */
  class Equal_2 {
  private:
    typename Base::Equal_2 m_object;
    unsigned int & m_counter1;
    unsigned int & m_counter2;
  public:
    Equal_2(const Base * base, unsigned int& counter1, unsigned int& counter2) :
      m_object(base->equal_2_object()),
      m_counter1(counter1), m_counter2(counter2)
    {}
    bool operator()(const X_monotone_curve_2 & xc1,
                    const X_monotone_curve_2 & xc2) const
    { ++m_counter1; return m_object(xc1, xc2); }
    bool operator()(const Point_2 & p1, const Point_2 & p2) const
    { ++m_counter2; return m_object(p1, p2); }    
  };

  /*! Compare the y value of two x-monotone curves immediately to the left of
   * their intersection point
   */
  class Compare_y_at_x_left_2 {
  private:
    typename Base::Compare_y_at_x_left_2 m_object;
    unsigned int & m_counter;
  public:
    Compare_y_at_x_left_2(const Base * base, unsigned int & counter) :
      m_object(base->compare_y_at_x_left_2_object()), m_counter(counter) {}
    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 const X_monotone_curve_2 & xc2,
                                 const Point_2 & p) const
    { ++m_counter; return m_object(xc1, xc2, p); }
  };

  /*! Compare the y value of two x-monotone curves immediately to the right of
   * their intersection point
   */
  class Compare_y_at_x_right_2 {
  private:
    typename Base::Compare_y_at_x_right_2 m_object;
    unsigned int & m_counter;
  public:
    Compare_y_at_x_right_2(const Base * base, unsigned int & counter) :
      m_object(base->compare_y_at_x_right_2_object()), m_counter(counter) {}
    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 const X_monotone_curve_2 & xc2,
                                 const Point_2 & p) const
    { ++m_counter; return m_object(xc1, xc2, p); }
  };
  
  /*! Split a curve into x-monotone pieces */
  class Make_x_monotone_2 {
  private:
    typename Base::Make_x_monotone_2 m_object;
    unsigned int & m_counter;
  public:
    Make_x_monotone_2(Base * base, unsigned int & counter) :
      m_object(base->make_x_monotone_2_object()), m_counter(counter) {}
    template<class OutputIterator>
    OutputIterator operator()(const Curve_2 & cv, OutputIterator oi)
    { ++m_counter; return m_object(cv, oi); }
  };

  /*! Split an x-monotone curve into two */
  class Split_2 {
  private:
    typename Base::Split_2 m_object;
    unsigned int & m_counter;
  public:
    Split_2(Base * base, unsigned int & counter) :
      m_object(base->split_2_object()), m_counter(counter) {}
    void operator()(const X_monotone_curve_2 & xc, const Point_2 & p,
                    X_monotone_curve_2 & xc1, X_monotone_curve_2 & xc2)
    { ++m_counter; m_object(xc, p, xc1, xc2); }
  };

  /*! compute intersections */
  class Intersect_2 {
  private:
    typename Base::Intersect_2 m_object;
    unsigned int & m_counter;
  public:
    Intersect_2(Base * base, unsigned int & counter) :
      m_object(base->intersect_2_object()), m_counter(counter) {}
    template<class OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2 & xc1,
                              const X_monotone_curve_2 & xc2,
                              OutputIterator oi)
    { ++m_counter; return m_object(xc1, xc2, oi); }
  };

  /*! Test whether two x-monotone curves are mergeable */
  class Are_mergeable_2 {
  private:
    typename Base::Are_mergeable_2 m_object;
    unsigned int & m_counter;
  public:
    Are_mergeable_2(const Base * base, unsigned int & counter) :
      m_object(base->are_mergeable_2_object()), m_counter(counter) {}
    bool operator()(const X_monotone_curve_2 & xc1,
                    const X_monotone_curve_2 & xc2) const
    { ++m_counter; return m_object(xc1, xc2); }
  };

  /*! Merge two x-monotone curves into one */
  class Merge_2 {
  private:
    typename Base::Merge_2 m_object;
    unsigned int & m_counter;
  public:
    Merge_2(Base * base, unsigned int & counter) :
      m_object(base->merge_2_object()), m_counter(counter) {}
    void operator()(const X_monotone_curve_2 & xc1,
                    const X_monotone_curve_2 & xc2,
                    X_monotone_curve_2 & xc)
    { ++m_counter; m_object(xc1, xc2, xc); }
  };

  /*! Construct an opposite x-monotone curve */
  class Construct_opposite_2 {
  private:
    typename Base::Construct_opposite_2 m_object;
    unsigned int & m_counter;
  public:
    Construct_opposite_2(const Base * base, unsigned int & counter) :
      m_object(base->construct_opposite_2_object()), m_counter(counter) {}
    X_monotone_curve_2 operator()(const X_monotone_curve_2 & xc)
    { ++m_counter; return m_object(xc); }
  };

  /*! Compare the two endpoints of a given curve lexigoraphically */
  class Compare_endpoints_xy_2 {
  private:
    typename Base::Compare_endpoints_xy_2 m_object;
    unsigned int & m_counter;
  public:
    Compare_endpoints_xy_2(const Base * base, unsigned int & counter) :
      m_object(base->compare_endpoints_xy_2_object()), m_counter(counter) {}
    Comparison_result operator()(const X_monotone_curve_2 & xc)
    { ++m_counter; return m_object(xc); }
  };

  //@}

  /// \name Obtain the appropriate functor
  //@{
  Compare_x_2 compare_x_2_object() const
  { return Compare_x_2(this, m_counters[COMPARE_X]); }
  
  Compare_xy_2 compare_xy_2_object() const
  { return Compare_xy_2(this, m_counters[COMPARE_XY]); }

  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(this, m_counters[CONSTRUCT_MIN_VERTEX]); }

  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(this, m_counters[CONSTRUCT_MAX_VERTEX]); }
  
  Boundary_in_x_2 boundary_in_x_2_object() const
  { return Boundary_in_x_2(this, m_counters[BOUNDARY_IN_X]); }

  Boundary_in_y_2 boundary_in_y_2_object() const
  { return Boundary_in_y_2(this, m_counters[BOUNDARY_IN_Y]); }  

  Is_vertical_2 is_vertical_2_object() const
  { return Is_vertical_2(this, m_counters[IS_VERTICAL]); }
  
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(this, m_counters[COMPARE_Y_AT_X]); }
  
  Equal_2 equal_2_object() const
  { return Equal_2(this, m_counters[EQUAL_POINTS], m_counters[EQUAL_CURVES]); }

  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(this, m_counters[COMPARE_Y_AT_X_LEFT]); }

  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(this, m_counters[COMPARE_Y_AT_X_RIGHT]); }
  
  Make_x_monotone_2 make_x_monotone_2_object()
  { return Make_x_monotone_2(this, m_counters[MAKE_X_MONOTONE]); }

  Split_2 split_2_object()
  { return Split_2(this, m_counters[SPLIT]); }

  Intersect_2 intersect_2_object()
  { return Intersect_2(this, m_counters[INTERSECT]); }

  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(this, m_counters[ARE_MERGEABLE]); }

  Merge_2 merge_2_object()
  { return Merge_2(this, m_counters[MERGE]); }

  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(this, m_counters[CONSTRUCT_OPPOSITE]); }

  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(this, m_counters[COMPARE_ENDPOINTS_XY]); }

  //@}

  static unsigned int increment(bool doit = true)
  {
    static unsigned int counter = 0;
    if (doit) ++counter;
    return counter;
  }
  
private:
  mutable unsigned int m_counters[NUMBER_OF_OPERATIONS];
};

template <class Out_stream, class Base_traits>
inline
Out_stream & operator<<(Out_stream & os,
                        const Arr_counting_traits_2<Base_traits> & traits)
{
  typedef Arr_counting_traits_2<Base_traits>            Traits;
  unsigned int sum = 0;
  unsigned int i;
  for (i = 0; i < Traits::NUMBER_OF_OPERATIONS; ++i)
    sum += traits.get_count(static_cast<typename Traits::Operation_id>(i));
  os << "count[COMPARE_X] = "
     << traits.get_count(Traits::COMPARE_X) << std::endl
     << "count[COMPARE_XY] = "
     << traits.get_count(Traits::COMPARE_XY) << std::endl
     << "count[CONSTRUCT_MIN_VERTEX] = "
     << traits.get_count(Traits::CONSTRUCT_MIN_VERTEX) << std::endl
     << "count[CONSTRUCT_MAX_VERTEX] = "
     << traits.get_count(Traits::CONSTRUCT_MAX_VERTEX) << std::endl
     << "count[IS_VERTICAL] = "
     << traits.get_count(Traits::IS_VERTICAL) << std::endl
     << "count[COMPARE_Y_AT_X] = "
     << traits.get_count(Traits::COMPARE_Y_AT_X) << std::endl
     << "count[EQUAL_POINTS] = "
     << traits.get_count(Traits::EQUAL_POINTS) << std::endl
     << "count[EQUAL_CURVES] = "
     << traits.get_count(Traits::EQUAL_CURVES) << std::endl
     << "count[COMPARE_Y_AT_X_LEFT] = "
     << traits.get_count(Traits::COMPARE_Y_AT_X_LEFT) << std::endl
     << "count[COMPARE_Y_AT_X_RIGHT] = "
     << traits.get_count(Traits::COMPARE_Y_AT_X_RIGHT) << std::endl
     << "count[MAKE_X_MONOTONE] = "
     << traits.get_count(Traits::MAKE_X_MONOTONE) << std::endl
     << "count[SPLIT] = "
     << traits.get_count(Traits::SPLIT) << std::endl
     << "count[INTERSECT] = "
     << traits.get_count(Traits::INTERSECT) << std::endl
     << "count[ARE_MERGEABLE] = "
     << traits.get_count(Traits::ARE_MERGEABLE) << std::endl
     << "count[MERGE] = "
     << traits.get_count(Traits::MERGE) << std::endl
     << "count[CONSTRUCT_OPPOSITE] = "
     << traits.get_count(Traits::CONSTRUCT_OPPOSITE) << std::endl
     << "count[COMPARE_ENDPOINTS_XY] = "
     << traits.get_count(Traits::COMPARE_ENDPOINTS_XY) << std::endl
     << "total = " << sum << std::endl
     << "No. of traits constructed = " << Traits::increment(false) << std::endl;
  return os;
}

CGAL_END_NAMESPACE

#endif
