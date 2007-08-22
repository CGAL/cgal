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

private:
  typedef Base_traits                           Base;
  typedef Arr_tracing_traits_2<Base>            Self;

  bool m_trace[NUMBER_OF_OPERATIONS];
  
public:
  /*! Default constructor */
  Arr_tracing_traits_2() :
    Base()
  {
    unsigned int i;
    for (i = 0; i < NUMBER_OF_OPERATIONS; ++i) m_trace[i] = true;
  }

  /*! Enable the trace of a traits operation
   * \param id the operation identifier
   */
  void enable_trace(Operation_id id) { m_trace[id] = true; }

  /*! Disable the trace of a traits operation
   * \param id the operation identifier
   */
  void disable_trace(Operation_id id) { m_trace[id] = false; }
  
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
  public:
    Compare_x_2(const Base * base) : m_object(base->compare_x_2_object()) {}

    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      std::cout << "compare_x 1" << std::endl
                << "  p1: " << p1 << std::endl
                << "  p2: " << p2 << std::endl;
      Comparison_result cr = m_object(p1, p2);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }

    Comparison_result operator()(const Point_2 & p,
                                 const X_monotone_curve_2 & xc, Curve_end ind)
      const
    {
      std::cout << "compare_x 2" << std::endl
                << "  p: " << p << std::endl
                << "  ind: " << ind << ", xc: " << xc << std::endl;
      Comparison_result cr = m_object(p, xc, ind);
      std::cout << "  result: " << std::endl;
      return cr;
    }

    Comparison_result operator()(const X_monotone_curve_2 & xc1, Curve_end ind1,
                                 const X_monotone_curve_2 & xc2, Curve_end ind2)
      const
    {
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
  public:
    Compare_xy_2(const Base * base) : m_object(base->compare_xy_2_object()) {}
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
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
  public:
    Boundary_type operator()(const X_monotone_curve_2 & xc, Curve_end ind) const
    {
      std::cout << "boundary_in_x" << std::endl
                << "  ind: " << ind << ", xc: " << xc << std::endl;
      Boundary_type bt = m_object(xc, ind);
      std::cout << "  result: " << bt << std::endl;
      return bt;
    }
  };

  /*! Determine whether an endpoint of an x-monotone curve lies on an
   * y-boundary.
   */
  class Boundary_in_y_2 {
  private:
    typename Base::Boundary_in_y_2 m_object;
  public:
    Boundary_type operator()(const X_monotone_curve_2 & xc, Curve_end ind) const
    {
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
  public:
    Construct_min_vertex_2(const Base * base) :
      m_object(base->construct_min_vertex_2_object()) {}
    const Point_2 operator()(const X_monotone_curve_2 & xc) const
    {
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
  public:
    Construct_max_vertex_2(const Base * base) :
      m_object(base->construct_max_vertex_2_object()) {}
    const Point_2 operator()(const X_monotone_curve_2 & xc) const
    {
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
  public:
    Is_vertical_2(const Base * base) : m_object(base->is_vertical_2_object()) {}
    bool operator()(const X_monotone_curve_2 & xc) const
    {
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
  public:
    Compare_y_at_x_2(const Base * base) :
      m_object(base->compare_y_at_x_2_object()) {}
    Comparison_result operator()(const Point_2 & p,
                                 const X_monotone_curve_2 & xc) const
    {
      std::cout << "compare_y_at_x 1" << std::endl
                << "  p: " << p << std::endl
                << "  xc: " << xc << std::endl;
      Comparison_result cr = m_object(p, xc);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }

    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 const X_monotone_curve_2 & xc2, 
                                 Curve_end ind) const
    {
      std::cout << "compare_y_at_x 2" << std::endl
                << "  ind: " << ind << std::endl
                << "  xc1: " << xc1 << std::endl
                << "  xc2: " << xc2 << std::endl;
      Comparison_result cr = m_object(xc1, xc2, ind);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }

    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 Curve_end ind1,
                                 const X_monotone_curve_2 & xc2, 
                                 Curve_end ind2) const
    {
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
  public:
    Equal_2(const Base * base) : m_object(base->equal_2_object()) {}
    bool operator()(const X_monotone_curve_2 & xc1,
                    const X_monotone_curve_2 & xc2) const
    {
      std::cout << "equal 1" << std::endl
                << "  xc1: " << xc1 << std::endl
                << "  xc1: " << xc1 << std::endl;
      bool equal = m_object(xc1, xc2);
      std::cout << "  result: " << equal << std::endl;
      return equal;
    }
    bool operator()(const Point_2 & p1, const Point_2 & p2) const
    {
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
  public:
    Compare_y_at_x_left_2(const Base * base) :
      m_object(base->compare_y_at_x_left_2_object()) {}
    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 const X_monotone_curve_2 & xc2,
                                 const Point_2 & p) const
    {
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
  public:
    Compare_y_at_x_right_2(const Base * base) :
      m_object(base->compare_y_at_x_right_2_object()) {}
    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 const X_monotone_curve_2 & xc2,
                                 const Point_2 & p) const
    {
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
  public:
    Make_x_monotone_2(Base * base) :
      m_object(base->make_x_monotone_2_object()) {}
    template<typename OutputIterator>
    OutputIterator operator()(const Curve_2 & cv, OutputIterator oi)
    {
      std::cout << "make_x_monotone" << std::endl
                << "  cv: " << cv << std::endl;        
      OutputIterator res_oi = m_object(cv, oi);
#if 1
      return res_oi;
#else
      /*! \todo the following code will fail to compile if OutputIterator is
       * a back_inserter (e.g., of a list container).
       * What is the right way to write generic code in this case?
       */
      if (res_oi == oi) return res_oi;

      unsigned int i = 0;
      OutputIterator it;
      for (it = oi; it != res_oi; ++it) {
        X_monotone_curve_2 xc;
        if (assign (xc, *it)) {
          std::cout << "  result[" << i++ << "]: xc: " << xc << std::endl;
          continue;
        }
      }
#endif
      return res_oi;
    }
  };

  /*! Split an x-monotone curve into two */
  class Split_2 {
  private:
    typename Base::Split_2 m_object;
  public:
    Split_2(Base * base) : m_object(base->split_2_object()) {}
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
  public:
    Intersect_2(Base * base) : m_object(base->intersect_2_object()) {}
    template<typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2 & xc1,
                              const X_monotone_curve_2 & xc2,
                              OutputIterator oi)
    {
      std::cout << "intersect" << std::endl
                << "  xc1: " << xc1 << std::endl
                << "  xc2: " << xc2 << std::endl;
      std::list<CGAL::Object> container;
      m_object(xc1, xc2, std::back_inserter(container));
      if (container.empty()) return oi;

      unsigned int i = 0;
      std::list<CGAL::Object>::iterator it;
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
  public:
    Are_mergeable_2(const Base * base) :
      m_object(base->are_mergeable_2_object()) {}
    bool operator()(const X_monotone_curve_2 & xc1,
                    const X_monotone_curve_2 & xc2) const
    {
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
  public:
    Merge_2(Base * base) : m_object(base->merge_2_object()) {}
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
  public:
    Construct_opposite_2(const Base * base) :
      m_object(base->construct_opposite_2_object()) {}
    X_monotone_curve_2 operator()(const X_monotone_curve_2 & xc)
    {
      return m_object(xc);
    }
  };

  /*! Compare the two endpoints of a given curve lexigoraphically */
  class Compare_endpoints_xy_2 {
  private:
    typename Base::Compare_endpoints_xy_2 m_object;
  public:
    Compare_endpoints_xy_2(const Base * base) :
      m_object(base->compare_endpoints_xy_2_object()) {}
    Comparison_result operator()(const X_monotone_curve_2 & xc)
    {
      std::cout << "compare_endpoints_xy" << std::endl
                << "  xc: " << xc << std::endl;
      Comparison_result cr = m_object(xc);
      std::cout << "  result: " << cr;
      return cr;
    }
  };

  //@}

  /// \name Obtain the appropriate functor
  //@{
  Compare_x_2 compare_x_2_object() const
  { return Compare_x_2(this); }
  
  Compare_xy_2 compare_xy_2_object() const
  { return Compare_xy_2(this); }

  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(this); }

  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(this); }
  
  Is_vertical_2 is_vertical_2_object() const
  { return Is_vertical_2(this); }
  
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(this); }
  
  Equal_2 equal_2_object() const
  { return Equal_2(this); }

  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(this); }

  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(this); }
  
  Make_x_monotone_2 make_x_monotone_2_object()
  { return Make_x_monotone_2(this); }

  Split_2 split_2_object()
  { return Split_2(this); }

  Intersect_2 intersect_2_object()
  { return Intersect_2(this); }

  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(this); }

  Merge_2 merge_2_object()
  { return Merge_2(this); }

  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(this); }

  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(this); }

  //@}
};

template <class OutputStream>
OutputStream & operator<<(OutputStream & os, Comparison_result cr)
{
  os << ((cr == SMALLER) ? "SMALLER" : (cr == EQUAL) ? "EQUAL" : "LARGER");
  return os;
}

/*! Inserter for the spherical_arc class used by the traits-class */
#if 0
template <class T_Kernel, class OutputStream>
OutputStream & operator<<(OutputStream & os,
                          const Arr_extended_direction_3<T_Kernel> & p)
{
  CGAL::To_double<typename T_Kernel::FT> todouble;
  os << "("
     << static_cast<float>(todouble(p.dx())) << ","
     << static_cast<float>(todouble(p.dy())) << ","
     << static_cast<float>(todouble(p.dz()))
     << ")"
     << ", "
     <<
    (p.is_min_boundary() ? "min" :
     p.is_max_boundary() ? "max" :
     p.is_mid_boundary() ? "dis" : "reg");
  return os;
}
#else
template <class T_Traits, class OutputStream>
OutputStream & operator<<(OutputStream & os,
                          const typename
                            Arr_tracing_traits_2<T_Traits>::Point_2 & p)
{
  typename T_Traits::Point_2 p_base = p;
  os << "(" << p_base << ")"
     << ", "
     <<
    (p.is_min_boundary() ? "min" :
     p.is_max_boundary() ? "max" :
     p.is_mid_boundary() ? "dis" : "reg");
  return os;
}
#endif

/*! Inserter for the spherical_arc class used by the traits-class */
#if 0
template <class T_Kernel, class OutputStream>
OutputStream & operator<<(OutputStream & os,
                          const Arr_spherical_arc_3<T_Kernel> & xc)
{
  os << "("
     << xc.left() << "," << xc.right()
     << ")"
     << ", " << (xc.is_vertical() ? " |" : "!|")
     << ", " << (xc.is_directed_right() ? "=>" : "<=");
  return os;
}
#else
template <class T_Traits, class OutputStream>
OutputStream & operator<<(OutputStream & os,
                          const typename
                            Arr_tracing_traits_2<T_Traits>::
                            X_monotone_curve_2 & xc)
{
  typename T_Traits::X_monotone_curve_2 xc_base = xc;
  os << "(" << xc << ")"
     << ", " << (xc.is_vertical() ? " |" : "!|")
     << ", " << (xc.is_directed_right() ? "=>" : "<=");
  return os;
}
#endif

CGAL_END_NAMESPACE

#endif
