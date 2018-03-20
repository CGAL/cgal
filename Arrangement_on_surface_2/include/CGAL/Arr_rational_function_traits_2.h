// Copyright (c) 2011 Tel-Aviv University (Israel), INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Oren Salzman <orenzalz@post.tau.ac.il >
//                 Michael Hemmer <Michael.Hemmer@sophia.inria.fr>

#ifndef CGAL_ARR_RATIONAL_ARC_TRAITS_D_1_H
#define CGAL_ARR_RATIONAL_ARC_TRAITS_D_1_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/assertions.h>
#include <CGAL/tags.h>
#include <CGAL/Fraction_traits.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Arr_rat_arc/Rational_arc_d_1.h>
#include <CGAL/Arr_rat_arc/Cache.h>



namespace CGAL {

/*! \class
 * A traits class for maintaining an arrangement of bounded arcs (segments) of
 * rational functions of arbitrary degree.
 *
 * The class is templated with two parameters: 
 * Alg_kernel A geometric kernel, where Alg_kernel::FT is the number type
 *            for the coordinates of arrangement vertices, which are algebraic
 *            numbers (defined by Nt_traits::Algebraic).
 * Nt_traits A traits class for performing various operations on the integer,
 *           rational and algebraic types. 
 */
 
template <typename AlgebraicKernel_d_1>
class Arr_rational_function_traits_2
{
public:
  typedef AlgebraicKernel_d_1                           Algebraic_kernel_d_1;
 
  typedef Arr_rational_function_traits_2<Algebraic_kernel_d_1>
                                                        Self;
  typedef Arr_rational_arc::Base_rational_arc_ds_1<Algebraic_kernel_d_1>
                                                        Base_rational_arc_ds_1;

  // Traits objects:
  typedef Arr_rational_arc::Base_rational_arc_d_1<Algebraic_kernel_d_1>
                                                              Base_curve_2;
  typedef Arr_rational_arc::Continuous_rational_arc_d_1<Algebraic_kernel_d_1>
                                                              X_monotone_curve_2;
  typedef Arr_rational_arc::Rational_arc_d_1<Algebraic_kernel_d_1>
                                                              Curve_2;
  typedef Arr_rational_arc::Algebraic_point_2<Algebraic_kernel_d_1>
                                                              Point_2;

  typedef typename Base_rational_arc_ds_1::Algebraic_real_1   Algebraic_real_1;
  typedef typename Base_rational_arc_ds_1::Multiplicity       Multiplicity;
  typedef typename Base_curve_2::Rat_vector                   Rat_vector;

  typedef typename Base_rational_arc_ds_1::Integer            Integer;
  typedef typename Base_rational_arc_ds_1::Rational           Rational; 
  typedef typename Base_rational_arc_ds_1::Polynomial_1       Polynomial_1; 
  typedef typename Base_rational_arc_ds_1::Coefficient        Coefficient; 

  typedef typename Base_rational_arc_ds_1::FT_rat_1           FT_rat_1; 
  typedef typename Base_rational_arc_ds_1::Polynomial_traits_1
    Polynomial_traits_1;
  
  typedef typename Algebraic_kernel_d_1::Bound                Bound; 
  typedef Bound
    Approximate_number_type; 
  
  typedef CGAL::Arr_rational_arc::Rational_function<Algebraic_kernel_d_1>
                                                              Rational_function;
  typedef CGAL::Arr_rational_arc::Cache<Algebraic_kernel_d_1> Cache;

  //Category tags:
  typedef Tag_true Has_left_category;
  typedef Tag_true Has_merge_category;
  typedef Tag_true Has_do_intersect_category;

  typedef Tag_true Has_vertical_segment_category;

  typedef Arr_open_side_tag          Left_side_category;
  typedef Arr_open_side_tag          Bottom_side_category;
  typedef Arr_open_side_tag          Top_side_category;
  typedef Arr_open_side_tag          Right_side_category;

private:
  mutable Cache                   _cache;
  mutable Algebraic_kernel_d_1*   _ak_ptr;
  bool                            delete_ak;

public:
  Algebraic_kernel_d_1* algebraic_kernel_d_1() const {return _ak_ptr;}

  bool delete_ak_internal_flag() const
  {
    return delete_ak;
  } 
  // Algebraic_kernel_d_1& algebraic_kernel_d_1()             {return _ak;}

public:
  const Cache& cache() const {return _cache;}

public:
  //------------
  //Constructors
  //------------

  //---------------------
  // Default constructor.
  Arr_rational_function_traits_2() : delete_ak(true)
  {
    _ak_ptr = new Algebraic_kernel_d_1;
    _cache.initialize(_ak_ptr);
  }

  Arr_rational_function_traits_2(Algebraic_kernel_d_1* ak_ptr) :
    _ak_ptr(ak_ptr),delete_ak(false)
  {
    _cache.initialize(_ak_ptr);
  }

  Arr_rational_function_traits_2(const Self& other)
    :delete_ak(other.delete_ak_internal_flag())
  {
    //copy kernel
    if (delete_ak)
      _ak_ptr = new Algebraic_kernel_d_1(*other.algebraic_kernel_d_1());
    else
      _ak_ptr = other.algebraic_kernel_d_1();    

    //copy cache
    _cache.initialize(other.cache(), _ak_ptr);
  }

  ~Arr_rational_function_traits_2()
  {
    if (delete_ak)
      delete (_ak_ptr);
  }

  /*! A functor that constructs an x_monotone curve */
  class Construct_x_monotone_curve_2
  {
  protected:
    typedef Arr_rational_function_traits_2<Algebraic_kernel_d_1> Traits;
    typedef CGAL::Arr_rational_arc::Cache<Algebraic_kernel_d_1>  Cache;

    /*! The traits */
    const Traits* _traits;

    /*! Constructor
     * \param traits the traits
     */
    Construct_x_monotone_curve_2(const Traits* traits) : _traits(traits) {}

    friend class Arr_rational_function_traits_2<Algebraic_kernel_d_1>;
    
  public:
    typedef typename Base_rational_arc_ds_1::Polynomial_1 Polynomial_1; 
    typedef typename Base_rational_arc_ds_1::Algebraic_real_1
                                                          Algebraic_real_1;
    typedef Arr_rational_arc::Continuous_rational_arc_d_1<Algebraic_kernel_d_1>
                                                          X_monotone_curve_2;
    typedef Polynomial_1                                  argument_type;
    typedef Polynomial_1                                  first_argument_type;
    typedef Polynomial_1                                  second_argument_type;
    typedef X_monotone_curve_2                            result_type;
    
    X_monotone_curve_2 operator()( const Polynomial_1& P) const
    {
      return X_monotone_curve_2(P, _traits->cache());
    }

    template <typename InputIterator>
    X_monotone_curve_2 operator()( InputIterator begin, InputIterator end) const
    {
      Rat_vector rat_vec(begin,end);
      return X_monotone_curve_2(rat_vec, _traits->cache());
    }

    X_monotone_curve_2 operator()(const Polynomial_1& P,
                                  const Algebraic_real_1& x_s,
                                  bool dir_right) const
    {
      return X_monotone_curve_2(P, x_s, dir_right, _traits->cache());
    }

    template <typename InputIterator>
    X_monotone_curve_2 operator()(InputIterator begin, InputIterator end,
                                  const Algebraic_real_1& x_s,
                                  bool dir_right) const
    {
      Rat_vector rat_vec(begin,end);
      return X_monotone_curve_2(rat_vec, x_s, dir_right, _traits->cache());
    }

    X_monotone_curve_2 operator()(const Polynomial_1& P,
                                  const Algebraic_real_1& x_s,
                                  const Algebraic_real_1& x_t) const
    {
      return X_monotone_curve_2(P, x_s, x_t, _traits->cache());
    }

    template <typename InputIterator>
    X_monotone_curve_2 operator()(InputIterator begin, InputIterator end,
                                  const Algebraic_real_1& x_s,
                                  const Algebraic_real_1& x_t) const
    {
      Rat_vector rat_vec(begin,end);
      return X_monotone_curve_2(rat_vec, x_s, x_t, _traits->cache());
    }

    X_monotone_curve_2 operator()(const Polynomial_1& P,
                                  const Polynomial_1& Q) const 
    {
      return X_monotone_curve_2(P, Q, _traits->cache());
    }

    template <typename InputIterator>
    X_monotone_curve_2 operator()(InputIterator begin_numer,
                                  InputIterator end_numer,
                                  InputIterator begin_denom,
                                  InputIterator end_denom) const 
    {
      Rat_vector rat_vec_numer(begin_numer,end_numer);
      Rat_vector rat_vec_denom(begin_denom,end_denom);
      return X_monotone_curve_2(rat_vec_numer, rat_vec_denom, _traits->cache());
    }

    X_monotone_curve_2 operator()(const Polynomial_1& P, const Polynomial_1& Q,
                                  const Algebraic_real_1& x_s,
                                  bool dir_right) const
    {
      return X_monotone_curve_2(P, Q, x_s, dir_right, _traits->cache());
    }

    template <typename InputIterator>
    X_monotone_curve_2 operator()(InputIterator begin_numer,
                                  InputIterator end_numer,
                                  InputIterator begin_denom,
                                  InputIterator end_denom,
                                  const Algebraic_real_1& x_s,
                                  bool dir_right) const
    {
      Rat_vector rat_vec_numer(begin_numer,end_numer);
      Rat_vector rat_vec_denom(begin_denom,end_denom);
      return X_monotone_curve_2(rat_vec_numer, rat_vec_denom, x_s,dir_right,
                                _traits->cache());
    }

    X_monotone_curve_2 operator()(const Polynomial_1& P,
                                  const Polynomial_1& Q,
                                  const Algebraic_real_1& x_s,
                                  const Algebraic_real_1& x_t) const
    {
      return X_monotone_curve_2(P, Q, x_s, x_t, _traits->cache());
    }

    template <typename InputIterator>
    X_monotone_curve_2 operator()(InputIterator begin_numer,
                                  InputIterator end_numer,
                                  InputIterator begin_denom,
                                  InputIterator end_denom,
                                  const Algebraic_real_1& x_s,
                                  const Algebraic_real_1& x_t) const
    {
      Rat_vector rat_vec_numer(begin_numer, end_numer);
      Rat_vector rat_vec_denom(begin_denom, end_denom);
      return X_monotone_curve_2(rat_vec_numer, rat_vec_denom, x_s, x_t,
                                _traits->cache());
    }
  };

  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const
  {
    return Construct_x_monotone_curve_2(this);
  }

  /*! A functor that constructs an arbitrary curve */
  class Construct_curve_2
  {
  protected:
    typedef Arr_rational_function_traits_2<Algebraic_kernel_d_1> Traits;
    typedef CGAL::Arr_rational_arc::Cache<Algebraic_kernel_d_1>  Cache;

    /*! The traits */
    const Traits* _traits;

    /*! Constructor
     * \param traits the traits
     */
    Construct_curve_2(const Traits* traits) : _traits(traits) {}

    friend class Arr_rational_function_traits_2<Algebraic_kernel_d_1>;
    
  public:
    typedef typename Base_rational_arc_ds_1::Polynomial_1 Polynomial_1; 
    typedef typename Base_rational_arc_ds_1::Algebraic_real_1
                                                          Algebraic_real_1;
    typedef Arr_rational_arc::Rational_arc_d_1<Algebraic_kernel_d_1>
                                                          Curve_2;
    typedef Polynomial_1                                  argument_type;
    typedef Polynomial_1                                  first_argument_type;
    typedef Polynomial_1                                  second_argument_type;
    typedef Curve_2                                       result_type;
    
    Curve_2 operator()(const Polynomial_1& P) const
    {
      return Curve_2(P, _traits->cache());
    }

    template <typename InputIterator>
    Curve_2 operator()(InputIterator begin, InputIterator end) const
    {
      Rat_vector rat_vec(begin, end);
      return Curve_2(rat_vec, _traits->cache());
    }

    Curve_2 operator()(const Polynomial_1& P,
                       const Algebraic_real_1& x_s, bool dir_right) const
    {
      return Curve_2(P, x_s, dir_right, _traits->cache());
    }

    template <typename InputIterator>
    Curve_2 operator()(InputIterator begin, InputIterator end,
                       const Algebraic_real_1& x_s, bool dir_right) const
    {
      Rat_vector rat_vec(begin, end);
      return Curve_2(rat_vec, x_s, dir_right, _traits->cache());
    }

    Curve_2 operator()(const Polynomial_1& P,
                       const Algebraic_real_1& x_s,
                       const Algebraic_real_1& x_t) const
    {
      return Curve_2(P, x_s, x_t, _traits->cache());
    }

    template <typename InputIterator>
    Curve_2 operator()(InputIterator begin, InputIterator end,
                       const Algebraic_real_1& x_s,
                       const Algebraic_real_1& x_t) const
    {
      Rat_vector rat_vec(begin,end);
      return Curve_2(rat_vec, x_s, x_t, _traits->cache());
    }

    Curve_2 operator()(const Polynomial_1& P, const Polynomial_1& Q) const 
    {
      return Curve_2(P, Q, _traits->cache());
    }

    template <typename InputIterator>
    Curve_2 operator()(InputIterator begin_numer, InputIterator end_numer,
                       InputIterator begin_denom, InputIterator end_denom) const 
    {
      Rat_vector rat_vec_numer(begin_numer, end_numer);
      Rat_vector rat_vec_denom(begin_denom, end_denom);
      return Curve_2(rat_vec_numer, rat_vec_denom, _traits->cache());
    }

    Curve_2 operator()(const Polynomial_1& P, const Polynomial_1& Q,
                       const Algebraic_real_1& x_s, bool dir_right) const
    {
      return Curve_2(P, Q, x_s, dir_right, _traits->cache());
    }

    template <typename InputIterator>
    Curve_2 operator()(InputIterator begin_numer, InputIterator end_numer,
                       InputIterator begin_denom, InputIterator end_denom,
                       const Algebraic_real_1& x_s, bool dir_right) const
    {
      Rat_vector rat_vec_numer(begin_numer,end_numer);
      Rat_vector rat_vec_denom(begin_denom,end_denom);
      return Curve_2(rat_vec_numer, rat_vec_denom, x_s, dir_right,
                     _traits->cache());
    }

    Curve_2 operator()(const Polynomial_1& P, const Polynomial_1& Q,
                       const Algebraic_real_1& x_s,
                       const Algebraic_real_1& x_t) const
    {
      return Curve_2(P, Q, x_s, x_t, _traits->cache());
    }

    template <typename InputIterator>
    Curve_2 operator()(InputIterator begin_numer, InputIterator end_numer,
                       InputIterator begin_denom, InputIterator end_denom,
                       const Algebraic_real_1& x_s,
                       const Algebraic_real_1& x_t) const
    {
      Rat_vector rat_vec_numer(begin_numer,end_numer);
      Rat_vector rat_vec_denom(begin_denom,end_denom);
      return Curve_2(rat_vec_numer, rat_vec_denom, x_s, x_t, _traits->cache());
    }
  };

  Construct_curve_2 construct_curve_2_object() const
  {
    return Construct_curve_2(this);
  }

  /*! Construct a point */
  class Construct_point_2
  {
  protected:
    typedef Arr_rational_function_traits_2<Algebraic_kernel_d_1> Traits;
    typedef CGAL::Arr_rational_arc::Cache<Algebraic_kernel_d_1>  Cache;

    /*! The traits */
    const Traits* _traits;

    /*! Constructor
     * \param traits the traits
     */
    Construct_point_2(const Traits* traits) : _traits(traits) {}

    friend class Arr_rational_function_traits_2<Algebraic_kernel_d_1>;
    
  public:
    Point_2 operator()(const Rational_function& rational_function,
                       const Algebraic_real_1& x_coordinate)
    { 
      return Point_2(rational_function, x_coordinate);
    }

    Point_2 operator()(const Rational& x, const Rational& y)
    { 
      Integer  y_numer,y_denom;
      typename FT_rat_1::Decompose()(y,y_numer,y_denom);
      
      return Point_2(_traits->cache().get_rational_function(Rational(y_numer,
                                                                     y_denom)),
                     _traits->algebraic_kernel_d_1()->
                       construct_algebraic_real_1_object()(x));
    }
    Point_2 operator()(const Algebraic_real_1& x, const Rational& y)
    {   
      Integer  y_numer;
      Integer  y_denom;
      typename FT_rat_1::Decompose()(y, y_numer, y_denom);
      return Point_2(_traits->cache().get_rational_function(Rational(y_numer,
                                                                     y_denom)),
                     x);
    }
  }; //Construct_point

  Construct_point_2 construct_point_2_object() const
  {
    return Construct_point_2(this);
  }

//   class Construct_vertical_segment
//   {
//   private:
//     Cache& _cache;

//   public:
//     Construct_vertical_segment(Cache& cache) : _cache(cache) {}

//     Vertical_segment operator()(const Point_2& p) const
//     { 
//       return Vertical_segment(p);
//     }

//     Vertical_segment operator()(const Point_2& p, bool is_directed_up) const
//     { 
//       return Vertical_segment(p, is_directed_up);
//     }

//     Vertical_segment operator()(const Point_2& p1,const Point_2& p2) const
//     {       
//       return Vertical_segment(p1, p2, _cache);
//     }
//   }; //Construct_vertical_segment

//   Construct_vertical_segment construct_vertical_segment_object() const
//   {
//     return Construct_vertical_segment(_cache);
//   }

  //------------------------
  //Functor definitions.
  //------------------------

  //---------------------------------------------------------------
  //A functor that compares the x-coordinates of two points 
  class Compare_x_2
  {
  public:
    /*!
     * Compare the x-coordinates of two points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2);
     *         SMALLER if x(p1) < x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      Comparison_result comp = CGAL::compare(p1.x(), p2.x());
      return (comp);
    }
  };

  /*! Obtain a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object() const
  {
    return Compare_x_2();
  }

  /*! A functor that compares two points lexigoraphically: by x, then by y. */
  class Compare_xy_2
  {
  protected:
    typedef Arr_rational_function_traits_2<Algebraic_kernel_d_1> Traits;
    typedef CGAL::Arr_rational_arc::Cache<Algebraic_kernel_d_1>  Cache;

    /*! The traits */
    const Traits* _traits;

    /*! Constructor
     * \param traits the traits
     */
    Compare_xy_2(const Traits* traits) : _traits(traits) {}

    friend class Arr_rational_function_traits_2<Algebraic_kernel_d_1>;

  public:
    /*!
     * Compares two points lexigoraphically: by x, then by y.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2), or if x(p1) = x(p2) and y(p1) > y(p2);
     *         SMALLER if x(p1) < x(p2), or if x(p1) = x(p2) and y(p1) < y(p2);
     *         EQUAL if the two points are equal.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const 
    {
      return p1.compare_xy_2(p2, _traits->cache());
    }
  };

  /*! Obtain a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object() const
  {
    return Compare_xy_2(this);
  }

  /*! A functor that obtains the left endpoint of a curve. */
  class Construct_min_vertex_2
  {
  public:
    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2 & cv) const
    {
      return (cv.left());
    }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  {
    return Construct_min_vertex_2();
  }

  /*! A functor that obtains the right endpoint of a curve. */
  class Construct_max_vertex_2
  {
  public:
    /*!
     * Get the right endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The right endpoint.
     */
    const Point_2& operator() (const X_monotone_curve_2& cv) const
    {
      return (cv.right());
    }
  };

  /*! Obtain a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  {
    return Construct_max_vertex_2();
  }

  /*! A functor that checks whether a given curve is vertical. */
  class Is_vertical_2
  {
  public:
    /*!
     * Check whether the given x-monotone curve is a vertical segment.
     * \param cv The curve.
     * \return (true) if the curve is a vertical segment; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2&) const
    {
      // A rational function can never be vertical.
      return false;
    }
  };

  /*! Obtain an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object() const
  {
    return Is_vertical_2();
  }

  /*! A functor that compares the y-coordinates of a point and a curve at
   * the point x-coordinate.
   */
  class Compare_y_at_x_2
  {
  private:
    Cache& _cache;
  public:
    Compare_y_at_x_2(Cache& cache) : _cache(cache) {}
    /*!
     * Return the location of the given point with respect to the input curve.
     * \param cv The curve.
     * \param p The point.
     * \pre p is in the x-range of cv.
     * \return SMALLER if y(p) < cv(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& cv) const
    {
      return (cv.point_position(p,_cache));
    }
  };

  /*! Obtain a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object () const
  {
    return Compare_y_at_x_2(_cache);
  }

  /*! A functor that compares compares the y-coordinates of two curves
   * immediately to the left of their intersection point.
   */
  class Compare_y_at_x_left_2
  {
  private:
    Cache& _cache;

  public:
    Compare_y_at_x_left_2(Cache& cache) :_cache(cache) {}
    /*!
     * Compares the y value of two x-monotone curves immediately to the left
     * of their intersection point.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its left.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the left of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  const Point_2& p) const
    {
      // Make sure that p lies on both curves, and that both are defined to its
      // left (so their left endpoint is lexicographically smaller than p).
      CGAL_precondition(cv1.point_position(p,_cache) == EQUAL &&
                        cv2.point_position(p,_cache) == EQUAL);

      CGAL_precondition((cv1.left_parameter_space_in_x() != ARR_INTERIOR ||
                         cv1.left_parameter_space_in_y() != ARR_INTERIOR ||
                         (p.x() > cv1.left().x())) &&
                        (cv2.left_parameter_space_in_x() != ARR_INTERIOR ||
                         cv2.left_parameter_space_in_y() != ARR_INTERIOR ||
                         (p.x() > cv2.left().x())));

      // Compare the two arcs.
      return cv1.compare_at_intersection (cv2,p,true,_cache);}
  };

  /*! Obtain a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  {
    return Compare_y_at_x_left_2(_cache);
  }

  /*! A functor that compares compares the y-coordinates of two curves
   * immediately to the right of their intersection point.
   */
  class Compare_y_at_x_right_2
  {
  private:
    Cache& _cache;

  public:
    Compare_y_at_x_right_2(Cache& cache) :_cache(cache) {}
    /*!
     * Compares the y value of two x-monotone curves immediately to the right
     * of their intersection point.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its right.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the right of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  const Point_2& p) const
    {
      // Make sure that p lies on both curves, and that both are defined to its
      // left (so their left endpoint is lexicographically smaller than p).
      CGAL_precondition(cv1.point_position (p,_cache) == EQUAL &&
                        cv2.point_position (p,_cache) == EQUAL);


      CGAL_precondition((cv1.right_parameter_space_in_x() != ARR_INTERIOR ||
                         cv1.right_parameter_space_in_y() != ARR_INTERIOR ||
                         (p.x() < cv1.right().x())) &&
                        (cv2.right_parameter_space_in_x() != ARR_INTERIOR ||
                         cv2.right_parameter_space_in_y() != ARR_INTERIOR ||
                         (p.x() < cv2.right().x())));

 
      // Compare the two arcs.
      return cv1.compare_at_intersection (cv2,p,false,_cache);
    }
  };

  /*! Obtain a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object () const
  {
    return Compare_y_at_x_right_2(_cache);
  }

  /*! A functor that checks whether two points and two curves are identical. */
  class Equal_2
  {
  protected:
    typedef Arr_rational_function_traits_2<Algebraic_kernel_d_1> Traits;
    typedef CGAL::Arr_rational_arc::Cache<Algebraic_kernel_d_1>  Cache;
    /*! The traits */
    const Traits* _traits;

    /*! Constructor
     * \param traits the traits
     */
    Equal_2(const Traits* traits) : _traits(traits) {}

    friend class Arr_rational_function_traits_2<Algebraic_kernel_d_1>;
    
  public:
    /*!
     * Check if the two x-monotone curves are the same (have the same graph).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2) const
    {
      if (&cv1 == &cv2)
        return true;

      return (cv1.equals(cv2));
    }

    /*!
     * Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const
    {
      if (&p1 == &p2)
        return true;

      return
        (p1.compare_xy_2(p2, _traits->cache()) == CGAL::EQUAL) ?
         true : false;
    }
  };

  /*! Obtain an Equal_2 functor object. */
  Equal_2 equal_2_object() const
  {
    return Equal_2(this);
  }

  /*! A functor that divides a curve into continues (x-monotone) curves. */
  class Make_x_monotone_2
  {
  public:

    /*!
     * Cut the given conic curve (or conic arc) into x-monotone subcurves 
     * and insert them to the given output iterator.
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. The returned
     *           objects is a wrapper for an X_monotone_curve_2 object.
     * \return The past-the-end iterator.
     */
    template<typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const 
    {
      // Make the rational arc continuous.
      std::list<X_monotone_curve_2>                           arcs;

      cv.make_continuous(std::back_inserter(arcs));

      // Create objects.
      typename std::list<X_monotone_curve_2>::const_iterator  iter;

      for (iter = arcs.begin(); iter != arcs.end(); ++iter)
      {
        *oi = make_object (*iter);
        ++oi;
      }

      return (oi);
    }
  };

  /*! Obtain a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object() const
  {
    return Make_x_monotone_2();
  }

  /*! A functor that splits a curve at a point. */
  class Split_2
  {
  private:
    Cache& _cache;

  public:
    Split_2(Cache& cache) : _cache(cache) {}
    /*!
     * Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve (p is its right endpoint).
     * \param c2 Output: The right resulting subcurve (p is its left endpoint).
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator()(const X_monotone_curve_2& cv, const Point_2 & p,
                    X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      cv.split(p, c1, c2, _cache);
    }
  };

  /*! Obtain a Split_2 functor object. */
  Split_2 split_2_object() const
  {
    return Split_2(_cache);
  }

  /*! A functor that computes intersections between two curves. */
  class Intersect_2
  {
  private:
    Cache& _cache;
  public:
    Intersect_2(Cache& cache) : _cache(cache) {}
    /*!
     * Find the intersections of the two given curves and insert them to the
     * given output iterator. As two segments may itersect only once, only a
     * single will be contained in the iterator.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              OutputIterator oi)  const
    {
      return (cv1.intersect (cv2, oi,_cache));
    }
  };

  /*! Obtain an Intersect_2 functor object. */
  Intersect_2 intersect_2_object() const
  {
    return Intersect_2(_cache);
  }

  /*! A functor that tests whether two curves can be merged. */
  class Are_mergeable_2
  {
  public:
    /*!
     * Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are mergeable - if they are supported
     *         by the same line and share a common endpoint; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const
    {
      return (cv1.can_merge_with(cv2));
    }
  };

  /*! Obtain an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object() const
  {
    return Are_mergeable_2();
  }

  /*! \class Merge_2
   * A functor that merges two x-monotone arcs into one.
   */
  class Merge_2
  {
  protected:
    typedef Arr_rational_function_traits_2<Algebraic_kernel_d_1>        Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;
    
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Merge_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_rational_function_traits_2<Algebraic_kernel_d_1>;

  public:
    /*!
     * Merge two given x-monotone curves into a single curve (segment).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable.
     */
    void operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2,
                    X_monotone_curve_2& c) const
    {
      CGAL_precondition(m_traits->are_mergeable_2_object()(cv2, cv1));

      c = cv1;
      c.merge(cv2);
    }
  };

  /*! Obtain a Merge_2 functor object. */
  Merge_2 merge_2_object() const
  {
    return Merge_2(this);
  }

  //@}

  /// \name Functor definitions to handle boundaries
  //@{

  /*! A function object that obtains the parameter space of a geometric
   * entity along the x-axis
   */
  class Parameter_space_in_x_2 {
  public:
    /*! Obtains the parameter space at the end of a line along the x-axis.
     * \param xcv the line
     * \param ce the line end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the line xcv.
     *   ARR_LEFT_BOUNDARY  - the line approaches the identification arc from
     *                        the right at the line left end.
     *   ARR_INTERIOR       - the line does not approache the identification arc.
     *   ARR_RIGHT_BOUNDARY - the line approaches the identification arc from
     *                        the left at the line right end.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & xcv,
        Arr_curve_end ce) const
    {
      return (ce == ARR_MIN_END) ?
        xcv.left_parameter_space_in_x() : xcv.right_parameter_space_in_x();
    }

    /*! Obtains the parameter space at a point along the x-axis.
     * \param p the point.
     * \return the parameter space at p.
     */
    Arr_parameter_space operator()(const Point_2 ) const
    {
      return ARR_INTERIOR;
    }
  };

  /*! Obtain a Parameter_space_in_x_2 function object */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  { return Parameter_space_in_x_2(); }
  
  /*! A function object that obtains the parameter space of a geometric
   * entity along the y-axis
   */
  class Parameter_space_in_y_2 {
  public:
    /*! Obtains the parameter space at the end of a line along the y-axis .
     * Note that if the line end coincides with a pole, then unless the line
     * coincides with the identification arc, the line end is considered to
     * be approaching the boundary, but not on the boundary.
     * If the line coincides with the identification arc, it is assumed to
     * be smaller than any other object.
     * \param xcv the line
     * \param ce the line end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the line xcv.
     *   ARR_BOTTOM_BOUNDARY  - the line approaches the south pole at the line
     *                          left end.
     *   ARR_INTERIOR         - the line does not approache a contraction point.
     *   ARR_TOP_BOUNDARY     - the line approaches the north pole at the line
     *                          right end.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & xcv,
        Arr_curve_end ce) const
    {
      return (ce == ARR_MIN_END) ?
        xcv.left_parameter_space_in_y() : xcv.right_parameter_space_in_y();
    }

    /*! Obtains the parameter space at a point along the y-axis.
     * \param p the point.
     * \return the parameter space at p.
     */
    Arr_parameter_space operator()(const Point_2 ) const
    {
      return ARR_INTERIOR;
    }
  };

  /*! Obtain a Parameter_space_in_y_2 function object */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(); }

  /*! A function object that compares the x-coordinates of arc ends near the
   * boundary of the parameter space
   */
  class Compare_x_near_boundary_2 {
  public:
    /*! Compare the x-coordinate of a point with the x-coordinate of
     * a line end near the boundary at y = +/- oo.
     * \param p the point direction.
     * \param xcv the line, the endpoint of which is compared.
     * \param ce the line-end indicator -
     *            ARR_MIN_END - the minimal end of xc or
     *            ARR_MAX_END - the maximal end of xc.
     * \return the comparison result:
     *         SMALLER - x(p) < x(xc, ce);
     *         EQUAL   - x(p) = x(xc, ce);
     *         LARGER  - x(p) > x(xc, ce).     
     * \pre p lies in the interior of the parameter space.
     * \pre the ce end of the line xcv lies on a boundary.
     */
    Comparison_result operator()(const Point_2 & p,
                                 const X_monotone_curve_2 & xcv,
                                 Arr_curve_end ce) const
    {
      Comparison_result r = xcv.compare_end(ce, p);
      if (r == EQUAL) 
        return EQUAL; 
      return (r == NEGATIVE) ? POSITIVE : NEGATIVE ; 
    }

    /*! Compare the x-coordinates of 2 arcs ends near the boundary of the
     * parameter space at y = +/- oo.
     * \param xcv1 the first arc.
     * \param ce1 the first arc end indicator -
     *            ARR_MIN_END - the minimal end of xcv1 or
     *            ARR_MAX_END - the maximal end of xcv1.
     * \param xcv2 the second arc.
     * \param ce2 the second arc end indicator -
     *            ARR_MIN_END - the minimal end of xcv2 or
     *            ARR_MAX_END - the maximal end of xcv2.
     * \return the second comparison result:
     *         SMALLER - x(xcv1, ce1) < x(xcv2, ce2);
     *         EQUAL   - x(xcv1, ce1) = x(xcv2, ce2);
     *         LARGER  - x(xcv1, ce1) > x(xcv2, ce2).
     * \pre the ce1 end of the line xcv1 lies on a boundary.
     * \pre the ce2 end of the line xcv2 lies on a boundary.
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce2) const
    {
      return xcv1.compare_ends(ce1, xcv2, ce2);
    }
  };

  /*! Obtain a Compare_x_near_boundary_2 function object */
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const
  { return Compare_x_near_boundary_2(); }
    

  /*! A function object that compares the y-coordinates of arc ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 
  {
  private:
    Cache& _cache;

  public:
    /*! Compare the y-coordinates of 2 lines at their ends near the boundary
     * of the parameter space at x = +/- oo.
     * \param xcv1 the first arc.
     * \param xcv2 the second arc.
     * \param ce the line end indicator.
     * \return the second comparison result.
     * \pre the ce ends of the lines xcv1 and xcv2 lie either on the left
     * boundary or on the right boundary of the parameter space.
     */
    Compare_y_near_boundary_2(Cache& cache) : _cache(cache) {}
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
        const X_monotone_curve_2 & xcv2,
        Arr_curve_end ce) const
    {
      return (ce == ARR_MIN_END) ?
        xcv1.compare_at_minus_infinity(xcv2,_cache) :
        xcv1.compare_at_plus_infinity(xcv2,_cache);
    }
  };

  /*! Obtain a Compare_y_near_boundary_2 function object */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
  { return Compare_y_near_boundary_2(_cache); }

 
  /*! A function object that compares at limit
   */
  //new functor
  class Compare_x_at_limit_2
  {
   public:
    /*! Compares the x coordinate of p with the curve end
     * of xcv that is defined by ce at its limit. 
     * Returns SMALLER, EQUAL, or LARGER accordingly.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2&  xcv, 
                                 Arr_curve_end ce)
    {
      CGAL_precondition(Parameter_space_in_x_2()(xcv,ce) == ARR_INTERIOR);
      CGAL_precondition(Parameter_space_in_y_2()(xcv,ce) != ARR_INTERIOR);
      return CGAL::compare(p.x(),
                           (ce == ARR_MIN_END) ? xcv.left_x() : xcv.right_x());
    }
    /*! Compares the curve end of  xcv1 that is defined by ce1 
     *  with the curve end of xcv2 that is defined by ce2
     * at their limits in x. 
     * Returns SMALLER, EQUAL, or LARGER accordingly.
     */
    Comparison_result operator()(const X_monotone_curve_2&  xcv1, 
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2&  xcv2, 
                                 Arr_curve_end ce2)
    {
      CGAL_precondition(Parameter_space_in_x_2()(xcv1,ce1) == ARR_INTERIOR);
      CGAL_precondition(Parameter_space_in_y_2()(xcv1,ce1) != ARR_INTERIOR);
      CGAL_precondition(Parameter_space_in_x_2()(xcv2,ce2) == ARR_INTERIOR);
      CGAL_precondition(Parameter_space_in_y_2()(xcv2,ce2) != ARR_INTERIOR);

      return CGAL::compare((ce1 == ARR_MIN_END) ? xcv1.left_x() : xcv1.right_x(),
                           (ce2 == ARR_MIN_END) ? xcv2.left_x() : xcv2.right_x());
    }

  }; //Compare_x_at_limit_2

  /*! Obtain a Compare_x_at_limit_2 function object */
  Compare_x_at_limit_2 compare_x_at_limit_2_object() const
  { return Compare_x_at_limit_2(); }
  //@}
  
  /// \name Functor definitions for the Boolean set-operation traits.
  //@{
 
  //new functor
  class Compare_x_near_limit_2
  {
  private:
    Cache& _cache;

  public:
    Compare_x_near_limit_2(Cache& cache) : _cache(cache) {}
    /*! Compares the curve end of  xcv1 that is defined by ce1 
     *  with the curve end of xcv2 that is defined by ce2
     * at their limits in x. 
     * Returns SMALLER, EQUAL, or LARGER accordingly.
     */
    Comparison_result operator()( const X_monotone_curve_2& xcv1, 
                                  const X_monotone_curve_2& xcv2, 
                                  Arr_curve_end ce) const
    {
      return xcv1.compare_near_end(xcv2,ce,_cache);
    }
  }; //Compare_x_near_limit_2

  /*! Obtain a Compare_x_near_limit_2 function object */
  Compare_x_near_limit_2 compare_x_near_limit_2_object() const
  { return Compare_x_near_limit_2(_cache); }

  class Compare_endpoints_xy_2
  {
  public:

    /*!
     * Compare the endpoints of an $x$-monotone curve lexicographically.
     * (assuming the curve has a designated source and target points).
     * \param cv The curve.
     * \return SMALLER if the curve is directed right;
     *         LARGER if the curve is directed left.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv)
    {
      if (cv.is_directed_right())
        return (SMALLER);
      else
        return (LARGER);
    }
  };

  /*! Obtain a Compare_endpoints_xy_2 functor object. */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  {
    return Compare_endpoints_xy_2();
  }

  class Construct_opposite_2
  {
  public:

    /*!
     * Construct an opposite x-monotone (with swapped source and target).
     * \param cv The curve.
     * \return The opposite curve.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& cv)
    {
      return (cv.flip());
    }
  };

  /*! Obtain a Construct_opposite_2 functor object. */
  Construct_opposite_2 construct_opposite_2_object() const
  {
    return Construct_opposite_2();
  }

  //@}

  class Approximate_2{
    Approximate_number_type approx_x(const Point_2& p){
      return Approximate_number_type(p.x().lower());
    } 
    Approximate_number_type approx_y(const Point_2& p){
      typedef typename Algebraic_kernel_d_1::Polynomial_1 Polynomial_1;
      typename CGAL::Coercion_traits<Polynomial_1,Bound>::Cast cast;  
      return
        cast(p.rational_function().numer()).evaluate(p.x().lower())/
        cast(p.rational_function().denom()).evaluate(p.x().lower());
    }
  public:
    Approximate_number_type operator()(const Point_2& p, int i){
      if(i==0) return approx_x(p); 
      if(i==1) return approx_y(p);
      CGAL_assertion(false);
      return Approximate_number_type(0);
    }
  };
  
  Approximate_2 approximate_2_object() const { return Approximate_2(); }

  void cleanup_cache() const
  {
    _cache.cleanup();
  }
}; //Arr_rational_function_traits_2

}   //namespace CGAL {

#include <CGAL/enable_warnings.h>

#endif  //CGAL_ARR_RATIONAL_ARC_TRAITS_D_1_H
