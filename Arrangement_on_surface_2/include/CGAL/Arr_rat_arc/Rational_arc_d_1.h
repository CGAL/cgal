// Copyright (c) 2011 Tel-Aviv University (Israel), INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Oren Salzman <orenzalz@post.tau.ac.il >
//                 Michael Hemmer <Michael.Hemmer@sophia.inria.fr>


#ifndef CGAL_RATIONAL_ARC_D_1_H
#define CGAL_RATIONAL_ARC_D_1_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <vector>
#include <list>
#include <ostream>
#include <CGAL/Arr_enums.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>

#include <CGAL/Fraction_traits.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d_1.h>

#include <CGAL/Arr_rat_arc/Base_rational_arc_ds_1.h>
#include <CGAL/Arr_rat_arc/Algebraic_point_2.h>
#include <CGAL/Arr_rat_arc/Cache.h>
#include <CGAL/Arr_rat_arc/Rational_function.h>
#include <CGAL/Arr_rat_arc/Rational_function_pair.h>

#include <boost/type_traits/is_same.hpp>

namespace CGAL {
namespace Arr_rational_arc {

//--------------------------------------------------------------------------//
//   class Base_rational_arc_d_1
//Representation of an segment of a rational function, given as:
//
//        Numerator(x)
//   y = -------------              x_l <= x <= x_r
//        Denominator(x)
//
// where Numerator and Denominator are polynomial with integer (or rational) coefficients.
// The class is templated with two parameters:
// Algebraic_kernel: An algebraic kernel for the intersection points of the curves
//
// This class serves as the base for the classes:
// Rational_arc_d_1 (a general, not necessarily continuous arc)
// Continuous_rational_arc_d_1 (a continuous portion of a rational function).
//--------------------------------------------------------------------------//

template <typename Algebraic_kernel_>
class Base_rational_arc_d_1
{
public:
  typedef Algebraic_kernel_                                 Algebraic_kernel;
  typedef Base_rational_arc_d_1<Algebraic_kernel>           Self;

  typedef CGAL::Arr_rational_arc::Base_rational_arc_ds_1<Algebraic_kernel>
    Base_rational_arc_ds_1;
  typedef CGAL::Arr_rational_arc::Rational_function<Algebraic_kernel>
                                                            Rational_function;
  typedef CGAL::Arr_rational_arc::Rational_function_pair<Algebraic_kernel>
    Rational_function_pair;
  typedef CGAL::Arr_rational_arc::Algebraic_point_2<Algebraic_kernel>
                                                            Algebraic_point_2;

  typedef CGAL::Arr_rational_arc::Cache<Algebraic_kernel>   Cache;

  typedef typename Base_rational_arc_ds_1::Multiplicity     Multiplicity;
  typedef typename Base_rational_arc_ds_1::Polynomial_1     Polynomial_1;
  typedef typename Base_rational_arc_ds_1::Coefficient      Coefficient;
  typedef typename Base_rational_arc_ds_1::Arithmetic_kernel
                                                            Arithmetic_kernel;
  typedef typename Base_rational_arc_ds_1::Rational         Rational;
  typedef typename Base_rational_arc_ds_1::Integer          Integer;
  typedef typename Base_rational_arc_ds_1::Algebraic_real_1 Algebraic_real_1;
  typedef typename Base_rational_arc_ds_1::Algebraic_vector Algebraic_vector;
  typedef typename Base_rational_arc_ds_1::Multiplicity_vector
                                                            Multiplicity_vector;

  typedef std::vector<Rational>                             Rat_vector;

  typedef Polynomial_traits_d<Polynomial_1>                 Polynomial_traits_1;
  typedef typename Base_rational_arc_ds_1::FT_rat_1         FT_rat_1;
  typedef typename Base_rational_arc_ds_1::Solve_1          Solve_1;
  typedef typename Algebraic_kernel::Bound                  Bound;
  typedef Algebraic_structure_traits<Polynomial_1>          AT_poly;

  typedef Polynomial<Rational>                              Poly_rat_1;
  typedef Polynomial_traits_d<Poly_rat_1>                   PT_rat_1;
  typedef Fraction_traits <Poly_rat_1>                      FT_poly_rat_1;

  typedef Algebraic_point_2                                 Point_2;

  CGAL_static_assertion((boost::is_same<Integer, Coefficient>::value));
  CGAL_static_assertion((boost::is_same<Polynomial_1,
                       typename FT_poly_rat_1::Numerator_type>::value));

public:
  const Rational_function& get_rational_function(const Polynomial_1& numerator,
                                                 const Polynomial_1& denominator,
                                                 const Cache& cache) const
  {
    return cache.get_rational_function(numerator, denominator);
  }
  const Rational_function& get_rational_function(const Rational& rat,
                                                 const Cache& cache) const
  {
    return cache.get_rational_function(rat);
  }
  const Rational_function_pair get_rational_pair(const Rational_function& f,
                                                 const Rational_function& g,
                                                 const Cache& cache) const
  {
    CGAL_precondition(f.id() != g.id());
    return cache.get_rational_pair(f, g);
  }

public:
  //------------------------------
  //Base_rational_arc_d_1 members
  //------------------------------
  enum
  {
    SRC_AT_X_MINUS_INFTY = 1,
    SRC_AT_X_PLUS_INFTY = 2,
    SRC_AT_Y_MINUS_INFTY = 4,
    SRC_AT_Y_PLUS_INFTY = 8,

    SRC_INFO_BITS = SRC_AT_X_MINUS_INFTY + SRC_AT_X_PLUS_INFTY +
    SRC_AT_Y_MINUS_INFTY + SRC_AT_Y_PLUS_INFTY,

    TRG_AT_X_MINUS_INFTY = 16,
    TRG_AT_X_PLUS_INFTY = 32,
    TRG_AT_Y_MINUS_INFTY = 64,
    TRG_AT_Y_PLUS_INFTY = 128,

    TRG_INFO_BITS = TRG_AT_X_MINUS_INFTY + TRG_AT_X_PLUS_INFTY +
    TRG_AT_Y_MINUS_INFTY + TRG_AT_Y_PLUS_INFTY,

    IS_DIRECTED_RIGHT = 256,
    IS_CONTINUOUS = 512,
    IS_VALID = 1024
  };

  Rational_function _f;    // The rational function
  Algebraic_point_2   _ps;   // The source point.
  Algebraic_point_2   _pt;   // The target point.
  int     _info;   // A set of Boolean flags.

public:
  //------------
  //Constructors
  //------------

  //---------------------------------------------------------------------------
  //default constructor
  Base_rational_arc_d_1() :
    _info(0)
  {}

  //---------------------------------------------------------------------------
  // Constructor of a whole polynomial curve defined by pcoeffs - the rational
  // coefficients of the polynomial p(x).

  Base_rational_arc_d_1(const Polynomial_1& P, const Cache& cache) :
    _info(0)
  {
    _init(P, Integer(1), cache);
  }

  Base_rational_arc_d_1(const Rat_vector& pcoeffs, const Cache& cache) :
    _info(0)
  {
    // Set the numerator & denominator polynomials.
    Polynomial_1 _numer;
    Poly_rat_1 num_rat(typename PT_rat_1::Construct_polynomial()(pcoeffs.begin(),
                                                                 pcoeffs.end()));
    Integer denom_int;
    typename FT_poly_rat_1::Decompose()(num_rat, _numer, denom_int);

    _init(_numer,denom_int,cache);
  }

  void _init(const Polynomial_1& P,const Integer& Q_int,const Cache& cache)
  {
    CGAL_precondition(!CGAL::is_zero(Q_int));
    //set rational function
    Polynomial_1 Q= typename Polynomial_traits_1::Construct_polynomial()(Q_int);
    _f = get_rational_function(P, Q, cache);

    //Mark that the endpoints of the polynomial are unbounded
    //(the source is at x = -oo and the target is at x = +oo).
    _info = (_info | SRC_AT_X_MINUS_INFTY);
    _info = (_info | TRG_AT_X_PLUS_INFTY);
    _info = (_info | IS_DIRECTED_RIGHT);

    // Check whether the end points lie at y = -oo or at y = +oo.
    const int    deg_num(CGAL::degree(P));
    Integer lead_coeff(CGAL::leading_coefficient(Q));
    CGAL::Sign   lead_sign(CGAL::sign(lead_coeff));

    if (deg_num > 0)
    {
      // Check if the degree is even or odd and check the sign of the leading
      // coefficient of the polynomial.

      CGAL_assertion(lead_sign != CGAL::ZERO);

      if (deg_num % 2 == 0)
      {
        // Polynomial of an even degree.
        if (lead_sign == CGAL::NEGATIVE)
          _info = (_info | SRC_AT_Y_MINUS_INFTY | TRG_AT_Y_MINUS_INFTY);
        else
          _info = (_info | SRC_AT_Y_PLUS_INFTY | TRG_AT_Y_PLUS_INFTY);
      }
      else
      {
        // Polynomial of an odd degree.
        if (lead_sign == CGAL::NEGATIVE)
          _info = (_info | SRC_AT_Y_PLUS_INFTY | TRG_AT_Y_MINUS_INFTY);
        else
          _info = (_info | SRC_AT_Y_MINUS_INFTY | TRG_AT_Y_PLUS_INFTY);
      }
    }
    else
    {
      // In the case of a constant polynomial it is possible to set a finite
      // y-coordinate for the source and target points.
      //x coordinate is 0 although in practice is +-oo
      _ps = Algebraic_point_2(_f,Algebraic_real_1());
      //x coordinate is 0 although in practice is +-oo
      _pt = Algebraic_point_2(_f,Algebraic_real_1());
    }

    // Mark that the arc is continuous and valid.
    _info = (_info | IS_CONTINUOUS);
    _info = (_info | IS_VALID);

  }

  //---------------------------------------------------------------------------
  //Constructor of a polynomial ray, defined by y = p(x),
  //for x_s <= x if the ray is directed to the right, or
  //for x_s >= x if it is directed to the left.
  //param pcoeffs The rational coefficients of the polynomial p(x).
  //param x_s The x-coordinate of the source point.
  //param dir_right Is the ray directed to the right (to +oo) or to the left
  //(to -oo).

  Base_rational_arc_d_1(const Polynomial_1& P, const Algebraic_real_1& x_s,
                        bool dir_right, const Cache& cache) :
    _info(0)
  {
    _init(P, Polynomial_1(1), x_s, dir_right, cache);
  }

  Base_rational_arc_d_1(const Rat_vector& pcoeffs,const Algebraic_real_1& x_s,
                        bool dir_right, const Cache& cache) :
    _info(0)
  {
    // Set the numerator & denominator polynomials.
    Polynomial_1 _numer;
    Poly_rat_1 num_rat(typename PT_rat_1::Construct_polynomial()(pcoeffs.begin(),
                                                                 pcoeffs.end()));
    Integer denom_int;
    typename FT_poly_rat_1::Decompose()(num_rat, _numer, denom_int);

    _init(_numer, denom_int, x_s, dir_right, cache);
  }

  void _init(const Polynomial_1& P,const Integer& Q_int,
             const Algebraic_real_1& x_s, bool dir_right,const Cache& cache)
  {
    CGAL_precondition(!CGAL::is_zero(Q_int));
    //set rational function
    Polynomial_1 Q = typename Polynomial_traits_1::Construct_polynomial()(Q_int);
    _f = get_rational_function(P,Q,cache);

    // Mark that the target points of the polynomial is unbounded.
    if (dir_right)
    {
      _info = (_info | TRG_AT_X_PLUS_INFTY);
      _info = (_info | IS_DIRECTED_RIGHT);
    }
    else
    {
      _info = (_info | TRG_AT_X_MINUS_INFTY);
    }

    // Set the source point.
    _ps=Algebraic_point_2(_f,x_s);

    // Check whether the target point lies at y = -oo or at y = +oo.
    const int   deg_num(CGAL::degree(P));
    Integer     lead_coeff(CGAL::leading_coefficient(P));
    CGAL::Sign  lead_sign(CGAL::sign(lead_coeff));

    if (deg_num > 0)
    {
      //Check if the degree is even or odd and check the sign of the leading
      // coefficient of the polynomial.
      CGAL_assertion(lead_sign != CGAL::ZERO);

      if (dir_right)
      {
        // The target is at x= +oo, thus:
        if (lead_sign == CGAL::POSITIVE)
          _info = (_info | TRG_AT_Y_PLUS_INFTY);
        else
          _info = (_info | TRG_AT_Y_MINUS_INFTY);
      }
      else
      {
        // The target is at x= -oo, thus:
        if ((deg_num % 2 == 0 && lead_sign == CGAL::POSITIVE) ||
            (deg_num % 2 == 1 && lead_sign == CGAL::NEGATIVE))
          _info = (_info | TRG_AT_Y_PLUS_INFTY);
        else
          _info = (_info | TRG_AT_Y_MINUS_INFTY);
      }
    }
    else
    {
      // In the case of a constant polynomial it is possible to set a finite
      // y-coordinate for the target point.
      // x coordinate is 0 although in practice is +-oo
      _pt = Algebraic_point_2(get_rational_function(P, Q, cache),
                              Algebraic_real_1());
    }

    // Mark that the arc is continuous and valid.
    _info = (_info | IS_CONTINUOUS);
    _info = (_info | IS_VALID);
  }

  //---------------------------------------------------------------------------
  //Constructor of a polynomial arc, defined by y = p(x), x_min <= x <= x_max.
  //for x_s <= x if the ray is directed to the right, or
  //for x_s >= x if it is directed to the left.
  //param pcoeffs The rational coefficients of the polynomial p(x).
  //param x_s The x-coordinate of the source point.
  //param x_t The x-coordinate of the target point.
  //precondition: The two x-coordinates must not be equal.
  Base_rational_arc_d_1(const Polynomial_1& P,
                        const Algebraic_real_1& x_s, const Algebraic_real_1& x_t,
                        const Cache& cache):
    _info(0)
  {
    _init(P,Integer(1), x_s, x_t, cache);
  }
  Base_rational_arc_d_1(const Rat_vector& pcoeffs,
                        const Algebraic_real_1& x_s,const Algebraic_real_1& x_t,
                        const Cache& cache):
    _info(0)
  {
    // Set the numerator & denominator polynomials.
    Polynomial_1 _numer;
    Poly_rat_1 num_rat(typename PT_rat_1::Construct_polynomial()(pcoeffs.begin(),
                                                                 pcoeffs.end()));
    Integer denom_int;
    typename FT_poly_rat_1::Decompose()(num_rat, _numer, denom_int);

    _init(_numer, denom_int, x_s, x_t, cache);
  }

  void _init(const Polynomial_1& P,const Integer& Q_int,
             const Algebraic_real_1& x_s,const Algebraic_real_1& x_t,
             const Cache& cache)
  {
    CGAL_precondition(!CGAL::is_zero(Q_int));
    //set rational function
    Polynomial_1 Q = typename Polynomial_traits_1::Construct_polynomial()(Q_int);
    _f = get_rational_function(P, Q, cache);

    // Compare the x-coordinates and determine the direction.
    Comparison_result   x_res = CGAL::compare(x_s, x_t);

    CGAL_precondition(x_res != EQUAL);

    if (x_res == SMALLER)
      _info = (_info | IS_DIRECTED_RIGHT);

    // Set the endpoints.
    _ps=Algebraic_point_2(_f,x_s);
    _pt=Algebraic_point_2(_f,x_t);

    // Mark that the arc is continuous and valid.
    _info = (_info | IS_CONTINUOUS);
    _info = (_info | IS_VALID);
    return;
  }

  //---------------------------------------------------------------------------
  //Constructor of a polynomial function, defined by y = p(x)/q(x) for any x.
  //param pcoeffs The rational coefficients of the polynomial p(x).
  //param qcoeffs The rational coefficients of the polynomial q(x).
  Base_rational_arc_d_1(const Polynomial_1& P, const Polynomial_1& Q,
                        const Cache& cache) :
    _info(0)
  {
    _init(P,Q,cache);
  }

  Base_rational_arc_d_1(const Rat_vector& pcoeffs, const Rat_vector& qcoeffs,
                        const Cache& cache) :
    _info(0)
  {
    Polynomial_1 _numer;
    Polynomial_1 _denom;
    Poly_rat_1 numer_rat(typename
                         PT_rat_1::Construct_polynomial()(pcoeffs.begin(),
                                                          pcoeffs.end()));
    Poly_rat_1 denom_rat(typename
                         PT_rat_1::Construct_polynomial()(qcoeffs.begin(),
                                                          qcoeffs.end()));
    Integer denom_numer_int,denom_denom_int;
    typename FT_poly_rat_1::Decompose()(numer_rat, _numer, denom_numer_int);
    typename FT_poly_rat_1::Decompose()(denom_rat, _denom, denom_denom_int);
    _numer *= denom_denom_int;
    _denom *= denom_numer_int;

    _init(_numer,_denom,cache);
  }
  void _init(const Polynomial_1& P_, const Polynomial_1& Q_, const Cache& cache)
  {
    CGAL_precondition(!CGAL::is_zero(Q_));
    //set rational function
    // Set the numerator & denominator polynomials.

    Polynomial_1 P;
    Polynomial_1 Q;
    _canonicalize(P_, Q_, P, Q);

    _f = get_rational_function(P, Q, cache);

    // Mark that the endpoints of the rational functions are unbounded (the
    // source is at x = -oo and the target is at x = +oo).
    _info = (_info | SRC_AT_X_MINUS_INFTY);
    _info = (_info | TRG_AT_X_PLUS_INFTY);
    _info = (_info | IS_DIRECTED_RIGHT);


    // Analyze the bahaviour of the rational function at x = -oo (the source).
    Algebraic_real_1          y0;
    const Arr_parameter_space inf_s = _analyze_at_minus_infinity(P, Q, y0);

    if (inf_s == ARR_BOTTOM_BOUNDARY)
      _info = (_info | SRC_AT_Y_MINUS_INFTY);
    else if (inf_s == ARR_TOP_BOUNDARY)
      _info = (_info | SRC_AT_Y_PLUS_INFTY);
    else // if (inf_s == ARR_INTERIOR)
      _ps = Algebraic_point_2();   //the point is a dummy
    //Analyze the bahaviour of the rational function at x = +oo (the target).
    const Arr_parameter_space inf_t = _analyze_at_plus_infinity(P, Q, y0);

    if (inf_t == ARR_BOTTOM_BOUNDARY)
      _info = (_info | TRG_AT_Y_MINUS_INFTY);
    else if (inf_t == ARR_TOP_BOUNDARY)
      _info = (_info | TRG_AT_Y_PLUS_INFTY);
    else // if (inf_t == ARR_INTERIOR)
      _pt =  Algebraic_point_2();   //the point is a dummy

    // Mark that the arc is valid. As it may have poles, we mark it
    // as continuous only if the denominator has no roots.
    _info = ( _info | ( this->_is_continuous() ?
            (IS_CONTINUOUS | IS_VALID) : IS_VALID ) );
  }

  //---------------------------------------------------------------------------
  //Constructor of a ray of a rational function, defined by y = p(x)/q(x),
  //for x_s <= x if the ray is directed to the right, or
  //for x_s >= x if the ray is directed to the left.
  //param pcoeffs The rational coefficients of the polynomial p(x).
  //param qcoeffs The rational coefficients of the polynomial q(x).
  //param x_s The x-coordinate of the source point.
  //param dir_right Is the ray directed to the right (to +oo) or to the left
  //(to -oo).
  Base_rational_arc_d_1(const Polynomial_1& P, const Polynomial_1& Q,
                        const Algebraic_real_1& x_s, bool dir_right,
                        const Cache& cache) :
    _info(0)
  {
    _init(P, Q, x_s, dir_right, cache);
  }
  Base_rational_arc_d_1(const Rat_vector& pcoeffs, const Rat_vector& qcoeffs,
                        const Algebraic_real_1& x_s, bool dir_right,
                        const Cache& cache) :
    _info(0)
  {
    // Set the numerator and denominator polynomials.
    Polynomial_1 _numer;
    Polynomial_1 _denom;
    Poly_rat_1 numer_rat(typename
                         PT_rat_1::Construct_polynomial()(pcoeffs.begin(),
                                                          pcoeffs.end()));
    Poly_rat_1 denom_rat(typename
                         PT_rat_1::Construct_polynomial()(qcoeffs.begin(),
                                                          qcoeffs.end()));
    Integer denom_numer_int,denom_denom_int;
    typename FT_poly_rat_1::Decompose()(numer_rat, _numer, denom_numer_int);
    typename FT_poly_rat_1::Decompose()(denom_rat, _denom, denom_denom_int);
    _numer *= denom_denom_int;
    _denom *= denom_numer_int;

    _init(_numer,_denom, x_s, dir_right, cache);
  }

  void _init(const Polynomial_1& P_, const Polynomial_1& Q_,
             const Algebraic_real_1& x_s, bool dir_right,
             const Cache& cache)
  {
    CGAL_precondition(!CGAL::is_zero(Q_));
    //set rational function
    Polynomial_1 P;
    Polynomial_1 Q;
    _canonicalize(P_,Q_,P,Q);
    _f = get_rational_function(P, Q, cache);

    // Mark that the target points of the polynomial is unbounded.
    if (dir_right)
    {
      _info = (_info | TRG_AT_X_PLUS_INFTY);
      _info = (_info | IS_DIRECTED_RIGHT);
    }
    else
    {
      _info = (_info | TRG_AT_X_MINUS_INFTY);
    }


    //The source point has a bounded x-coordinate.
    _ps = Algebraic_point_2(_f, x_s);
    //check if the source point lies next to a pole.
    if (typename Algebraic_kernel::Sign_at_1()(Q, x_s) != CGAL::ZERO)
    {
      // We have a nomral endpoint.
      //nothing to do....
    }
    else
    {
      // The y-coodinate is unbounded, but we can set its sign.
      std::pair<CGAL::Sign, CGAL::Sign>  signs = _analyze_near_pole(x_s);
      const CGAL::Sign sign_s = (dir_right ? signs.second : signs.first);

      _info = (sign_s == CGAL::NEGATIVE) ?
        (_info | SRC_AT_Y_MINUS_INFTY) : (_info | SRC_AT_Y_PLUS_INFTY);
    }

    // Set the properties of the target.
    Algebraic_real_1    y0;
    Arr_parameter_space inf_t;
    if (dir_right)
    {
      // The target point is at x = +oo.
      inf_t=_analyze_at_plus_infinity(P, Q, y0);
    }
    else
    {
      // The target point is at x = -oo.
      inf_t =_analyze_at_minus_infinity(P, Q, y0);
    }

    if (inf_t == ARR_BOTTOM_BOUNDARY)
      _info = (_info | TRG_AT_Y_MINUS_INFTY);
    else if (inf_t == ARR_TOP_BOUNDARY)
      _info = (_info | TRG_AT_Y_PLUS_INFTY);
    else // if (inf_t == ARR_INTERIOR)
      _pt = Algebraic_point_2( ); //the point is a dummy

    // Mark that the arc is valid. As it may have poles, we mark it
    // as continuous only if the denominator has no roots.
    _info = ( _info | ( this->_is_continuous() ?
            (IS_CONTINUOUS | IS_VALID) : IS_VALID ) );
  }
  //---------------------------------------------------------------------------
  //Constructor of a ray of a rational function, defined by y = p(x)/q(x),
  //where: x_min <= x <= x_max
  //param pcoeffs The rational coefficients of the polynomial p(x).
  //param qcoeffs The rational coefficients of the polynomial q(x).
  //param x_s The x-coordinate of the source point.
  //param x_t The x-coordinate of the target point.
  //precondition: The two x-coordinates must not be equal.
  Base_rational_arc_d_1(const Polynomial_1& P, const Polynomial_1& Q,
                        const Algebraic_real_1& x_s, const Algebraic_real_1& x_t,
                        const Cache& cache):
    _info(0)
  {
    _init(P, Q, x_s, x_t, cache);
  }

  Base_rational_arc_d_1(const Rat_vector& pcoeffs, const Rat_vector& qcoeffs,
                        const Algebraic_real_1& x_s, const Algebraic_real_1& x_t,
                        const Cache& cache):
    _info(0)
  {
    Polynomial_1 _numer;
    Polynomial_1 _denom;
    Poly_rat_1 numer_rat(typename
                         PT_rat_1::Construct_polynomial()(pcoeffs.begin(),
                                                          pcoeffs.end()));
    Poly_rat_1 denom_rat(typename
                         PT_rat_1::Construct_polynomial()(qcoeffs.begin(),
                                                          qcoeffs.end()));
    Integer denom_numer_int,denom_denom_int;
    typename FT_poly_rat_1::Decompose()(numer_rat,_numer,denom_numer_int);
    typename FT_poly_rat_1::Decompose()(denom_rat,_denom,denom_denom_int);
    _numer *= denom_denom_int;
    _denom *= denom_numer_int;

    _init(_numer,_denom,x_s,x_t,cache);
  }

  void _init(const Polynomial_1& P_, const Polynomial_1& Q_,
             const Algebraic_real_1& x_s, const Algebraic_real_1& x_t,
             const Cache& cache)
  {
    CGAL_precondition(!CGAL::is_zero(Q_));
    //set rational function
    Polynomial_1 P;
    Polynomial_1 Q;
    _canonicalize(P_, Q_, P, Q);
    _f = get_rational_function(P, Q, cache);

    // Compare the x-coordinates and determine the direction.
    Comparison_result   x_res = CGAL::compare(x_s, x_t);
    CGAL_precondition(x_res != EQUAL);

    if (x_res == SMALLER)
      _info = (_info | IS_DIRECTED_RIGHT);

        //Set the source point and check if it lies next to a pole.
    _ps = Algebraic_point_2(_f, x_s);

    //check if source point lies next to a pole.
    if (typename Algebraic_kernel::Sign_at_1()(Q,x_s) != CGAL::ZERO)
    {
      // We have a nomral endpoint.
      //nothing to do ..
    }
    else
    {
      // The y-coodinate is unbounded, but we can set its sign.
      std::pair<CGAL::Sign, CGAL::Sign>  signs = _analyze_near_pole(x_s);
      const CGAL::Sign sign_s =
        ((_info & IS_DIRECTED_RIGHT) != 0) ? signs.second : signs.first;

      if (sign_s == CGAL::NEGATIVE)
        _info = (_info | SRC_AT_Y_MINUS_INFTY);
      else
        _info = (_info | SRC_AT_Y_PLUS_INFTY);
    }

    //Set the target point and check if it lies next to a pole.
    _pt=Algebraic_point_2(_f,x_t);

    //check if target point lies next to a pole.
    if (typename Algebraic_kernel::Sign_at_1()(Q,x_t) != CGAL::ZERO)
    {
      // We have a nomral endpoint.
      //nothing to do ..
    }
    else
    {
      // The y-coodinate is unbounded, but we can set its sign.
      std::pair<CGAL::Sign, CGAL::Sign>  signs = _analyze_near_pole(x_t);
      const CGAL::Sign sign_t =
        ((_info & IS_DIRECTED_RIGHT) != 0) ? signs.first : signs.second;

      if (sign_t == CGAL::NEGATIVE)
        _info = (_info | TRG_AT_Y_MINUS_INFTY);
      else
        _info = (_info | TRG_AT_Y_PLUS_INFTY);
    }

    //Mark that the arc is valid. As it may have poles, we mark it
    //as continuous only if the denominator has no roots.
    _info = ( _info | ( this->_is_continuous() ?
            (IS_CONTINUOUS | IS_VALID) : IS_VALID ) );

  }

  //-----------------------------
  // Accessing the arc properties
  //-----------------------------

  //-----------------------------------------------------------------
  //Get the numerator polynomial of the underlying rational function.
  const Polynomial_1& numerator() const
  {
    return(_f.numer());
  }

  //------------------------------------------------------------------
  //Get the denominator polynomial of the underlying rational function
  const Polynomial_1& denominator() const
  {
    return(_f.denom());
  }

  //---------------------------------------------------------
  //Check if the x-coordinate of the source point is infinite
  Arr_parameter_space source_parameter_space_in_x() const
  {
    return
      ((_info & SRC_AT_X_MINUS_INFTY) != 0) ? ARR_LEFT_BOUNDARY :
      ((_info & SRC_AT_X_PLUS_INFTY) != 0) ? ARR_RIGHT_BOUNDARY :
      ARR_INTERIOR;
  }

  //---------------------------------------------------------
  //Check if the y-coordinate of the source point is infinite
  Arr_parameter_space source_parameter_space_in_y() const
  {
    return
      ((_info & SRC_AT_Y_MINUS_INFTY) != 0) ? ARR_BOTTOM_BOUNDARY :
      ((_info & SRC_AT_Y_PLUS_INFTY) != 0) ? ARR_TOP_BOUNDARY :
      ARR_INTERIOR;
  }

  //---------------------------------------------------------
  //Check if the x-coordinate of the target point is infinite
  Arr_parameter_space target_boundary_in_x() const
  {
    return
      ((_info & TRG_AT_X_MINUS_INFTY) != 0) ? ARR_LEFT_BOUNDARY :
      ((_info & TRG_AT_X_PLUS_INFTY) != 0) ? ARR_RIGHT_BOUNDARY :
      ARR_INTERIOR;
  }

  //---------------------------------------------------------
  //Check if the y-coordinate of the target point is infinite
  Arr_parameter_space target_boundary_in_y() const
  {
    return ((_info & TRG_AT_Y_MINUS_INFTY) != 0) ? ARR_BOTTOM_BOUNDARY :
      ((_info & TRG_AT_Y_PLUS_INFTY) != 0) ? ARR_TOP_BOUNDARY :
      ARR_INTERIOR;
  }

  //--------------------
  //Get the source point
  const Algebraic_point_2& source() const
  {
    CGAL_precondition((_info & IS_VALID) != 0 &&
                      source_parameter_space_in_x() == ARR_INTERIOR &&
                      source_parameter_space_in_y() == ARR_INTERIOR);
    return (_ps);
  }

  //----------------------------------------
  //Get the x-coordinate of the source point
  Algebraic_real_1 source_x() const
  {
    CGAL_precondition((_info & IS_VALID) != 0 &&
                      source_parameter_space_in_x() == ARR_INTERIOR);
    return (_ps.x());
  }

  //----------------------------------------
  //Get the y-coordinate of the source point
  //TODO: should we eliminate the function???
  //Algebraic source_y () const
  //{
  // CGAL_precondition ((_info & IS_VALID) != 0 &&
  //        source_parameter_space_in_y() == ARR_INTERIOR);
  // return (_ps.y());
  //}

  //--------------------
  //Get the target point
  const Algebraic_point_2& target() const
  {
    CGAL_precondition((_info & IS_VALID) != 0 &&
                      target_boundary_in_x() == ARR_INTERIOR &&
                      target_boundary_in_y() == ARR_INTERIOR);
    return (_pt);
  }

  //----------------------------------------
  //Get the x-coordinate of the target point
  Algebraic_real_1 target_x() const
  {
    CGAL_precondition((_info & IS_VALID) != 0 &&
                      target_boundary_in_x() == ARR_INTERIOR);
    return (_pt.x());
  }

  //----------------------------------------
  //Get the y-coordinate of the target point
  //TODO: should we eliminate the function???
  //Algebraic target_y () const
  //{
  // CGAL_precondition ((_info & IS_VALID) != 0 &&
  //        target_boundary_in_y() == ARR_INTERIOR);
  // return (_pt.y());
  //}

  //-------------------------------------------------------
  //Check if the x-coordinate of the left point is infinite
  Arr_parameter_space left_parameter_space_in_x() const
  {
    return ((_info & IS_DIRECTED_RIGHT) != 0) ?
      source_parameter_space_in_x() : target_boundary_in_x();
  }

  //-------------------------------------------------------
  //Check if the y-coordinate of the left point is infinite
  Arr_parameter_space left_parameter_space_in_y() const
  {
    return ((_info & IS_DIRECTED_RIGHT) != 0) ?
      source_parameter_space_in_y() : target_boundary_in_y();
  }
  //--------------------------------------------------------
  //Check if the x-coordinate of the right point is infinite
  Arr_parameter_space right_parameter_space_in_x() const
  {
    return ((_info & IS_DIRECTED_RIGHT) != 0) ?
      target_boundary_in_x() : source_parameter_space_in_x();
  }

  //--------------------------------------------------------
  //Check if the y-coordinate of the right point is infinite
  Arr_parameter_space right_parameter_space_in_y() const
  {
    return ((_info & IS_DIRECTED_RIGHT) != 0) ?
      target_boundary_in_y() : source_parameter_space_in_y();
  }

  //------------------------------------
  //Get the x_value of the left endpoint
  const Algebraic_real_1& left_x() const
  {
    CGAL_precondition(left_parameter_space_in_x() == ARR_INTERIOR);
    return ((_info & IS_DIRECTED_RIGHT) ? _ps.x() : _pt.x());
  }

  //------------------------------------
  //Get the x_value of the right endpoint
  const Algebraic_real_1& right_x() const
  {
    CGAL_precondition((_info & IS_VALID) != 0 &&
        right_parameter_space_in_x() == ARR_INTERIOR );
    return ((_info & IS_DIRECTED_RIGHT) ? _pt.x() : _ps.x());
  }
  //---------------------
  //Get the left endpoint
  const Algebraic_point_2& left() const
  {
    CGAL_precondition(left_parameter_space_in_x() == ARR_INTERIOR &&
        left_parameter_space_in_y() == ARR_INTERIOR);
    return ((_info & IS_DIRECTED_RIGHT) ? _ps : _pt);
  }

  //----------------------
  //Get the right endpoint
  const Algebraic_point_2& right() const
  {
    CGAL_precondition(right_parameter_space_in_x() == ARR_INTERIOR &&
        right_parameter_space_in_y() == ARR_INTERIOR);
    return ((_info & IS_DIRECTED_RIGHT) ? _pt : _ps);
  }

  //-------------------------
  //Check if the arc is valid
  bool is_valid() const
  {
    return ((_info & IS_VALID) != 0);
  }

  //------------------------------
  //Check if the arc is continuous
  bool is_continuous() const
  {
    return ((_info & IS_CONTINUOUS) != 0);
  }

  //----------------------------------
  //Check if the arc is directed right
  bool is_directed_right() const
  {
    return ((_info & IS_DIRECTED_RIGHT) != 0);
  }

  //--------------
  //name Modifiers
  //---------------

  //--------------------------------
  //Mark the arc as being continuous
  void set_continuous()
  {
    _info = (_info | IS_CONTINUOUS);
  }

  //-----------------------------
  //Mark the arc as being invalid
  void set_invalid()
  {
    _info = (_info & ~IS_VALID);
  }

//   bool is_intersecting(Self& arc)
//   {
//     Arr_parameter_space left_parameter_space = ( (left_parameter_space_in_x()!= ARR_INTERIOR) &&(arc.left_parameter_space_in_x() != ARR_INTERIOR) ) ?
//       ARR_LEFT_BOUNDARY :
//       ARR_INTERIOR;
//     Arr_parameter_space right_parameter_space = ( (right_parameter_space_in_x()!= ARR_INTERIOR) &&(arc.right_parameter_space_in_x()!= ARR_INTERIOR) ) ?
//       ARR_RIGHT_BOUNDARY :
//       ARR_INTERIOR;
//     Algebraic_real_1 left,right;
//     if (left_parameter_space == ARR_INTERIOR)
//       left =
//         (left_parameter_space_in_x()!= ARR_INTERIOR)     ? arc.left_x():
//         (arc.left_parameter_space_in_x()!= ARR_INTERIOR) ? left_x() :
//         (arc.left_x() < left_x())                 ? left_x() :
//         arc.left_x();
//     if (right_parameter_space == ARR_INTERIOR)
//       right = (right_parameter_space_in_x()!= ARR_INTERIOR)  ? arc.right_x():
//         (arc.right_parameter_space_in_x()!= ARR_INTERIOR) ? right_x() :
//         (arc.right_x() < right_x())     ? right_x() :
//         arc.right_x();
//     if (left > right)
//       return false;

//     //check if the base functions are equal
//     if (_has_same_base (arc))
//       return true;
//     Rational_function_pair rat_pair = get_rational_pair(get_rational_function(this->_numer,this->_denom),
//         get_rational_function(arc._numer,arc._denom));
//     return rat_pair.is_intersecting_in_range(left_parameter_space,left,right_parameter_space,right);
//   }

  //--------------------------------------------------------
  //Split the arc into two at a given pole.
  //The function returns the sub-arc to the left of the pole
  //and sets (*this) to be the right sub-arc.
  //param x0 The x-coordinate of the pole.
  //precondition x0 lies in the interior of the arc.
  //return The sub-arc to the left of the pole.

  Self split_at_pole(const Algebraic_real_1& x0)
  {
    // Analyze the behaviour of the function near the given pole.
    const std::pair<CGAL::Sign, CGAL::Sign>  signs = _analyze_near_pole(x0);
    const CGAL::Sign    sign_left = signs.first;
    const CGAL::Sign    sign_right = signs.second;

    // Create a fictitious point that represents the x-coordinate of the pole.
    Algebraic_point_2 p0(_f, x0);

    // Make a copy of the current arc.
    Self       c1 = *this;

    // Split the arc, such that c1 lies to the left of the pole and (*this)
    // to its right.
    if ((_info & IS_DIRECTED_RIGHT) != 0)
    {
      c1._pt = p0;
      c1._info = (c1._info & ~TRG_INFO_BITS);
      if (sign_left == CGAL::NEGATIVE)
        c1._info = (c1._info | TRG_AT_Y_MINUS_INFTY);
      else
        c1._info = (c1._info | TRG_AT_Y_PLUS_INFTY);

      this->_ps = p0;
      this->_info = (this->_info & ~SRC_INFO_BITS);
      if (sign_right == CGAL::NEGATIVE)
        this->_info = (this->_info | SRC_AT_Y_MINUS_INFTY);
      else
        this->_info = (this->_info | SRC_AT_Y_PLUS_INFTY);
    }
    else
    {
      c1._ps = p0;
      c1._info = (c1._info & ~SRC_INFO_BITS);
      if (sign_left == CGAL::NEGATIVE)
        c1._info = (c1._info | SRC_AT_Y_MINUS_INFTY);
      else
        c1._info = (c1._info | SRC_AT_Y_PLUS_INFTY);

      this->_pt = p0;
      this->_info = (this->_info & ~TRG_INFO_BITS);
      if (sign_right == CGAL::NEGATIVE)
        this->_info = (this->_info | TRG_AT_Y_MINUS_INFTY);
      else
        this->_info = (this->_info | TRG_AT_Y_PLUS_INFTY);
    }

    // Mark the sub-arc c1 as continuous.
    c1._info = (c1._info | IS_CONTINUOUS);

    return (c1);
  }

  //---------------
  //name Predicates
  //---------------


  //---------------------------------------------------------------------------
  //Get the relative position of the point with respect to the rational arc.
  //param p The query point.
  //precondition: p is in the x-range of the arc.
  //    both p's supporting curve and the rational arc are continous
  //return SMALLER if the point is below the arc;
  //       LARGER if the point is above the arc;
  //       EQUAL if p lies on the arc.
  Comparison_result point_position(const Algebraic_point_2& p,
                                   const Cache& cache) const
  {
    // Make sure that p is in the x-range of the arc and check whether it
    // has the same x-coordinate as one of the endpoints.
    CGAL_precondition(is_continuous());
    CGAL_precondition(_is_in_true_x_range(p.x()));
    if (p.rational_function() == _f)
      return EQUAL;
    Rational_function_pair rat_pair(get_rational_pair(p.rational_function(), _f,
                                                      cache));
    return rat_pair.compare_f_g_at(p.x());
  }
  //---------------------------------------------------------------------------
  //Compare the x-coordinate of a vertical asymptote of the arc
  //(one of its ends) and the given point.
  Comparison_result compare_end(Arr_curve_end ce,
                                const Algebraic_point_2& p) const
  {
    Algebraic_real_1 x0;

    if (ce == ARR_MIN_END) {
      CGAL_assertion(left_parameter_space_in_x() == ARR_INTERIOR &&
                     left_parameter_space_in_y() != ARR_INTERIOR);
      if ((_info & IS_DIRECTED_RIGHT) != 0)
        x0 = _ps.x();
      else
        x0 = _pt.x();
    }
    else
    {
      CGAL_assertion(right_parameter_space_in_x() == ARR_INTERIOR &&
                     right_parameter_space_in_y() != ARR_INTERIOR);
      if ((_info & IS_DIRECTED_RIGHT) != 0)
        x0 = _pt.x();
      else
        x0 = _ps.x();
    }

    // Compare the x-coordinates.
    const Comparison_result  res = CGAL::compare(x0, p.x());

    if (res != EQUAL)
      return (res);

    return ((ce == ARR_MIN_END) ?  LARGER : SMALLER);
  }

  //------------------------------------------------------------------
  //Compare the x-coordinate of a vertical asymptotes of the two arcs.
  //approaching from the same direction
  Comparison_result compare_near_end(const Self& arc, Arr_curve_end ce,
                                     const Cache& cache) const
  {
    CGAL_precondition_code(
      if (ce == ARR_MIN_END)
      {
        CGAL_precondition(this->left_parameter_space_in_y() != ARR_INTERIOR);
        CGAL_precondition(arc.left_parameter_space_in_y() != ARR_INTERIOR);
        CGAL_precondition(this->left_parameter_space_in_y() ==
                          arc.left_parameter_space_in_y());
      }
      else // (ce == ARR_MAX_END)
      {
        CGAL_precondition(this->right_parameter_space_in_y() != ARR_INTERIOR);
        CGAL_precondition(arc.right_parameter_space_in_y() != ARR_INTERIOR);
        CGAL_precondition(this->right_parameter_space_in_y() ==
                          arc.right_parameter_space_in_y());
      }
      );

    // Get the x-coordinates of the vertical asymptote.
    Algebraic_real_1 x((ce == ARR_MIN_END)? this->left_x() : this->right_x());
    CGAL_assertion(CGAL::compare(x, (ce == ARR_MIN_END)?
                                 arc.left_x () : arc.right_x()) == EQUAL);

    //both arcs have vertical asymptotes and come from the same side of the
    //x-axis compare value of functions close to the root on the correct side
    if (_has_same_base(arc))
      return CGAL::EQUAL;
    Rational_function_pair rat_pair = get_rational_pair(_f,arc._f,cache);

    CGAL::Comparison_result comp_f_g_y =
      rat_pair.compare_f_g_at(x,ce == ARR_MAX_END ?
                              CGAL::NEGATIVE : CGAL::POSITIVE);
    if (ce == ARR_MAX_END)
    {
      return (right_parameter_space_in_y() == ARR_BOTTOM_BOUNDARY) ?
        comp_f_g_y : -comp_f_g_y;
    }
    else
    {
      return (left_parameter_space_in_y() == ARR_BOTTOM_BOUNDARY) ?
        -comp_f_g_y : comp_f_g_y;
    }
  }

  //------------------------------------------------------------------
  //Compare the x-coordinate of a vertical asymptotes of the two arcs.
  Comparison_result compare_ends(Arr_curve_end ind1,const Self& arc,
                                 Arr_curve_end ind2, const Cache& cache) const
  {
    // Get the x-coordinates of the first vertical asymptote.
    Algebraic_real_1             x1;

    if (ind1 == ARR_MIN_END)
    {
      CGAL_assertion(left_parameter_space_in_x() == ARR_INTERIOR &&
                     left_parameter_space_in_y() != ARR_INTERIOR);
      if ((_info & IS_DIRECTED_RIGHT) != 0)
        x1 = _ps.x();
      else
        x1 = _pt.x();
    }
    else
    {
      CGAL_assertion(right_parameter_space_in_x() == ARR_INTERIOR &&
                     right_parameter_space_in_y() != ARR_INTERIOR);
      if ((_info & IS_DIRECTED_RIGHT) != 0)
        x1 = _pt.x();
      else
        x1 = _ps.x();
    }

    // Get the x-coordinates of the second vertical asymptote.
    Algebraic_real_1                x2;

    if (ind2 == ARR_MIN_END)
    {
      CGAL_assertion(arc.left_parameter_space_in_x() == ARR_INTERIOR &&
                     arc.left_parameter_space_in_y() != ARR_INTERIOR);
      if ((arc._info & IS_DIRECTED_RIGHT) != 0)
        x2 = arc._ps.x();
      else
        x2 = arc._pt.x();
    }
    else
    {
      CGAL_assertion(arc.right_parameter_space_in_x() == ARR_INTERIOR &&
                     arc.right_parameter_space_in_y() != ARR_INTERIOR);
      if ((arc._info & IS_DIRECTED_RIGHT) != 0)
        x2 = arc._pt.x();
      else
        x2 = arc._ps.x();
    }

    // Compare the x-coordinates. In case they are not equal we are done.
    const Comparison_result  res = CGAL::compare(x1, x2);

    if (res != EQUAL)
      return (res);

    // If the x-coordinates of the asymptote are equal, but one arc is
    // defined to the left of the vertical asymptote and the other to its
    // right, we can easily determine the comparison result.
    if (ind1 == ARR_MAX_END && ind2 == ARR_MIN_END)
      return (SMALLER);
    else if (ind1 == ARR_MIN_END && ind2 == ARR_MAX_END)
      return (LARGER);

    //both arcs have vertical asymptotes and come from the same side of the x-axis
    //compare value of functions close to the root on the correct side
    if (_has_same_base(arc))
      return CGAL::EQUAL;
    Rational_function_pair rat_pair = get_rational_pair(_f,arc._f,cache);

    CGAL::Comparison_result comp_f_g_y =
      rat_pair.compare_f_g_at(x1, ind1 == ARR_MAX_END ?
                              CGAL::NEGATIVE : CGAL::POSITIVE);
    if( ind1 == ARR_MAX_END)
    {
      CGAL_postcondition(ind2 == ARR_MAX_END);
      CGAL_postcondition(right_parameter_space_in_y() ==
                         arc.right_parameter_space_in_y());
      return (right_parameter_space_in_y() == ARR_BOTTOM_BOUNDARY) ?
        comp_f_g_y : -comp_f_g_y;
    }
    else
    {
      CGAL_postcondition(ind2 == ARR_MIN_END);
      CGAL_postcondition(left_parameter_space_in_y() ==
                         arc.left_parameter_space_in_y());
      return (left_parameter_space_in_y() == ARR_BOTTOM_BOUNDARY) ?
        -comp_f_g_y : comp_f_g_y;
    }
  }

  //------------------------------------------------------------------
  // Compare the slopes of the arc with another given arc
  // at their given intersection point.
  // param cv The given arc.
  // param p The intersection point.
  // param mult Output: The mutiplicity of the intersection point.
  // return SMALLER if (*this) slope is less than cv's;
  //      EQUAL if the two slopes are equal;
  //         LARGER if (*this) slope is greater than cv's.
  //Comparison_result compare_slopes (const Self& arc,const Algebraic_point_2& p,unsigned int& mult) const
  //
  // deleted!!!

  //------------------------------------------------------------------
  // Compare the two arcs at a given intersection point
  // param arc The given arc
  // param p  The intersection point
  // param to_left to check to the left or to the right of intersection point
  // precondition: Both arcs intersect at p
  // return SMALLER if (*this) lies below the other arc;
  //   EQUAL if the two supporting functions are equal;
  //   LARGER if (*this) lies above the other arc.

  Comparison_result compare_at_intersection(const Self& arc,
                                            const Algebraic_point_2& p,
                                            bool to_left, const Cache& cache) const
  {
    CGAL_precondition(this->point_position(p,cache) == CGAL::EQUAL &&
                      arc.point_position(p,cache)   == CGAL::EQUAL);

    //check if the base functions are equal
    if (_has_same_base(arc))
      return CGAL::EQUAL;

    Rational_function_pair rat_pair = get_rational_pair(this->_f, arc._f,cache);
    return rat_pair.compare_f_g_at(p.x(),
                                   to_left ? CGAL::SMALLER : CGAL::LARGER);
  }

  //------------------------------------------------------------------
  // Compare the two arcs at x = -oo.
  // param arc The given arc
  // precondition: Both arcs have a left end which is unbounded in x.
  // return SMALLER if (*this) lies below the other arc;
  //   EQUAL if the two supporting functions are equal;
  //   LARGER if (*this) lies above the other arc.

  Comparison_result compare_at_minus_infinity(const Self& arc,
                                              const Cache& cache) const
  {
    CGAL_precondition(left_parameter_space_in_x() == ARR_LEFT_BOUNDARY &&
                      arc.left_parameter_space_in_x() == ARR_LEFT_BOUNDARY);

    //check if the base functions are equal
    if (_has_same_base(arc))
      return CGAL::EQUAL;

    Rational_function_pair rat_pair =
      get_rational_pair(this->_f, arc._f, cache);
    return rat_pair.compare_f_g_at(ARR_LEFT_BOUNDARY);
  }

  //------------------------------------------------------------------
  //Compare the two arcs at x = +oo.
  //param arc The given arc.
  //pre Both arcs are have a right end which is unbounded in x.
  //return SMALLER if (*this) lies below the other arc;
  //       EQUAL if the two supporting functions are equal;
  //       LARGER if (*this) lies above the other arc.

  Comparison_result compare_at_plus_infinity(const Self& arc,
                                             const Cache& cache) const
  {
    CGAL_precondition(right_parameter_space_in_x() == ARR_RIGHT_BOUNDARY &&
                      arc.right_parameter_space_in_x() == ARR_RIGHT_BOUNDARY);

    //check if the base functions are equal
    if (_has_same_base(arc))
      return CGAL::EQUAL;

    Rational_function_pair rat_pair = get_rational_pair(this->_f, arc._f, cache);
    return rat_pair.compare_f_g_at(ARR_RIGHT_BOUNDARY);
  }

  //----------------------------------------------------------
  // Check whether the two arcs are equal (have the same graph).
  //param arc The compared arc.
  //return true if the two arcs have the same graph; false otherwise.
  bool equals(const Self& arc) const
  {
    // The two arc must have the same base rational function.
    CGAL_precondition(is_valid());
    CGAL_precondition(arc.is_valid());
    if (! _has_same_base(arc))
      return false;

    // Check that the arc left endpoints are the same.
    Arr_parameter_space inf1 = left_parameter_space_in_x();
    Arr_parameter_space inf2 = arc.left_parameter_space_in_x();

    if (inf1 != inf2)
      return false;

    if (inf1 == ARR_INTERIOR)
    {
      inf1 = left_parameter_space_in_y();
      inf2 = arc.left_parameter_space_in_y();

      if (inf1 != inf2)
        return false;
    }

    if (inf1 == ARR_INTERIOR &&
        CGAL::compare(left().x(), arc.left().x()) != EQUAL)
    {
      return false;
    }

    // Check that the arc right endpoints are the same.
    inf1 = right_parameter_space_in_x();
    inf2 = arc.right_parameter_space_in_x();

    if (inf1 != inf2)
      return false;

    if (inf1 == ARR_INTERIOR)
    {
        inf1 = right_parameter_space_in_y();
        inf2 = arc.right_parameter_space_in_y();

        if (inf1 != inf2)
          return false;
    }

    if (inf1 == ARR_INTERIOR &&
        CGAL::compare(right().x(), arc.right().x()) != EQUAL)
      return false;

    // If we reached here, the two arc are equal:
    return true;
  }

  //----------------------------------------------------------
  //Check whether it is possible to merge the arc with the given arc.
  //param arc The query arc.
  //return true if it is possible to merge the two arcs;
  //   false otherwise.
  bool can_merge_with(const Self& arc) const
  {
    // In order to merge the two arcs, they should have the same base
    // rational function.
    CGAL_precondition(is_valid());
    CGAL_precondition(arc.is_valid());
    if (! _has_same_base(arc))
      return false;

    // Check if the left endpoint of one curve is the right endpoint of
    // the other.

    return ((right_parameter_space_in_x() == ARR_INTERIOR &&
             right_parameter_space_in_y() == ARR_INTERIOR &&
             arc.left_parameter_space_in_x() == ARR_INTERIOR &&
             arc.left_parameter_space_in_y() == ARR_INTERIOR &&
             //CGAL::equal(right_x() ,arc.left_x() )) ||
             (right_x() == arc.left_x() )) ||
            (left_parameter_space_in_x() == ARR_INTERIOR &&
             left_parameter_space_in_y() == ARR_INTERIOR &&
             arc.right_parameter_space_in_x() == ARR_INTERIOR &&
             arc.right_parameter_space_in_y() == ARR_INTERIOR &&
             //CGAL::equal(left_x(), arc.right_x())));
             (left_x() == arc.right_x())));
  }
  //------------------------------------
  //Constructions of points and curves
  //------------------------------------

  //------------------------------------
  // Flip the arc (swap its source and target).
  // return The flipped arc.
  Self flip() const
  {
    CGAL_precondition(is_valid());

    // Create the flipped arc.
    Self   arc;

    arc._f = _f;
    arc._ps = _pt;
    arc._pt = _ps;

    // Manipulate the information bits.
    int    src_info = (_info & SRC_INFO_BITS);
    int    trg_info = (_info & TRG_INFO_BITS);
    arc._info = (src_info << 4) | (trg_info >> 4) | IS_VALID;

    if ((_info & IS_DIRECTED_RIGHT) == 0)
      arc._info = (arc._info | IS_DIRECTED_RIGHT);

    if ((_info & IS_CONTINUOUS) != 0)
      arc._info = (arc._info | IS_CONTINUOUS);

    return (arc);
  }

  //------------------------
  // Print the rational arc.
  std::ostream& print(std::ostream& os) const
  {
    // Print y as a rational function of x.
    os << "y = (";
    Base_rational_arc_ds_1::print_polynomial(os, this->numerator(), 'x');
    os << ") / (";
    Base_rational_arc_ds_1::print_polynomial(os, this->denominator(), 'x');
    os << ") on ";

    // Print the definition range.
    Arr_parameter_space      inf_x = source_parameter_space_in_x();
    if (inf_x == ARR_LEFT_BOUNDARY)
      os << "(-oo";
    else if (inf_x == ARR_RIGHT_BOUNDARY)
      os << "(+oo";
    else if (source_parameter_space_in_y() != ARR_INTERIOR)
      os << '(' << CGAL::to_double(source_x());
    else
      os << '[' << CGAL::to_double(source().x());

    os << ", ";

    inf_x = target_boundary_in_x();
    if (inf_x == ARR_LEFT_BOUNDARY)
      os << "-oo)";
    else if (inf_x == ARR_RIGHT_BOUNDARY)
      os << "+oo)";
    else if (target_boundary_in_y() != ARR_INTERIOR)
      os << CGAL::to_double(target_x()) << ')';
    else
      os << CGAL::to_double(target().x()) << ']';

    return (os);
  }

protected:


  //-------------------------------
  //Auxiliary (protected) functions.
  //-------------------------------

  //--------------------------------------------------------------------------
  // Cannonicalize numerator and denominator such that:
  //  There are no common devisor
  //  If negative sign exists, it is in the numerator
  void _canonicalize(const Polynomial_1& P,const Polynomial_1& Q,
                     Polynomial_1& P_new  ,Polynomial_1& Q_new)
  {
    Polynomial_1 gcd = CGAL::gcd(P,Q);
    if (gcd != 1)
    {
      P_new=CGAL::integral_division(P,gcd);
      Q_new=CGAL::integral_division(Q,gcd);
    }
    else
    {
      P_new=P;
      Q_new=Q;
    }

    if (typename AT_poly::Unit_part()(Q_new) == -1) //leading coefficient sign
    {
      P_new*=-1;
      Q_new*=-1;
    }
    return;
  }

  //--------------------------------------------------------------------------
  // Check if the given x-value is in the x-range of the arc.
  // param x The x-value.
  // param eq_src Output: Is this value equal to the x-coordinate of the
  //                       source point.
  // param eq_trg Output: Is this value equal to the x-coordinate of the
  //                       target point.

  bool _is_in_x_range(const Algebraic_real_1& x, bool& eq_src,
                      bool& eq_trg) const
  {
    Comparison_result  res1;

    eq_src = eq_trg = false;
    if ((_info & IS_DIRECTED_RIGHT) != 0)
    {
      // Compare to the left endpoint (the source in this case).
      if ((_info & SRC_AT_X_MINUS_INFTY) != 0)
      {
        res1 = LARGER;
      }
      else
      {
        res1 = CGAL::compare(x, _ps.x());

        if (res1 == SMALLER)
          return false;

        if (res1 == EQUAL)
        {
          eq_src = true;
          return true;
        }
      }

      // Compare to the right endpoint (the target in this case).
      if ((_info & TRG_AT_X_PLUS_INFTY) != 0)
        return true;

      const Comparison_result  res2 = CGAL::compare(x, _pt.x());

      if (res2 == LARGER)
        return false;

      if (res2 == EQUAL)
        eq_trg = true;

      return true;
    }

    // Compare to the left endpoint (the target in this case).
    if ((_info & TRG_AT_X_MINUS_INFTY) != 0)
    {
      res1 = LARGER;
    }
    else
    {
      res1 = CGAL::compare(x, _pt.x());

      if (res1 == SMALLER)
        return false;

      if (res1 == EQUAL)
      {
        eq_trg = true;
        return true;
      }
    }

    // Compare to the right endpoint (the source in this case).
    if ((_info & SRC_AT_X_PLUS_INFTY) != 0)
      return true;

    const Comparison_result  res2 = CGAL::compare(x, _ps.x());

    if (res2 == LARGER)
      return false;

    if (res2 == EQUAL)
      eq_src = true;

    return true;
  }

  //----------------------------------------------------------
  // Check if the given x-value is in the x-range of the arc,
  // excluding its open ends.
  bool _is_in_true_x_range(const Algebraic_real_1& x) const
  {
    bool          eq_src, eq_trg;
    const bool    is_in_x_range_closure = _is_in_x_range(x, eq_src, eq_trg);

    if (! is_in_x_range_closure)
      return false;

    // Check if we have a vertical asymptote at the source point.
    if (eq_src && (_info & (SRC_AT_Y_MINUS_INFTY | SRC_AT_Y_PLUS_INFTY)) != 0)
      return false;

    // Check if we have a vertical asymptote at the target point.
    if (eq_trg && (_info & (TRG_AT_Y_MINUS_INFTY | TRG_AT_Y_PLUS_INFTY)) != 0)
      return false;

    // If we reached here, the value is in the true x-range of the arc.
    return true;
  }

  //------------------------------------------------------------------------
  // Check if the underlying rational function is the same in the given arc.
  //  param arc The given arc.
  //   return true if arc's underlying rational function is the same
  //       as of *this;
  //   false otherwise.

  bool _has_same_base(const Self& arc) const
  {
    return (this->_f == arc._f);
  }
  //bool operator == (const Self& arc) const
  //{
  //  if (this == &arc)
  //    return true;
  //  if ((_info ==arc._info) &&(_has_same_base(arc) ))
  //    {
  //      bool same_source(true);
  //      bool same_target(true);
  //      if (  (this->source_parameter_space_in_x () == ARR_INTERIOR) &&
  //          (this->source_parameter_space_in_y () == ARR_INTERIOR) )
  //        same_source = (this->source() == arc.source());
  //      if (  (this->target_boundary_in_x () == ARR_INTERIOR) &&
  //          (this->target_boundary_in_y () == ARR_INTERIOR) )
  //        same_target = (this->target() == arc.target());
  //      return (same_source && same_target);
  //    }
  //  return false;
  //}

  //---------------------------------------------------------------------
  //Compute infinity type of the rational function P(x)/Q(x) at x = -oo.
  // param y Output: The value of the horizontal asymptote (if exists).
  // return The infinity type for the y-coordinate at x = -oo.

  Arr_parameter_space _analyze_at_minus_infinity(const Polynomial_1& P,
                                                 const Polynomial_1& Q,
                                                 Algebraic_real_1& y) const
  {
    // Get the degree of the polynomials.
    const int    deg_p(CGAL::degree(P));
    const int    deg_q(CGAL::degree(Q));

    if (deg_p <= 0 || deg_p <= deg_q)
      {
        // We have a zero polynomial or a zero asymptote.
        y = 0;
        return (ARR_INTERIOR);
      }
    // Get the leading coefficients.
    Integer p_lead(CGAL::leading_coefficient(P));
    Integer q_lead(CGAL::leading_coefficient(Q));

    if (deg_p == deg_q)
      {
        // We have a horizontal asymptote.
        y = Algebraic_real_1(Rational(p_lead) / Rational(q_lead));
        return (ARR_INTERIOR);
      }

    // We have a tendency to infinity.
    const int    def_diff = deg_p - deg_q;

    return (CGAL::sign(p_lead) == CGAL::sign (q_lead)) ?
      ((def_diff % 2 == 0) ? ARR_TOP_BOUNDARY : ARR_BOTTOM_BOUNDARY) :
      ((def_diff % 2 == 0) ? ARR_BOTTOM_BOUNDARY : ARR_TOP_BOUNDARY);
  }

  //---------------------------------------------------------------------
  //Compute infinity type of the rational function P(x)/Q(x) at x = +oo.
  // param y Output: The value of the horizontal asymptote (if exists).
  // return The infinity type for the y-coordinate at x = +oo.

  Arr_parameter_space _analyze_at_plus_infinity(const Polynomial_1& P,
                                                const Polynomial_1& Q,
                                                Algebraic_real_1& y) const
  {
    // Get the degree of the polynomials.
    const int    deg_p(CGAL::degree(P));
    const int    deg_q(CGAL::degree(Q));

    if (deg_p <= 0 || deg_p <= deg_q)
      {
        // We have a zero polynomial or a zero asymptote.
        y = 0;
        return (ARR_INTERIOR);
      }

    // Get the leading coefficients.
    Integer p_lead(CGAL::leading_coefficient(P));
    Integer q_lead(CGAL::leading_coefficient(Q));

    if (deg_p == deg_q)
      {
        // We have a horizontal asymptote.
        y = Algebraic_real_1(Rational(p_lead) / Rational(q_lead));
        return (ARR_INTERIOR);
      }

    // We have a tendency to infinity.
    return (CGAL::sign(p_lead) == CGAL::sign(q_lead)) ?
      ARR_TOP_BOUNDARY : ARR_BOTTOM_BOUNDARY;
  }

  //---------------------------------------------------------------------
  //Compute all zeros of the denominator polynomial that lie within the
  // x-range of the arc

  template <typename OutputIterator>
  OutputIterator _denominator_roots(OutputIterator oi, bool& root_at_ps,
                                    bool& root_at_pt) const
  {
    root_at_ps = root_at_pt = false;

    if (CGAL::degree(this->denominator()) == 0)
      return (oi);

    // Compute the roots of the denominator polynomial.
    std::list<Algebraic_real_1>                    q_roots;
    bool                                           eq_src, eq_trg;
    typename std::list<Algebraic_real_1>::const_iterator  x_iter;

    //solve for roots without caring for multiplicity
    //hence the usage of the bool var

    std::copy(_f.poles().begin(),_f.poles().end(),std::back_inserter (q_roots));

    // Go over the roots and check whether they lie in the x-range of the arc.
    for (x_iter = q_roots.begin(); x_iter != q_roots.end(); ++x_iter)
    {
      if (_is_in_x_range (*x_iter, eq_src, eq_trg))
      {
        if (eq_src)
        {
          root_at_ps = true;
        }
        else if (eq_trg)
        {
          root_at_pt = true;
        }
        else
        {
          // The root lies in the interior of the arc.
          *oi++ = *x_iter;
        }
      }
    }

    return (oi);
  }

  //---------------------------------------------------------------------
  // Check whether the arc is continuous.
  bool _is_continuous()
  {
    // Compute the roots of the denominator polynomial, and make sure
    // there are none in the range of definition.
    std::list<Algebraic_real_1>          q_roots;
    bool                                 root_at_ps, root_at_pt;

    this->_denominator_roots(std::back_inserter(q_roots),
        root_at_ps, root_at_pt);

    return (q_roots.empty());
  }

  //---------------------------------------------------------------------
  //Determine the signs of the rational functions infinitisimally to the left
  //   and to the right of the given pole.
  //   param x0 The x-coordinate of the pole.
  //   pre x0 lies in the interior of the arc.
  //   return The signs to the left and to the right of x0.
  std::pair<CGAL::Sign, CGAL::Sign>
  _analyze_near_pole(const Algebraic_real_1& x0) const
  {
    return std::make_pair( _f.sign_at(x0,CGAL::NEGATIVE),
        _f.sign_at(x0,CGAL::POSITIVE));
  }

  //---------------------------------------------------------------------
  // Print a polynomial nicely.

  //std::ostream& _print_polynomial (std::ostream& os,
  //    const Polynomial_1& poly,
  //    char var) const
  //{
  //  // Get the degree.
  //  const int    deg = CGAL::degree(poly);
  //
  //  Integer     coeff;
  //  CGAL::Sign  sgn;
  //  int         k;

  //  if (deg < 0)
  //    {
  //      os << '0';
  //      return (os);
  //    }

  //  for (k = deg; k >= 0; k--)
  //    {
  //      //coeff = pt::Get_coefficient()(poly, k);
  //      coeff = CGAL::get_coefficient(poly, k);
  //
  //      if (k == deg)
  //        os << coeff;
  //      else if ((sgn = CGAL::sign (coeff)) == POSITIVE)
  //        os << " + " << coeff;
  //      else if (sgn == NEGATIVE)
  //        os << " - " << -coeff;
  //      else
  //        continue;
  //
  //      if (k > 1)
  //        os << '*' << var << '^' << k;
  //      else if (k == 1)
  //        os << '*' << var;
  //    }

  //  return (os);
  //}

};

//-------------------------------
//! Exporter for rational arcs.
template <typename Algebraic_kernel_>
std::ostream& operator<<(std::ostream& os,
                         const Base_rational_arc_d_1<Algebraic_kernel_> & arc)
{
  return (arc.print(os));
}

/*! \class Continuous_rational_arc_d_1
 * Representation of a continuous portion of a rational function.
 */
template <typename Algebraic_kernel_>
class Continuous_rational_arc_d_1:
    public Base_rational_arc_d_1<Algebraic_kernel_>
{
public:
  bool is_left_to_right() const
  { return (this->_info & Base::IS_DIRECTED_RIGHT) != 0; }

public:
  typedef Algebraic_kernel_                             Algebraic_kernel;

  typedef Continuous_rational_arc_d_1<Algebraic_kernel> Self;
  typedef Base_rational_arc_d_1<Algebraic_kernel>       Base;

  typedef typename Base::Integer                        Integer;
  typedef typename Base::Rational                       Rational;
  typedef typename Base::Algebraic_real_1               Algebraic_real_1;
  typedef typename Base::Algebraic_point_2              Algebraic_point_2;
  typedef typename Base::Polynomial_1                   Polynomial_1;
  typedef typename Base::Multiplicity                   Multiplicity;
  typedef typename Base::Rational_function              Rational_function;
  typedef typename Base::Rational_function_pair         Rational_function_pair;

  typedef typename Base::Rat_vector                     Rat_vector;
  typedef typename Base::Algebraic_vector               Algebraic_vector;
  typedef typename Base::Multiplicity_vector            Multiplicity_vector;

  typedef typename Base::Cache                          Cache;

  typedef std::pair<Algebraic_point_2, Multiplicity>    Intersection_point;
  //typedef std::pair<Algebraic_point_2, unsigned int>  Intersection_point;


  /// \name Constrcution methods.
  //@{

  /*!
   * Default constructor.
   */
  Continuous_rational_arc_d_1() :
    Base()
  {}

  /*!
   * Constrcutor from a base arc.
   */
  Continuous_rational_arc_d_1(const Base& arc) :
    Base(arc)
  {
    CGAL_precondition(arc.is_continuous());
  }

  /*!
   * Constructor of a whole polynomial curve.
   * \param pcoeffs The rational coefficients of the polynomial p(x).
   */
  Continuous_rational_arc_d_1(const Polynomial_1& P, const Cache& cache) :
    Base(P, cache)
  {}

  Continuous_rational_arc_d_1(const Rat_vector& pcoeffs, const Cache& cache) :
    Base(pcoeffs, cache)
  {}

  /*!
   * Constructor of a polynomial ray, defined by y = p(x), for x_s <= x if the
   * ray is directed to the right, or for x_s >= x if it is directed to the
   * left.
   * \param pcoeffs The rational coefficients of the polynomial p(x).
   * \param x_s The x-coordinate of the source point.
   * \param dir_right Is the ray directed to the right (to +oo)
   *                  or to the left (to -oo).
   */
  Continuous_rational_arc_d_1(const Polynomial_1& P,
                              const Algebraic_real_1& x_s,
                              bool dir_right, const Cache& cache) :
    Base(P, x_s, dir_right,cache)
  {}

  Continuous_rational_arc_d_1(const Rat_vector& pcoeffs,
                              const Algebraic_real_1& x_s,
                              bool dir_right, const Cache& cache) :
    Base(pcoeffs, x_s, dir_right,cache)
  {}


  /*!
   * Constructor of a polynomial arc, defined by y = p(x), x_min <= x <= x_max.
   * \param pcoeffs The rational coefficients of the polynomial p(x).
   * \param x_s The x-coordinate of the source point.
   * \param x_t The x-coordinate of the target point.
   * \pre The two x-coordinates must not be equal.
   */
  Continuous_rational_arc_d_1(const Polynomial_1& P,
                              const Algebraic_real_1& x_s,
                              const Algebraic_real_1& x_t, const Cache& cache) :
    Base(P, x_s, x_t,cache)
  {}

  Continuous_rational_arc_d_1(const Rat_vector& pcoeffs,
                              const Algebraic_real_1& x_s,
                              const Algebraic_real_1& x_t, const Cache& cache) :
    Base(pcoeffs, x_s, x_t,cache)
  {}

  /*!
   * Constructor of a polynomial function, defined by y = p(x)/q(x) for any x.
   * \param pcoeffs The rational coefficients of the polynomial p(x).
   * \param qcoeffs The rational coefficients of the polynomial q(x).
   * \pre The denominator polynomial q(x) does not have any roots.
   */
  Continuous_rational_arc_d_1(const Polynomial_1& P, const Polynomial_1& Q,
                              const Cache& cache) :
    Base(P, Q,cache)
  {
    if (!this->_is_continuous())
    {
      // Invalid arc, as it is not continuous.
      this->set_invalid();
    }
  }

  Continuous_rational_arc_d_1(const Rat_vector& pcoeffs,
                              const Rat_vector& qcoeffs, const Cache& cache) :
    Base(pcoeffs, qcoeffs,cache)
  {
    if (!this->_is_continuous())
    {
      // Invalid arc, as it is not continuous.
      this->set_invalid();
    }
  }

  /*!
   * Constructor of a ray of a rational function, defined by y = p(x)/q(x),
   * for x_s <= x if the ray is directed to the right, or for x_s >= x if it
   * is directed to the left.
   * \param pcoeffs The rational coefficients of the polynomial p(x).
   * \param qcoeffs The rational coefficients of the polynomial q(x).
   * \param x_s The x-coordinate of the source point.
   * \param dir_right Is the ray directed to the right (to +oo)
   *                  or to the left (to -oo).
   * \pre The denominator polynomial q(x) does not have any roots in the
   *      x-range of definition.
   */
  Continuous_rational_arc_d_1(const Polynomial_1& P,const Polynomial_1& Q,
                              const Algebraic_real_1& x_s, bool dir_right,
                              const Cache& cache) :
    Base(P, Q, x_s, dir_right,cache)
  {
    if (!this->_is_continuous())
    {
      // Invalid arc, as it is not continuous.
      this->set_invalid();
    }
  }

  Continuous_rational_arc_d_1(const Rat_vector& pcoeffs,
                              const Rat_vector& qcoeffs,
                              const Algebraic_real_1& x_s, bool dir_right,
                              const Cache& cache) :
    Base(pcoeffs, qcoeffs, x_s, dir_right,cache)
  {
    if (!this->_is_continuous())
    {
      // Invalid arc, as it is not continuous.
      this->set_invalid();
    }
  }

  /*!
   * Constructor of a bounded rational arc, defined by y = p(x)/q(x),
   * where: x_min <= x <= x_max.
   * \param pcoeffs The rational coefficients of the polynomial p(x).
   * \param qcoeffs The rational coefficients of the polynomial q(x).
   * \param x_s The x-coordinate of the source point.
   * \param x_t The x-coordinate of the target point.
   * \pre The two x-coordinates must not be equal.
   * \pre The denominator polynomial q(x) does not have any roots in the
   *      x-range of definition (x_min, x_max).
   */
  Continuous_rational_arc_d_1(const Polynomial_1& P,const Polynomial_1& Q,
                              const Algebraic_real_1& x_s,
                              const Algebraic_real_1& x_t,const Cache& cache) :
    Base(P, Q, x_s, x_t,cache)
  {
    if (!this->_is_continuous())
    {
      // Invalid arc, as it is not continuous.
      this->set_invalid();
    }
  }
  Continuous_rational_arc_d_1(const Rat_vector& pcoeffs,
                              const Rat_vector& qcoeffs,
                              const Algebraic_real_1& x_s,
                              const Algebraic_real_1& x_t, const Cache& cache) :
    Base(pcoeffs, qcoeffs, x_s, x_t,cache)
  {
    if (!this->_is_continuous())
    {
      // Invalid arc, as it is not continuous.
      this->set_invalid();
    }
  }

  //@}

  /// \name Constructions of points and curves.
  //@{

  /*! Compute the intersections with the given arc.
   * \param arc The given intersecting arc.
   * \param oi The output iterator.
   * \return The past-the-end iterator.
   */
  template <typename OutputIterator>
  OutputIterator intersect(const Self& arc, OutputIterator oi,
                           const Cache& cache) const
  {
    typedef boost::variant<Intersection_point, Self>    Intersection_result;

    CGAL_precondition(this->is_valid() && this->is_continuous());
    CGAL_precondition(arc.is_valid() && arc.is_continuous());

    if (this->equals(arc)) {
      Self overlap_arc(*this);
      *oi++ = Intersection_result(overlap_arc);
      return oi;
    }

    if (this->_has_same_base(arc)) {
      // Get the left and right endpoints of (*this) and their information
      // bits.
      const Algebraic_point_2& left1 =
        (this->is_directed_right() ? this->_ps : this->_pt);
      const Algebraic_point_2& right1 =
        (this->is_directed_right() ? this->_pt : this->_ps);
      int info_left1, info_right1;

      if (this->is_directed_right()) {
        info_left1 = (this->_info & this->SRC_INFO_BITS);
        info_right1 = ((this->_info & this->TRG_INFO_BITS) >> 4);
      }
      else {
        info_right1 = (this->_info & this->SRC_INFO_BITS);
        info_left1 = ((this->_info & this->TRG_INFO_BITS) >> 4);
      }

      // Get the left and right endpoints of the other arc and their
      // information bits.
      const Algebraic_point_2& left2 =
        (arc.is_directed_right() ? arc._ps : arc._pt);
      const Algebraic_point_2& right2 =
        (arc.is_directed_right() ? arc._pt : arc._ps);
      int info_left2, info_right2;

      if (arc.is_directed_right()) {
        info_left2 = (arc._info & this->SRC_INFO_BITS);
        info_right2 = ((arc._info & this->TRG_INFO_BITS) >> 4);
      }
      else {
        info_right2 = (arc._info & this->SRC_INFO_BITS);
        info_left2 = ((arc._info & this->TRG_INFO_BITS) >> 4);
      }

      // Locate the left curve-end with larger x-coordinate.
      bool at_minus_infinity = false;
      Arr_parameter_space inf_l1 = this->left_parameter_space_in_x();
      Arr_parameter_space inf_l2 = arc.left_parameter_space_in_x();
      Algebraic_point_2 p_left;
      int info_left;

      if (inf_l1 == ARR_INTERIOR && inf_l2 == ARR_INTERIOR) {
        // Let p_left be the rightmost of the two left endpoints.
        if (left1.x() > left2.x()) {
          p_left = left1;
          info_left = info_left1;
        }
        else {
          p_left = left2;
          info_left = info_left2;
        }
      }
      else if (inf_l1 == ARR_INTERIOR) {
        // Let p_left be the left endpoint of (*this).
        p_left = left1;
        info_left = info_left1;
      }
      else if (inf_l2 == ARR_INTERIOR) {
        // Let p_left be the left endpoint of the other arc.
        p_left = left2;
        info_left = info_left2;
      }
      else {
        // Both arcs are defined at x = -oo.
        at_minus_infinity = true;
        info_left = info_left1;
      }

      // Locate the right curve-end with smaller x-coordinate.
      bool at_plus_infinity = false;
      Arr_parameter_space inf_r1 = this->right_parameter_space_in_x();
      Arr_parameter_space inf_r2 = arc.right_parameter_space_in_x();
      Algebraic_point_2 p_right;
      int info_right;

      if (inf_r1 == ARR_INTERIOR && inf_r2 == ARR_INTERIOR) {
        // Let p_right be the rightmost of the two right endpoints.
        if (right1.x() < right2.x()) {
          p_right = right1;
          info_right = info_right1;
        }
        else {
          p_right = right2;
          info_right = info_right2;
        }
      }
      else if (inf_r1 == ARR_INTERIOR) {
        // Let p_right be the right endpoint of (*this).
        p_right = right1;
        info_right = info_right1;
      }
      else if (inf_r2 == ARR_INTERIOR) {
        // Let p_right be the right endpoint of the other arc.
        p_right = right2;
        info_right = info_right2;
      }
      else {
        // Both arcs are defined at x = +oo.
        at_plus_infinity = true;
        info_right = info_right2;
      }

      // Check the case of two bounded (in x) ends.
      if (! at_minus_infinity && ! at_plus_infinity) {
        Comparison_result res = CGAL::compare(p_left.x(), p_right.x());

        // The x-range of the overlap is empty, so there is no overlap.
        if (res == LARGER) return oi;

        if (res == EQUAL) {
          // We have a single overlapping point. Just make sure this point
          // is not at y = -/+ oo.
          if (info_left &&
              (this->SRC_AT_Y_MINUS_INFTY | this->SRC_AT_Y_PLUS_INFTY) == 0 &&
              info_right &&
              (this->SRC_AT_Y_MINUS_INFTY | this->SRC_AT_Y_PLUS_INFTY) == 0)
          {
            Intersection_point ip(p_left, 0);
            *oi++ = Intersection_result(ip);
          }

          return oi;
        }
      }

      // Create the overlapping portion of the rational arc by properly setting
      // the source (left) and target (right) endpoints and their information
      // bits.
      Self overlap_arc(*this);

      overlap_arc._ps = p_left;
      overlap_arc._pt = p_right;

      overlap_arc._info = ((info_left) | (info_right << 4) |
                           this->IS_DIRECTED_RIGHT | this->IS_CONTINUOUS |
                           this->IS_VALID);

      *oi++ = Intersection_result(overlap_arc);
      return oi;
    }

    // We wish to find the intersection points between:
    //
    //   y = p1(x)/q1(x)    and     y = p2(x)/q2(x)
    //
    // It is clear that the x-coordinates of the intersection points are
    // the roots of the polynomial: ip(x) = p1(x)*q2(x) - p2(x)*q1(x).

    Rational_function_pair rat_pair =
      this->get_rational_pair(this->_f, arc._f, cache);

    typename Algebraic_vector::const_iterator  x_iter;
    typename Multiplicity_vector::const_iterator  m_iter;

    // Go over the x-values we obtained. For each value produce an
    // intersection point if it is contained in the x-range of both curves.
    CGAL_precondition(rat_pair.roots().size() ==
                      rat_pair.multiplicities().size());
    for (x_iter = rat_pair.roots().begin(),
           m_iter = rat_pair.multiplicities().begin();
         x_iter != rat_pair.roots().end();
         ++x_iter, ++m_iter)
    {
      if (this->_is_in_true_x_range(*x_iter) && arc._is_in_true_x_range(*x_iter))
      {
        // Compute the intersection point and obtain its multiplicity.
        Algebraic_point_2 p(this->_f, *x_iter);
        // Output the intersection point:
        Intersection_point ip(p, *m_iter);
        *oi++ = Intersection_result(ip);
      }
    }

    return oi;
  }

  /*!
   * Split the arc into two at a given split point.
   * \param p The split point.
   * \param c1 Output: The first resulting arc, lying to the left of p.
   * \param c2 Output: The first resulting arc, lying to the right of p.
   * \pre p lies in the interior of the arc (not one of its endpoints).
   */
  void split(const Algebraic_point_2& p, Self& c1, Self& c2,
             const Cache& CGAL_assertion_code(cache)) const
  {
    CGAL_precondition(this->is_valid() && this->is_continuous());

    // Make sure that p lies on the interior of the arc.
    CGAL_precondition(this->point_position(p,cache) == EQUAL &&
                      (this->source_parameter_space_in_x() != ARR_INTERIOR ||
                       this->source_parameter_space_in_y() != ARR_INTERIOR ||
                       (p.x() != this->_ps.x()))  &&
                      (this->target_boundary_in_x() != ARR_INTERIOR ||
                       this->target_boundary_in_y() != ARR_INTERIOR ||
                       (p.x() != this->_pt.x())));

    // Make copies of the current arc.
    c1 = *this;
    c2 = *this;

    // Split the arc, such that c1 lies to the left of c2.
    if ((this->_info & this->IS_DIRECTED_RIGHT) != 0)
    {
      c1._pt = p;
      c1._info = (c1._info & ~this->TRG_INFO_BITS);
      c2._ps = p;
      c2._info = (c2._info & ~this->SRC_INFO_BITS);
    }
    else
    {
      c1._ps = p;
      c1._info = (c1._info & ~this->SRC_INFO_BITS);
      c2._pt = p;
      c2._info = (c2._info & ~this->TRG_INFO_BITS);
    }
  }

  /*!
   * Merge the current arc with the given arc.
   * \param arc The arc to merge with.
   * \pre The two arcs are mergeable.
   */
  void merge(const Self& arc)
  {
    CGAL_precondition(this->is_valid() && this->is_continuous());
    CGAL_precondition(arc.is_valid() && arc.is_continuous());
    CGAL_precondition(this->can_merge_with(arc));

    // Check if we should extend the arc to the left or to the right.
    if (this->right_parameter_space_in_x() == ARR_INTERIOR &&
        this->right_parameter_space_in_y() == ARR_INTERIOR &&
        arc.left_parameter_space_in_x() == ARR_INTERIOR &&
        arc.left_parameter_space_in_y() == ARR_INTERIOR &&
        (this->right().x() == arc.left().x()))
    {
      // Extend the arc to the right.
      if ((this->_info & this->IS_DIRECTED_RIGHT) != 0)
      {
        if (arc.right_parameter_space_in_x() == ARR_INTERIOR &&
            arc.right_parameter_space_in_y() == ARR_INTERIOR)
        {
          this->_pt = arc.right();
        }
        else
        {
          if (arc.right_parameter_space_in_x() == ARR_LEFT_BOUNDARY)
            this->_info = (this->_info | this->TRG_AT_X_MINUS_INFTY);
          else if (arc.right_parameter_space_in_x() == ARR_RIGHT_BOUNDARY)
            this->_info = (this->_info | this->TRG_AT_X_PLUS_INFTY);

          if (arc.right_parameter_space_in_y() == ARR_BOTTOM_BOUNDARY)
            this->_info = (this->_info | this->TRG_AT_Y_MINUS_INFTY);
          else if (arc.right_parameter_space_in_y() == ARR_TOP_BOUNDARY)
            this->_info = (this->_info | this->TRG_AT_Y_PLUS_INFTY);

          this->_pt = (arc._info & this->IS_DIRECTED_RIGHT) ? arc._pt : arc._ps;
        }
      }
      else
      {
        if (arc.right_parameter_space_in_x() == ARR_INTERIOR &&
            arc.right_parameter_space_in_y() == ARR_INTERIOR)
        {
          this->_ps = arc.right();
        }
        else
        {
          if (arc.right_parameter_space_in_x() == ARR_LEFT_BOUNDARY)
            this->_info = (this->_info | this->SRC_AT_X_MINUS_INFTY);
          else if (arc.right_parameter_space_in_x() == ARR_RIGHT_BOUNDARY)
            this->_info = (this->_info | this->SRC_AT_X_PLUS_INFTY);

          if (arc.right_parameter_space_in_y() == ARR_BOTTOM_BOUNDARY)
            this->_info = (this->_info | this->SRC_AT_Y_MINUS_INFTY);
          else if (arc.right_parameter_space_in_y() == ARR_TOP_BOUNDARY)
            this->_info = (this->_info | this->SRC_AT_Y_PLUS_INFTY);

          this->_ps = (arc._info & this->IS_DIRECTED_RIGHT) ? arc._pt : arc._ps;
        }
      }
    }
    else
    {
      CGAL_precondition(this->left_parameter_space_in_x() == ARR_INTERIOR &&
                        this->left_parameter_space_in_y() == ARR_INTERIOR &&
                        arc.right_parameter_space_in_x() == ARR_INTERIOR &&
                        arc.right_parameter_space_in_y() == ARR_INTERIOR &&
                        (this->left().x() == arc.right().x()));

      // Extend the arc to the left.
      if ((this->_info & this->IS_DIRECTED_RIGHT) != 0)
      {
        if (arc.left_parameter_space_in_x() == ARR_INTERIOR &&
            arc.left_parameter_space_in_y() == ARR_INTERIOR)
        {
          this->_ps = arc.left();
        }
        else
        {
          if (arc.left_parameter_space_in_x() == ARR_LEFT_BOUNDARY)
            this->_info = (this->_info | this->SRC_AT_X_MINUS_INFTY);
          else if (arc.left_parameter_space_in_x() == ARR_RIGHT_BOUNDARY)
            this->_info = (this->_info | this->SRC_AT_X_PLUS_INFTY);

          if (arc.left_parameter_space_in_y() == ARR_BOTTOM_BOUNDARY)
            this->_info = (this->_info | this->SRC_AT_Y_MINUS_INFTY);
          else if (arc.left_parameter_space_in_y() == ARR_TOP_BOUNDARY)
            this->_info = (this->_info | this->SRC_AT_Y_PLUS_INFTY);

          this->_ps = (arc._info & this->IS_DIRECTED_RIGHT) ? arc._ps : arc._pt;
        }
      }
      else
      {
        if (arc.left_parameter_space_in_x() == ARR_INTERIOR &&
            arc.left_parameter_space_in_y() == ARR_INTERIOR)
        {
          this->_pt = arc.left();
        }
        else
        {
          if (arc.left_parameter_space_in_x() == ARR_LEFT_BOUNDARY)
            this->_info = (this->_info | this->TRG_AT_X_MINUS_INFTY);
          else if (arc.left_parameter_space_in_x() == ARR_RIGHT_BOUNDARY)
            this->_info = (this->_info | this->TRG_AT_X_PLUS_INFTY);

          if (arc.left_parameter_space_in_y() == ARR_BOTTOM_BOUNDARY)
            this->_info = (this->_info | this->TRG_AT_Y_MINUS_INFTY);
          else if (arc.left_parameter_space_in_y() == ARR_TOP_BOUNDARY)
            this->_info = (this->_info | this->TRG_AT_Y_PLUS_INFTY);

          this->_pt = (arc._info & this->IS_DIRECTED_RIGHT) ? arc._ps : arc._pt;
        }
      }
    }
  }
  //@}


};

//*! \class Rational_arc_2
  // * Representation of a generic, not necessarily continuous, portion of a
  // * rational function.
  // */
template <typename Algebraic_kernel_>
class Rational_arc_d_1 : public Base_rational_arc_d_1<Algebraic_kernel_>
{
public:
  typedef Algebraic_kernel_                             Algebraic_kernel;

  typedef Rational_arc_d_1<Algebraic_kernel>            Self;
  typedef Base_rational_arc_d_1<Algebraic_kernel>       Base;
  typedef Continuous_rational_arc_d_1<Algebraic_kernel> Continuous_arc;

  typedef typename Base::Integer                        Integer;
  typedef typename Base::Rational                       Rational;
  typedef typename Base::Algebraic_real_1               Algebraic_real_1;
  typedef typename Base::Algebraic_point_2              Algebraic_point_2;
  typedef typename Base::Polynomial_1                   Polynomial_1;

  typedef typename Base::Rat_vector                     Rat_vector;

  typedef typename Base::Cache                          Cache;

  /// \name Constrcution methods.
  //@{

  /*!
   * Default constructor.
   */
  Rational_arc_d_1() :
    Base()
  {}

  /*!
   * Constructor of a whole polynomial curve.
   * \param pcoeffs The rational coefficients of the polynomial p(x).
   */
  Rational_arc_d_1(const Polynomial_1& P, const Cache& cache) :
    Base(P, cache)
  {}

  Rational_arc_d_1(const Rat_vector& pcoeffs, const Cache& cache) :
    Base(pcoeffs, cache)
  {}

  /*!
   * Constructor of a polynomial ray, defined by y = p(x), for x_s <= x if the
   * ray is directed to the right, or for x_s >= x if it is directed to the
   * left.
   * \param pcoeffs The rational coefficients of the polynomial p(x).
   * \param x_s The x-coordinate of the source point.
   * \param dir_right Is the ray directed to the right (to +oo)
   *                  or to the left (to -oo).
   */
  Rational_arc_d_1(const Polynomial_1& P, const Algebraic_real_1& x_s,
                   bool dir_right, const Cache& cache) :
    Base(P, x_s, dir_right,cache)
  {}

  Rational_arc_d_1(const Rat_vector& pcoeffs, const Algebraic_real_1& x_s,
                   bool dir_right, const Cache& cache) :
    Base(pcoeffs, x_s, dir_right,cache)
  {}


  /*!
   * Constructor of a polynomial arc, defined by y = p(x), x_min <= x <= x_max.
   * \param pcoeffs The rational coefficients of the polynomial p(x).
   * \param x_s The x-coordinate of the source point.
   * \param x_t The x-coordinate of the target point.
   * \pre The two x-coordinates must not be equal.
   */
  Rational_arc_d_1(const Polynomial_1& P, const Algebraic_real_1& x_s,
                   const Algebraic_real_1& x_t, const Cache& cache) :
    Base(P, x_s, x_t, cache)
  {}

  Rational_arc_d_1(const Rat_vector& pcoeffs, const Algebraic_real_1& x_s,
                   const Algebraic_real_1& x_t,const Cache& cache) :
    Base(pcoeffs, x_s, x_t, cache)
  {}

  /*!
   * Constructor of a polynomial function, defined by y = p(x)/q(x) for any x.
   * \param pcoeffs The rational coefficients of the polynomial p(x).
   * \param qcoeffs The rational coefficients of the polynomial q(x).
   */
  Rational_arc_d_1(const Polynomial_1& P,const Polynomial_1& Q, const Cache& cache) :
    Base(P, Q,cache)
  {}

  Rational_arc_d_1(const Rat_vector& pcoeffs, const Rat_vector& qcoeffs,
                   const Cache& cache) :
    Base(pcoeffs, qcoeffs, cache)
  {}

  /*!
   * Constructor of a ray of a rational function, defined by y = p(x)/q(x),
   * for x_s <= x if the ray is directed to the right, or for x_s >= x if it
   * is directed to the left.
   * \param pcoeffs The rational coefficients of the polynomial p(x).
   * \param qcoeffs The rational coefficients of the polynomial q(x).
   * \param x_s The x-coordinate of the source point.
   * \param dir_right Is the ray directed to the right (to +oo)
   *                  or to the left (to -oo).
   */
  Rational_arc_d_1(const Polynomial_1& P,const Polynomial_1& Q,
                   const Algebraic_real_1& x_s, bool dir_right, const Cache& cache) :
    Base(P, Q, x_s, dir_right, cache)
  {}

  Rational_arc_d_1(const Rat_vector& pcoeffs, const Rat_vector& qcoeffs,
                   const Algebraic_real_1& x_s, bool dir_right, const Cache& cache) :
    Base(pcoeffs, qcoeffs, x_s, dir_right, cache)
  {}

  /*!
   * Constructor of a bounded rational arc, defined by y = p(x)/q(x),
   * where: x_min <= x <= x_max.
   * \param pcoeffs The rational coefficients of the polynomial p(x).
   * \param qcoeffs The rational coefficients of the polynomial q(x).
   * \param x_s The x-coordinate of the source point.
   * \param x_t The x-coordinate of the target point.
   * \pre The two x-coordinates must not be equal.
   */
  Rational_arc_d_1(const Polynomial_1& P,const Polynomial_1& Q,
                   const Algebraic_real_1& x_s, const Algebraic_real_1& x_t,
                   const Cache& cache) :
    Base(P, Q, x_s, x_t, cache)
  {}

  Rational_arc_d_1(const Rat_vector& pcoeffs, const Rat_vector& qcoeffs,
                   const Algebraic_real_1& x_s, const Algebraic_real_1& x_t,
                   const Cache& cache) :
    Base(pcoeffs, qcoeffs, x_s, x_t, cache)
  {}
  //@}

  /*!
   * Subdivide the given portion of a rational function into continuous
   * sub-arcs, splitting it at the roots of the denominator polynomial.
   * \param oi An output iterator of Continuous_rational_arc_d_1 objects.
   */
  template <typename OutputIterator>
  OutputIterator make_continuous(OutputIterator oi) const
  {
    // Compute the roots of the denominator polynomial.
    std::list<Algebraic_real_1>          q_roots;
    bool                                 root_at_ps, root_at_pt;

    if ((this->_info & this->IS_CONTINUOUS) == 0)
      this->_denominator_roots(std::back_inserter(q_roots),
                               root_at_ps, root_at_pt);

    // Check the case of a continuous arc:
    Base    arc = *this;

    if (q_roots.empty())
    {
      arc.set_continuous();
      *oi++ = Continuous_arc(arc);
      return (oi);
    }

    // The denominator has roots: split the arc accordingly.
    typename std::list<Algebraic_real_1>::const_iterator iter;

    for (iter = q_roots.begin(); iter != q_roots.end(); ++iter)
    {
      *oi++ = Continuous_arc(arc.split_at_pole(*iter));
    }

    // Add the final x-monotone sub-arc.
    arc.set_continuous();
    *oi++ = Continuous_arc(arc);

    return (oi);
  }
};

}   //name_space Arr_rational_arc
}       //namespace CGAL {

#endif //CGAL_RATIONAL_ARC_D_1_H
