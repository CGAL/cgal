/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CORE (http://cs.nyu.edu/exact/core/).
 * You can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * File: Expr.h
 * Synopsis: a class of Expression in Level 3
 *
 * Written by
 *       Koji Ouchi <ouchi@simulation.nyu.edu>
 *       Chee Yap <yap@cs.nyu.edu>
 *       Igor Pechtchanski <pechtcha@cs.nyu.edu>
 *       Vijay Karamcheti <vijayk@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *       Sylvain Pion <pion@cs.nyu.edu>
 *       Vikram Sharma<sharma@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Source: /home/exact/cvsroot/exact/corelib/inc/CORE/Promote.h,v $
 * $Revision$ $Date$
 ***************************************************************************/

#ifndef __PROMOTE_H__
#define __PROMOTE_H__

#include <CGAL/CORE/Impl.h>

namespace CORE { 

/// hasExactDivision()
///   CHECKING if NT has exact division
///   NOTE: it is important that the compiler does not try to
///   prove theorems about arithmetic identities like "x*(y/x) == y"
///   USAGE:  If you want to check if a number type NT has exact division, do for example,
///            if (hasExactDivision< NT >::check()) ...
///   		We use this in Polynomial<NT> class.
template < class NT >
struct hasExactDivision {
  static bool check() {		// This default function is supposed to work for NT other than BigRat or Expr
     return false;
  }
};

template<> struct hasExactDivision<Expr> {
  static bool check() {
     return true;
  }
};
template<> struct hasExactDivision<BigRat> {
  static bool check() {
     return true;
  }
};

template<typename T1, typename T2>
class Promotion;

template<typename T>
class Promotion<T, T> {
  public:
    typedef T ResultT;
};

#define MAX_TYPE(T1, T2)		\
  typename Promotion<T1, T2>::ResultT

#define DEFINE_MAX_TYPE(T1, T2, Tr)	\
  template<> class Promotion<T1, T2> {	\
    public:				\
      typedef Tr ResultT;		\
  };					\
  template<> class Promotion<T2, T1> {	\
    public:				\
      typedef Tr ResultT;		\
  };

/*
 * For example:
 *
 * DEFINE_MAX_TYPE(BigInt, BigRat, BigRat)   	// define the promotion
 *
 * template<typename T1, typename T2> 		// define function f with type templates
 *   MAX_TYPE(T1, T2) f(T1& , T2& );
 *
 * or
 *
 * template<typename T1, typename T2> 		// define function f with type templates
 *   const MAX_TYPE(T1, T2)& f(T1& , T2& );
 *
 * BigInt  a  =  1;
 * BigRat  b  = "1/3";
 * BigRat  c  =  f(a, b);			// or, typename Promotion<BigInt, BigRat>::ResultT c = f(a,b);
 *
 * REMARK: this mechanism is used by the eval function for polynomial evaluation (see Poly.tcc)
 * where the two types are NT (type of coefficients) and N (type of evaluation point).
 */

/* 
 * primary types: (11)
 *
 * 	bool, 
 *	char, unsigned char, 
 *	short, unsigned short,
 * 	int, unsigned int, 
 *	long, unsigned long, 
 *	float, double
 *
 * CORE types: (5)
 *
 * 	BigInt < BigFloat < BigRat < Real < Expr
 *
 *      (NOTE: BigFloat here must be error-free)
 *
 */

class BigInt;
class BigFloat;
class BigRat;
class Expr;

DEFINE_MAX_TYPE(long, BigInt, BigInt)
DEFINE_MAX_TYPE(long, BigFloat, BigFloat)
DEFINE_MAX_TYPE(long, BigRat, BigRat)
DEFINE_MAX_TYPE(long, Expr, Expr)

DEFINE_MAX_TYPE(int, BigInt, BigInt)
DEFINE_MAX_TYPE(int, BigFloat, BigFloat)
DEFINE_MAX_TYPE(int, BigRat, BigRat)
DEFINE_MAX_TYPE(int, Expr, Expr)

DEFINE_MAX_TYPE(BigInt, BigFloat, BigFloat)
DEFINE_MAX_TYPE(BigInt, BigRat, BigRat)
DEFINE_MAX_TYPE(BigInt, Expr, Expr)

DEFINE_MAX_TYPE(BigFloat, BigRat, BigRat)
DEFINE_MAX_TYPE(BigFloat, Expr, Expr)

DEFINE_MAX_TYPE(BigRat, Expr, Expr)

} //namespace CORE

#endif //__PROMOTE_H__
