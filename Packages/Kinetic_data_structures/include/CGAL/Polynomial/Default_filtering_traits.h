// Copyright (c) 2005  Stanford University (USA).
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_DEFAULT_FILTERING_TRAITS_H
#define CGAL_POLYNOMIAL_DEFAULT_FILTERING_TRAITS_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/Interval_polynomial.h>
#include <CGAL/Polynomial/polynomial_converters.h>

#ifdef CGAL_POLYNOMIAL_USE_CGAL
#include <CGAL/Gmpq.h>
#define CGAL_DEFAULT_FILTERING_DEFAULT_NT =CGAL::Gmpq
#endif

CGAL_POLYNOMIAL_BEGIN_NAMESPACE
template <class NT CGAL_DEFAULT_FILTERING_DEFAULT_NT>
struct Default_filtering_traits {
  typedef Polynomial<NT> Exact_function;
  typedef Interval_polynomial Interval_function;
  typedef Polynomial_converter<Exact_function, 
			       Interval_function, 
			       To_interval<NT> > Exact_to_interval_converter;
};

CGAL_POLYNOMIAL_END_NAMESPACE

#ifdef CGAL_DEFAULT_FILTERING_DEFAULT_NT
#undef CGAL_DEFAULT_FILTERING_DEFAULT_NT
#endif
#endif
