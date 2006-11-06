// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Stefan Schirra
 

#ifndef CGAL_NUMBER_TYPE_BASIC_H
#define CGAL_NUMBER_TYPE_BASIC_H

#define CGAL_PI 3.14159265358979323846


#ifdef CGAL_USE_NTS_NAMESPACE

#define CGAL_NTS_BEGIN_NAMESPACE namespace NTS {
#define CGAL_NTS_END_NAMESPACE }
#define CGAL_NTS ::CGAL::NTS::

#else

#define CGAL_NTS_BEGIN_NAMESPACE
#define CGAL_NTS_END_NAMESPACE
#define CGAL_NTS ::CGAL::

#endif

// CGAL uses std::min and std::max
#include <algorithm>

#include <CGAL/config.h>

CGAL_BEGIN_NAMESPACE

using std::min;
using std::max;

CGAL_END_NAMESPACE

#include <CGAL/functional_base.h> // Unary_function, Binary_function

#include <CGAL/Kernel/mpl.h>
#include <CGAL/known_bit_size_integers.h>
#include <CGAL/enum.h>

#include <CGAL/Coercion_traits.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/Number_type_traits.h> 

#include <CGAL/Fraction_traits.h>
//#include <CGAL/Scalar_factor_traits.h>
//#include <CGAL/Algebraic_number_traits.h>


#include <CGAL/utils_classes.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/tags.h>
#include <cmath>

//#include <CGAL/Algebraic_structure_traits.h>
//#include <CGAL/Real_embeddable_traits.h>
//#include <CGAL/number_utils_fwd.h>

#include <CGAL/float.h>
#include <CGAL/double.h>
#include <CGAL/long_double.h>
#include <CGAL/int.h>
#ifdef CGAL_USE_LONG_LONG
#include <CGAL/long_long.h>
#endif

#include <CGAL/FPU.h>

#include <CGAL/Interval_nt_fwd.h>

// Including all number type files is necessary for compilers implementing
// two-stage name lookup (like g++ >= 3.4).
// A nicer solution needs more thought.

#ifndef CGAL_CFG_NO_TWO_STAGE_NAME_LOOKUP

#include <CGAL/Interval_nt_fwd.h>
#include <CGAL/Lazy_exact_nt_fwd.h>
//#include <CGAL/MP_Float_fwd.h>
#include <CGAL/Nef_polynomial_fwd.h>
#include <CGAL/Number_type_checker_fwd.h>

//#ifdef CGAL_USE_GMP
//#include <CGAL/Gmpzq_fwd.h>
//#endif

#ifdef CGAL_USE_GMPXX
#  include <CGAL/gmpxx_fwd.h>
#endif


#include <CGAL/Quotient_fwd.h>
#include <CGAL/Root_of_2_fwd.h>

#include <CGAL/make_root_of_2.h>

// We must also include the following three, because of the overloadings
// for Quotient<MP_Float> and Quotient<Gmpz/Gmpzf>, which triggers their
// instantiation, even if only to_double(double) is called, at least
// when Quotient is defined...

//#include <CGAL/MP_Float.h>
#ifdef CGAL_USE_GMP
//#include <CGAL/Gmpz.h>
//#include <CGAL/Gmpq.h>
//#include <CGAL/Gmpzf.h>
#endif

CGAL_BEGIN_NAMESPACE

#if 0
// Polynomial

template <typename T> class Polynomial;

template <typename ET>
double to_double(const Polynomial<ET> &);

template <typename ET>
std::pair<double,double> to_interval(const Polynomial<ET> &);

template <typename ET>
Sign sign(const Polynomial<ET> &);


template <typename ET>
Polynomial<ET> abs(const Polynomial<ET> &);

template <typename ET>
bool is_finite(const Polynomial<ET> &);

template <typename ET>
bool is_valid(const Polynomial<ET> &);

template <typename ET>
Polynomial<ET> gcd(const Polynomial<ET> &, const Polynomial<ET> &);

// Nef_polynomial

template <typename T> class Nef_polynomial;

template <typename ET>
double to_double(const Nef_polynomial<ET> &);

template <typename ET>
Nef_polynomial<ET> gcd(const Nef_polynomial<ET> &, const Nef_polynomial<ET> &);
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CFG_NO_TWO_STAGE_NAME_LOOKUP


#include <CGAL/number_utils.h>
#include <CGAL/number_utils_classes.h>

#endif // CGAL_NUMBER_TYPE_BASIC_H
