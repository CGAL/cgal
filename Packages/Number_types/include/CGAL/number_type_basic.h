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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Stefan Schirra
 

#ifndef CGAL_NUMBER_TYPE_BASIC_H
#define CGAL_NUMBER_TYPE_BASIC_H

#define CGAL_PI 3.14159265358979323846

#ifdef CGAL_USE_ADL_FOR_NT
// Attempt at using Koenig lookup for the NT interface.
#  define CGAL_NTS
#else
#  define CGAL_NTS CGAL::
#endif
// #define CGAL_NTS CGAL::NTS::

#if ((__GNUC__ == 2) && (__GNUC_MINOR__ == 95))
#include <cmath>
#endif

// CGAL uses std::min and std::max

#include <algorithm>

CGAL_BEGIN_NAMESPACE

using std::min;
using std::max;

CGAL_END_NAMESPACE

#include <CGAL/Number_type_traits.h>
#include <CGAL/number_utils.h>
#include <CGAL/double.h>
#include <CGAL/float.h>
#include <CGAL/int.h>

// Including all number type files is necessary for compilers implementing
// two-stage name lookup (like g++ >= 3.4).
// A nicer solution needs more thought.

#ifdef CGAL_CFG_HAS_TWO_STAGE_NAME_LOOKUP

#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Fixed_precision_nt.h>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#endif

#ifdef CGAL_USE_GMPXX
#include <CGAL/gmpxx.h>
#endif

#ifdef CGAL_USE_CORE
#include <CGAL/CORE_Expr.h>
#endif

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_bigfloat.h>
#include <CGAL/leda_real.h>
#endif

#endif // CGAL_CFG_HAS_TWO_STAGE_NAME_LOOKUP

#include <CGAL/number_utils_classes.h>

#endif // CGAL_NUMBER_TYPE_BASIC_H
