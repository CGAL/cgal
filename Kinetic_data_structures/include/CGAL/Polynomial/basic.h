// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_BASIC_H
#define CGAL_POLYNOMIAL_BASIC_H

#include <CGAL/Polynomial/internal/config.h>
#include <CGAL/Polynomial/internal/Statistics.h>
#include <CGAL/assertions.h>
/*!

  \file CGAL/Polynomial/basic.h The file which defines the basic
  things needed throught the module.  All other headers should include
  this first.

*/

#ifdef CGAL_USE_GMP
#ifndef CGAL_POLYNOMIAL_USE_CGAL
#define CGAL_POLYNOMIAL_USE_CGAL
#endif
#endif

#include <CGAL/Polynomial/internal/macros.h>

#ifdef CGAL_POLYNOMIAL_USE_CGAL
/*
  When CGAL is present
*/
#include <CGAL/basic.h>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpq.h>
namespace CGAL { namespace POLYNOMIAL {
typedef CGAL::Gmpq Default_field_nt;
} } //namespace CGAL::POLYNOMIAL
#else
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
namespace CGAL { namespace POLYNOMIAL {
typedef CGAL::Quotient<CGAL::MP_Float> Default_field_nt;
} } //namespace CGAL::POLYNOMIAL
#endif

namespace CGAL { namespace POLYNOMIAL {

/*typedef CGAL::Sign Sign;
static const Sign ZERO= CGAL::ZERO;
static const Sign POSITIVE= CGAL::POSITIVE;
static const Sign NEGATIVE= CGAL::NEGATIVE;*/
/*typedef int Comparison_result;
static const int EQUAL= CGAL::EQUAL;
static const int SMALLER= CGAL::SMALLER;
static const int LARGER = CGAL::LARGER;*/
static const int UNKNOWN = -3;
typedef int Order;
static const int STRICTLY_BELOW = -3;
static const int BELOW=-2;
static const int CONTAINED=-1;
static const int CONTAINS=1;
static const int ABOVE=2;
static const int STRICTLY_ABOVE=3;
//typedef enum Sign {ZERO=CGAL::ZERO, POSITIVE=CGAL::POSITIVE, NEGATIVE=CGAL::NEGATIVE} Sign;

//typedef CGAL::Comparison_result Comparison_result;

/*template <class NT>
Sign sign(const NT &nt)
{
    return CGAL::sign(nt);
    }*/


/*typedef ::CGAL::Integral_domain_without_division_tag             Integral_domain_without_division_tag;
typedef ::CGAL::Euclidean_ring_tag   Euclidean_ring_tag;
typedef ::CGAL::Field_tag            Field_tag;
typedef ::CGAL::Field_with_sqrt_tag       Field_with_sqrt_tag;*/

} } //namespace CGAL::POLYNOMIAL

#define CGAL_POLYNOMIAL_TO_DOUBLE(d) CGAL::to_double(d)

#define CGAL_POLYNOMIAL_TO_INTERVAL(d) CGAL::to_interval(d)

#else
/*
  When no CGAL is present
*/

namespace CGAL { namespace POLYNOMIAL {

typedef int Sign;
static const int ZERO =0;
static const int POSITIVE=1;
static const int NEGATIVE=-1;
typedef int Comparison_result;
static const int EQUAL= 0;
static const int SMALLER= -1;
static const int LARGER = 1;

template <class NT>
Sign sign(const NT &nt)
{
    if (nt >0) return POSITIVE;
    else if (nt <0) return NEGATIVE;
    else return ZERO;
}


struct Integral_domain_without_division_tag {};
struct Euclidean_ring_tag {};
struct Field_tag {};
struct Field_with_sqrt_tag {};

} } //namespace CGAL::POLYNOMIAL
#endif

#include <limits>

/*
  Shared
*/

namespace CGAL { namespace POLYNOMIAL {

//! Used for signs when filtering is involved
typedef enum Extended_sign
{
    EXTENDED_NEGATIVE=NEGATIVE, EXTENDED_ZERO=ZERO,
    EXTENDED_POSITIVE=POSITIVE, EXTENDED_UNKNOWN=2
} Extended_sign;

//! Calculate the extened sign of a number.
/*!  For a regular number, the extened sign is teh same as the sign.

There is a specializated version for filtering included in
  CGAL/Polynomial/Tools/interval_arithmetic.h
*/
template <class NT>
inline Extended_sign extended_sign(const NT &nt)
{
  // for VC
  switch(CGAL::sign(nt)) {
        case ZERO: return EXTENDED_ZERO;
        case POSITIVE: return EXTENDED_POSITIVE;
        default: return EXTENDED_NEGATIVE;
    };
}


template <class Rt>
inline Rt infinity()
{
  //CGAL_static_assertion(std::numeric_limits<Rt>::is_specialized);
    if (std::numeric_limits<Rt>::has_infinity) return std::numeric_limits<Rt>::infinity();
    else return (std::numeric_limits<Rt>::max)();
}


} } //namespace CGAL::POLYNOMIAL

/*!  \namespace CGAL::POLYNOMIAL This is the namespace where all the
  directly exposed members of the framework are defined.

  \todo put specialized templates first in header

  \todo wrapper support for boost::interval

  \namespace CGAL::POLYNOMIAL::internal The classes in this namespace
  should not be accessed directly by the users of the
  library. However, many of them are exposed in other ways (for
example as members of a kernel), so things need to be documented.

\namespace CGAL We don't put much in CGAL. Just to_double,
to_interval, sign and a couple of other methods.
*/

/*
namespace CGAL { namespace POLYNOMIAL {

//!\todo are these needed? The tags in the polynomial NS?
struct Descartes_tag {};
struct Sturm_tag {};
struct Bezier_tag {};

} } //namespace CGAL::POLYNOMIAL
*/

#include <limits>

#endif // CGAL_POLYNOMIAL_BASIC_H
