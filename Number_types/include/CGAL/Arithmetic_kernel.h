// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
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
// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>
//
// ============================================================================
//
//    \brief provide class Arithmetic_kernel, a collection of number types. 
//

// 

/*! \file CGAL/Arithmetic_kernel.h
 *  \brief Declarations pertaining to CGAL::Arithmetic_kernel
 */

#ifndef CGAL_ARITHMETIC_KERNEL_H
#define CGAL_ARITHMETIC_KERNEL_H

#include <CGAL/basic.h>

//#include <CGAL/Polynomial.h> // TODO: Not available yet.
//#include <CGAL/Algebraic_real.h> // TODO: Not available yet.

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_bigfloat.h>
#include <CGAL/leda_real.h>
#include <CGAL/leda_bigfloat_interval.h>

#endif // CGAL_USE_LEDA

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpzf.h>
#include <CGAL/Gmpq.h>

#ifdef CGAL_USE_GMPXX
#include <CGAL/gmpxx.h>
#endif // CGAL_USE_GMPXX
#endif // CGAL_USE_GMP

#ifdef CGAL_USE_CORE
#include <CGAL/CORE_BigFloat.h>
#include <CGAL/CORE_BigInt.h>
#include <CGAL/CORE_BigRat.h>
#include <CGAL/CORE_Expr.h>
#endif // CGAL_USE_CORE

#include <CGAL/Sqrt_extension.h>


CGAL_BEGIN_NAMESPACE

namespace CGALi {

class Arithmetic_kernel_base{
public:
    typedef CGAL::Null_tag Integer;
    typedef CGAL::Null_tag Rational;
    typedef CGAL::Null_tag Field_with_sqrt;
    typedef CGAL::Null_tag Field_with_kth_root;
    typedef CGAL::Null_tag Field_with_root_of;
    typedef CGAL::Null_tag Bigfloat;
    typedef CGAL::Null_tag Bigfloat_interval;
};

}
#ifdef CGAL_USE_LEDA 
/*! \ingroup CGAL_Arithmetic_kernel
 *  \brief  The LEDA set of exact number types
 */
class LEDA_arithmetic_kernel : public CGALi::Arithmetic_kernel_base {
public:
    //! exact integers
    typedef leda_integer Integer;
    //! exact float nummber
    typedef leda_bigfloat Exact_float_number;
    //! exact rationals, constructible from integers
    typedef leda_rational Rational;
    //! exact root expressions, constructible from integers and rationals
    typedef leda_real Field_with_sqrt;

    // undocumented
    typedef leda_bigfloat          Bigfloat;
    typedef leda_bigfloat_interval Bigfloat_interval;

};
#endif // CGAL_USE_LEDA

#ifdef CGAL_USE_CORE
/*! \ingroup CGAL_Arithmetic_kernel
 *  \brief  The CORE set of exact number types
 */
class CORE_arithmetic_kernel : public CGALi::Arithmetic_kernel_base {
public:
    //! exact integers
    typedef CORE::BigInt Integer;
    //! exact float nummber
    typedef CORE::BigRat Exact_float_number;
    //! exact rationals, constructible from integers
    typedef CORE::BigRat Rational;
    //! exact root expressions, constructible from integers and rationals
    typedef CORE::Expr Field_with_sqrt;
    // undocumented
    //typedef CORE::BigFloat          Bigfloat;
    typedef CORE::BigFloat          Bigfloat_interval;

};
#endif // CGAL_USE_CORE

#ifdef CGAL_USE_GMP
/*! \ingroup CGAL_Arithmetic_kernel
 *  \brief  The GMP set of exact number types
 */
class GMP_arithmetic_kernel : public CGALi::Arithmetic_kernel_base {
public:
    //! exact integers
    typedef CGAL::Gmpz Integer;
    //! exact float nummber
    typedef CGAL::Gmpq Exact_float_number;
    //! exact rationals, constructible from integers
    typedef CGAL::Gmpq Rational;
    //! exact root expressions, constructible from integers and rationals
    typedef CGAL::Null_tag  Field_with_sqrt;
    typedef CGAL::Null_tag  Field_with_kth_root;
    typedef CGAL::Null_tag  Field_with_root_of;
};
#endif // CGAL_USE_GMP

#ifdef CGAL_USE_GMPXX
/*! \ingroup CGAL_Arithmetic_kernel
 *  \brief  The GMP set of exact number types
 */
class GMPXX_arithmetic_kernel : public CGALi::Arithmetic_kernel_base {
public:
    //! exact integers
    typedef mpz_class Integer;
    //! exact float nummber
    typedef mpq_class Exact_float_number;
    //! exact rationals, constructible from integers
    typedef mpq_class Rational;
    //! exact root expressions, constructible from integers and rationals
    typedef CGAL::Null_tag  Field_with_sqrt;
};
#endif // CGAL_USE_GMPXX

// Define a default Arithmetic_kernel
#if defined(CGAL_USE_CORE)
typedef CORE_arithmetic_kernel Arithmetic_kernel;
#else
#if defined(CGAL_USE_LEDA)
typedef LEDA_arithmetic_kernel Arithmetic_kernel;
#endif
#endif

#if defined(CGAL_USE_LEDA) || defined(CGAL_USE_CORE)
#ifndef CGAL_HAVE_DEFAULT_ARITHMETIC_KERNEL
#define CGAL_HAVE_DEFAULT_ARITHMETIC_KERNEL 1
#endif
#endif // defined(CGAL_USE_LEDA) || defined(CGAL_USE_CORE)

// Macro to snap typedefs in Arithmetic_kernel
#define CGAL_SNAP_ARITHMETIC_KERNEL_TYPEDEFS(AT) \
  typedef typename AT::Integer Integer; \
  typedef typename AT::Rational Rational; \
  typedef typename AT::Field_with_sqrt Field_with_sqrt; 

// end #define


template< class NT > struct Get_arithmetic_kernel;

#ifdef CGAL_USE_LEDA
    
    template <>
    struct Get_arithmetic_kernel<leda::integer> {
        typedef LEDA_arithmetic_kernel Arithmetic_kernel;
    };
    template <>
    struct Get_arithmetic_kernel<leda::rational>{
        typedef LEDA_arithmetic_kernel Arithmetic_kernel;
    };
    template <>
    struct Get_arithmetic_kernel<leda::real>{
        typedef LEDA_arithmetic_kernel Arithmetic_kernel;
    };
    template <>
    struct Get_arithmetic_kernel<leda::bigfloat>{
        typedef LEDA_arithmetic_kernel Arithmetic_kernel;
    };
    template <>
    struct Get_arithmetic_kernel<CGAL::leda_bigfloat_interval>{
        typedef LEDA_arithmetic_kernel Arithmetic_kernel;
    };
#endif //CGAL_USE_LEDA
    
    
#ifdef CGAL_USE_CORE
    
    template <>
    struct Get_arithmetic_kernel<CORE::BigInt>{
        typedef CORE_arithmetic_kernel Arithmetic_kernel;
    };
    template <>
    struct Get_arithmetic_kernel<CORE::BigRat>{
        typedef CORE_arithmetic_kernel Arithmetic_kernel;
    };
    template <>
    struct Get_arithmetic_kernel<CORE::Expr>{
        typedef CORE_arithmetic_kernel Arithmetic_kernel;
    };
    template <>
    struct Get_arithmetic_kernel<CORE::BigFloat>{
        typedef CORE_arithmetic_kernel Arithmetic_kernel;
    };
#endif //CGAL_USE_CORE
    
    template <class COEFF, class ROOT>
    struct Get_arithmetic_kernel<Sqrt_extension<COEFF,ROOT> >{
        typedef Get_arithmetic_kernel<COEFF> GET;
        typedef typename GET::Arithmetic_kernel Arithmetic_kernel;
    };
    

CGAL_END_NAMESPACE

#endif // CGAL_ARITHMETIC_KERNEL_H
// EOF
