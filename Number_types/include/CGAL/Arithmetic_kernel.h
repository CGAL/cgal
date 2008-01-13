// TODO: Add licence
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

#ifdef CGAL_INTERN_USE_BFI
#include <CGAL/Algebraic_kernel_d/leda_interval_support.h>
#endif //CGAL_INTERN_USE_BFI

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

#ifdef CGAL_USE_LEDA 
/*! \ingroup CGAL_Arithmetic_kernel
 *  \brief  The LEDA set of exact number types
 */
class LEDA_arithmetic_kernel {
public:
    //! exact integers
    typedef leda_integer Integer;
    //! exact float nummber
    typedef leda_bigfloat Exact_float_number;
    //! exact rationals, constructible from integers
    typedef leda_rational Rational;
    //! exact root expressions, constructible from integers and rationals
    typedef leda_real Field_with_sqrt;

#ifdef CGAL_INTERN_USE_BFI
    // undocumented
    typedef leda_bigfloat          Bigfloat;
    typedef CGALi::leda_bigfloat_interval Bigfloat_interval;
#endif //CGAL_INTERN_USE_BFI    

};
#endif // CGAL_USE_LEDA

#ifdef CGAL_USE_CORE
/*! \ingroup CGAL_Arithmetic_kernel
 *  \brief  The CORE set of exact number types
 */
class CORE_arithmetic_kernel {
public:
    //! exact integers
    typedef CORE::BigInt Integer;
    //! exact float nummber
    typedef CORE::BigRat Exact_float_number;
    //! exact rationals, constructible from integers
    typedef CORE::BigRat Rational;
    //! exact root expressions, constructible from integers and rationals
    typedef CORE::Expr Field_with_sqrt;
#ifdef CGAL_INTERN_USE_BFI
    // undocumented
    typedef CORE::BigFloat          Bigfloat;
    typedef CORE::BigFloat          Bigfloat_interval;
#endif //CGAL_INTERN_USE_BFI

};
#endif // CGAL_USE_CORE

#ifdef CGAL_USE_GMP
/*! \ingroup CGAL_Arithmetic_kernel
 *  \brief  The GMP set of exact number types
 */
class GMP_arithmetic_kernel {
public:
    //! exact integers
    typedef CGAL::Gmpz Integer;
    //! exact float nummber
    typedef CGAL::Gmpq Exact_float_number;
    //! exact rationals, constructible from integers
    typedef CGAL::Gmpq Rational;
    //! exact root expressions, constructible from integers and rationals
    typedef CGAL::Null_tag  Field_with_sqrt;
};
#endif // CGAL_USE_GMP

#ifdef CGAL_USE_GMPXX
/*! \ingroup CGAL_Arithmetic_kernel
 *  \brief  The GMP set of exact number types
 */
class GMPXX_arithmetic_kernel {
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
//   typedef typename AT::Poly_int1 Poly_int1; 
//   typedef typename AT::Poly_int2 Poly_int2; 
//   typedef typename AT::Poly_int3 Poly_int3;
//   typedef typename AT::Poly_rat1 Poly_rat1; 
//   typedef typename AT::Poly_rat2 Poly_rat2; 
//   typedef typename AT::Poly_rat3 Poly_rat3;

// end #define

template <class NT>
class Lazy_exact_nt;


namespace INTERN_AK{
template <class NT>
struct Lazy_exact_type{
    typedef CGAL::Lazy_exact_nt<NT> type;
};
template <>
struct Lazy_exact_type<CGAL::Null_tag>{
    typedef CGAL::Null_tag type;
};
} // namespace INTERN_AK

template <class AT>
struct Lazy_exact_arithmetic_kernel{
    typedef typename INTERN_AK::Lazy_exact_type<typename AT::Integer>::type Integer;
    typedef typename INTERN_AK::Lazy_exact_type<typename AT::Exact_float_number>::type Exact_float_number;
    typedef typename INTERN_AK::Lazy_exact_type<typename AT::Rational>::type Rational;
    typedef typename INTERN_AK::Lazy_exact_type<typename AT::Field_with_sqrt>::type Field_with_sqrt;
};


#ifdef CGAL_INTERN_USE_BFI

namespace CGALi {
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
#endif //CGAL_USE_CORE
    
    template <class COEFF, class ROOT>
    struct Get_arithmetic_kernel<Sqrt_extension<COEFF,ROOT> >{
        typedef Get_arithmetic_kernel<COEFF> GET;
        typedef typename GET::Arithmetic_kernel Arithmetic_kernel;
    };
    
    
} // namespace CGALi
#endif //CGAL_INTERN_USE_BFI


CGAL_END_NAMESPACE

#endif // CGAL_ARITHMETIC_KERNEL_H
// EOF
