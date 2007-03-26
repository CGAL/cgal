// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     :  
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.


/*! \file NiX/interval_support.h
  This is experimental 
*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_INTERVAL_SUPPORT_H
#define CGAL_ALGEBRAIC_KERNEL_D_INTERVAL_SUPPORT_H

#include <CGAL/basic.h> 
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Sqrt_extension.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/Algebraic_kernel_d/leda_interval_support.h>
#endif //  LiS_HAVE_LEDA

#ifdef CGAL_USE_CORE
#include <CGAL/Algebraic_kernel_d/core_interval_support.h>
#endif //  LiS_HAVE_LEDA


///////// covnersion tools: 
CGAL_BEGIN_NAMESPACE

namespace CGALi {

// Forward declaration of Algebraic_real_pure
template< class COEFF, class RAT, class POLICY >
class Algebraic_real_pure;

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


template <class NTX>
typename Get_arithmetic_kernel<NTX>::Arithmetic_kernel::Bigfloat_interval
convert_to_bfi(const NTX& x) {
    typedef typename Get_arithmetic_kernel<NTX>::Arithmetic_kernel AT;
    typedef typename AT::Field_with_sqrt FWS;
    typename CGAL::Coercion_traits< FWS , NTX >::Cast cast;
    return convert_to_bfi(cast(x));
}


template <class COEFF, class RAT, class POLICY >
typename Get_arithmetic_kernel<COEFF>::Arithmetic_kernel::Bigfloat_interval
inline
convert_to_bfi(const Algebraic_real_pure< COEFF, RAT, POLICY >& x){
    typedef typename Get_arithmetic_kernel<COEFF>::Arithmetic_kernel AT;
    typedef typename AT::Bigfloat BF;
    typedef typename AT::Bigfloat_interval BFI; 
    typedef Algebraic_real_pure< COEFF, RAT, POLICY > ALG;

    if (x.is_rational()) return convert_to_bfi(x.rational());
    if(CGAL::sign(x) == CGAL::ZERO) return (BF(0));
    
    CGAL_postcondition(CGAL::sign(x.low()) == CGAL::sign(x.high()));
    long final_prec = set_precision(BF(),get_precision(BF())+4);
  
    BFI bfi = CGALi::hull(convert_to_bfi(x.low()), convert_to_bfi(x.high()));
    
    while( !singleton(bfi) &&  get_significant_bits(bfi) < final_prec  ){
        x.refine();
        bfi = CGALi::hull(
                convert_to_bfi(x.low()), 
                convert_to_bfi(x.high()));
    }

    set_precision(BF(),final_prec);
    return bfi; 
}


} // namespace CGALi

CGAL_END_NAMESPACE

#endif // NiX_INTERVAL_SUPPORT_H
