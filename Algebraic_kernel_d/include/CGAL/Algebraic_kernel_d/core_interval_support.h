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

/*! \file NiX/CORE/interval_support.h
  This is experimental 
  use CORE::BigFloat as interval type. 
*/


#ifndef CGAL_ALGEBRAIC_KERNEL_D_CORE_INTERVAL_SUPPORT_H
#define CGAL_ALGEBRAIC_KERNEL_D_CORE_INTERVAL_SUPPORT_H

#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE
#warning This header file needs CORE installed in order to work properly.
#else // LiS_HAVE_CORE

#include <CGAL/CORE_BigFloat.h>

/*#include <NiX/CORE/BigFloat.h>
#include <NiX/CORE/Expr.h>    
#include <NiX/CORE/BigInt.h>
#include <NiX/CORE/BigRat.h>
#include <NiX/Algebraic_real.h>*/

CGAL_BEGIN_NAMESPACE

namespace CGALi {

CORE::BigFloat
inline
round(const CORE::BigFloat& x, long rel_prec = CORE::defRelPrec.toLong() ){
    // since there is not rel prec defined if in_zero(x)
    if (x.isZeroIn()) return x;

    typedef CORE::BigFloat BF;
    typedef CORE::BigFloat BFI;
    typedef CORE::BigInt Integer;
    BF xr;

    CORE::BigInt m = x.m();
    long         err = x.err();
    long         exp = x.exp();

    long bits  = ::CORE::bitLength(m);
    long shift = bits - rel_prec;

    if( shift > 0 ){
        Integer new_m   = m >> shift ;
        long    new_err = (err >> shift)+2;
        //NiX_postcond( (x.m()<0) == (new_m<0));
        xr = BF(new_m,new_err,0)*BF::exp2(exp*14+shift);
    }else{
        // noting to do
        xr = x;
    }

    CGAL_postcondition(BF(xr.m()-xr.err(),0,xr.exp()) <= BF(x.m()-x.err(),0,x.exp()));
    CGAL_postcondition(BF(xr.m()+xr.err(),0,xr.exp()) >= BF(x.m()+x.err(),0,x.exp()));

    return xr;
}



long  
inline 
get_significant_bits(const CORE::BigFloat& x){
    //if (x.isExact()) return 1000000 ; // CORE::CORE_posInfty.toLong(); 
    if(x.err() == 0 ) 
        return ::CORE::bitLength(x.m()); 
    else
        return ::CORE::bitLength(x.m()) - ::CORE::bitLength(x.err());
}

long 
inline 
set_precision(CORE::BigFloat, long prec){
    long result =  ::CORE::defRelPrec.toLong();
    ::CORE::defRelPrec = prec; 
    return result; 
}
long 
inline
get_precision(CORE::BigFloat){
    return  ::CORE::defRelPrec.toLong(); 
}

CORE::BigFloat 
inline
upper(CORE::BigFloat x){
    return ::CORE::BigFloat(x.m()+x.err(),0,x.exp());
}

CORE::BigFloat 
inline 
lower(CORE::BigFloat x){
    return ::CORE::BigFloat(x.m()-x.err(),0,x.exp());
}

bool 
inline 
in_zero(CORE::BigFloat x){
    return x.isZeroIn(); 
}

bool 
inline 
overlap(CORE::BigFloat x, CORE::BigFloat y){
    bool result = in_zero(x-y);
    return result;
}

CORE::BigFloat
inline 
hull(CORE::BigFloat x, CORE::BigFloat y){
    return CORE::centerize(x,y);
}
bool
inline
singleton(CORE::BigFloat x){
    return x.isExact();
}


CORE::BigFloat 
inline
convert_to_bfi( ::CORE::Expr x){
    return CGALi::round(CORE::BigFloat(x));
}

CORE::BigFloat
inline
convert_to_bfi(const ::CORE::BigInt& x){
    return CGALi::round(CORE::BigFloat(x));
}

CORE::BigFloat 
inline
convert_to_bfi(const ::CORE::BigRat& x){
    return CGALi::round(CORE::BigFloat(x));
}


std::pair<double, double> to_interval( const CORE::BigFloat& x ) {
    // TODO: is the precission correct when using this?
    return CGAL::to_interval( x );
}

} // namespace CGALi 

CGAL_END_NAMESPACE

#endif // CGAL_USE_CORE
#endif // CGAL_ALGEBRAIC_KERNEL_D_CORE_INTERVAL_SUPPORT_H
