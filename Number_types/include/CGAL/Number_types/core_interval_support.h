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


#ifndef CGAL_NUMBER_TYPES_CORE_INTERVAL_SUPPORT_H
#define CGAL_NUMBER_TYPES_CORE_INTERVAL_SUPPORT_H

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

// Forward declarations
long inline get_significant_bits(const CORE::BigFloat& x);

CORE::BigFloat 
inline 
round(const CORE::BigFloat& x, long rel_prec = CORE::defRelPrec.toLong() ){
    CGAL_postcondition(rel_prec >= 0);   
    // since there is not rel prec defined if in_zero(x)
    if (x.isZeroIn()) return x; 
    
    if (CGALi::get_significant_bits(x) <= rel_prec) return x;
    
    typedef CORE::BigFloat BF; 
    typedef CORE::BigFloat BFI; 
    typedef CORE::BigInt Integer;
    BF xr;
    
    CORE::BigInt m = x.m();
    long         err = x.err();
    long         exp = x.exp(); 
    
    long shift = ::CORE::bitLength(m) - rel_prec - 1;
    if( shift > 0 ){
        Integer new_m   = m >> shift ; 
        if(err == 0){
            xr = BF(new_m,1,0)*BF::exp2(exp*14+shift);
        }else{
            xr = BF(new_m,2,0)*BF::exp2(exp*14+shift);
        }
    }else{
        // noting to do
        xr = x; 
    }

    CGAL_postcondition(CGAL::abs(CGALi::get_significant_bits(xr) - rel_prec) <= 1);   
    CGAL_postcondition(BF(xr.m()-xr.err(),0,xr.exp()) <= BF(x.m()-x.err(),0,x.exp()));
    CGAL_postcondition(BF(xr.m()+xr.err(),0,xr.exp()) >= BF(x.m()+x.err(),0,x.exp()));
    return xr;     
}

long  
inline 
get_significant_bits(const CORE::BigFloat& x){
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
    CORE::BigFloat result = ::CORE::BigFloat(x.m()+x.err(),0,x.exp());
    CGAL_postcondition(result >= x);
    return result; 
}

CORE::BigFloat 
inline 
lower(CORE::BigFloat x){
    CORE::BigFloat result = ::CORE::BigFloat(x.m()-x.err(),0,x.exp());
    CGAL_postcondition(result <= x);
    return result; 
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
    //std::cout <<" x:" << x.m()<<"+-"<< x.err()<<"*2^14*"<<x.exp()<<std::endl;
    //std::cout <<" y:" << y.m()<<"+-"<< y.err()<<"*2^14*"<<y.exp()<<std::endl;
    CORE::BigFloat result;
    
    // Unfortunaltely, CORE::centerize(x,y) has bugs. 
    if ((x.m() == y.m()) && (x.err() == y.err()) && (x.exp() == y.exp())) {
        //rep(x) == rep(y) 
        return x;
    }
    CORE::BigFloat lower = std::min(CGALi::lower(x),CGALi::lower(y));
    CORE::BigFloat upper = std::max(CGALi::upper(x),CGALi::upper(y));
    CORE::BigFloat mid = (lower + upper)/2;
    CORE::BigFloat err = CGALi::round((upper - lower)/2,0);
    err = CORE::BigFloat(0,err.m().longValue()+err.err(),err.exp());
    result = mid + err;  

    CGAL_postcondition(CGALi::lower(result) <= std::min(CGALi::lower(x),CGALi::lower(y)));
    CGAL_postcondition(CGALi::upper(result) >= std::max(CGALi::upper(x),CGALi::upper(y)));
    return result ;
}
bool
inline
singleton(CORE::BigFloat x){
    return (x.err() == 0);
}

CORE::BigFloat 
inline
convert_to_bfi(const ::CORE::Expr& x){
    return round(CORE::BigFloat(x));
}

CORE::BigFloat
inline
convert_to_bfi(const ::CORE::BigInt& x){
    return round(CORE::BigFloat(x));
}

CORE::BigFloat 
inline
convert_to_bfi(const ::CORE::BigRat& x){
    return round(CORE::BigFloat(x));
}

// We need this function, because we need to_interval in namespace CGAL
std::pair<double, double> 
inline to_interval( const CORE::BigFloat& x ) {
    // TODO: is the precission correct when using this?
    return CGAL::to_interval( x );
}

CORE::BigFloat
inline sqrt( const CORE::BigFloat& x ) {
    return CGAL::sqrt( x );
}


/*CORE::BigFloat
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
        //CGAL_postcondition( (x.m()<0) == (new_m<0));
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
}*/


} // namespace CGALi 

CGAL_END_NAMESPACE

#endif // CGAL_USE_CORE
#endif // CGAL_NUMBER_TYPES_CORE_INTERVAL_SUPPORT_H
