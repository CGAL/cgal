// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: 4515 rschindl$
// 
//
// Author(s)     : Ralf Schindlbeck <rschindl@mpi-inf.mpg.de>
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.



/*! \file CGAL/core_interval_support.h
  This is experimental 
  use CORE::BigFloat as interval type. 
*/


#ifndef CGAL_CORE_INTERVAL_SUPPORT_H
#define CGAL_CORE_INTERVAL_SUPPORT_H

#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE
#warning This header file needs CORE installed in order to work properly.
#else // CGAL_USE_CORE

//#include <NiX/CORE/BigFloat.h>
//#include <NiX/CORE/Expr.h>    
//#include <NiX/CORE/BigInt.h>
//#include <NiX/CORE/BigRat.h>

#include <CGAL/CORE/BigFloat.h>

// namespace NiX{
CGAL_BEGIN_NAMESPACE

template<typename BFI> long get_significant_bits(BFI bfi);

   CORE::BigFloat 
   inline 
   round(const CORE::BigFloat& x, long rel_prec = CORE::defRelPrec.toLong() ){
   CGAL_postcondition(rel_prec >= 0);   
   // since there is not rel prec defined if in_zero(x)
   if (x.isZeroIn()) return x; 
   if (CGAL::get_significant_bits(x) <= rel_prec) return x;
   
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

    CGAL_postcondition(CGAL::abs(CGAL::get_significant_bits(xr) - rel_prec) <= 1);   
    CGAL_postcondition(BF(xr.m()-xr.err(),0,xr.exp()) <= BF(x.m()-x.err(),0,x.exp()));
    CGAL_postcondition(BF(xr.m()+xr.err(),0,xr.exp()) >= BF(x.m()+x.err(),0,x.exp()));
    return xr;     
}

template<typename BFI> class  Bigfloat_interval_traits;

template<> class Bigfloat_interval_traits<CORE::BigFloat> 

   {
   public:
     typedef CORE::BigFloat NT;
     typedef CORE::BigFloat BF;

     typedef Bigfloat_interval_traits<NT> Self;

     class Get_significant_bits {
     public:
         // type for the \c AdaptableUnaryFunction concept.
         typedef NT  argument_type;
         // type for the \c AdaptableUnaryFunction concept.
         typedef long  result_type;

         long operator()( NT x) const {
           //std::cout << ":::::::::::::::::::: "<< std::endl;
           if(x.err() == 0 ) {
             //std::cout << "Exact: \n" << x.m() << std::endl;
                 return ::CORE::bitLength(x.m()); 
           }
           else {
             //std::cout << "With error:\n" << x.m() << std::endl << x.err() << std::endl;
             //std::cout << "bitlength m: " << ::CORE::bitLength(x.m()) << "\nbitlength e: " << ::CORE::bitLength(x.err()) << std::endl;
             return ::CORE::bitLength(x.m()) - ::CORE::bitLength(x.err());
           }

         }
     };

     class Set_precision {
     public:
         // type for the \c AdaptableUnaryFunction concept.
         typedef long  argument_type;
         // type for the \c AdaptableUnaryFunction concept.
         typedef long  result_type;  
     
         long operator() ( long prec ) const {
             long result =  ::CORE::defRelPrec.toLong();
             ::CORE::defRelPrec = prec; 
             ::CORE::defBFdivRelPrec = prec;
             return result; 
         }
     };
     
     class Get_precision {
     public:
         // type for the \c AdaptableGenerator concept.
         typedef long  result_type;  
     
         long operator() () const {
             return  ::CORE::defRelPrec.toLong(); 
         }
     };

     class Upper {
     public:
         // type for the \c AdaptableUnaryFunction concept.
         typedef NT  argument_type;
         // type for the \c AdaptableUnaryFunction concept.
         typedef BF  result_type;
         
         BF operator() ( NT x ) const {
             CORE::BigFloat result = ::CORE::BigFloat(x.m()+x.err(),0,x.exp());
             CGAL_postcondition(result >= x);
             return result; 
         }
     };

     class Lower {
     public:
         // type for the \c AdaptableUnaryFunction concept.
         typedef NT  argument_type;
         // type for the \c AdaptableUnaryFunction concept.
         typedef BF  result_type;
         
         BF operator() ( NT x ) const {
             CORE::BigFloat result = ::CORE::BigFloat(x.m()-x.err(),0,x.exp());
             CGAL_postcondition(result <= x);
             return result; 
         }
     };

     
     class In_zero {
     public:
         // type for the \c AdaptableUnaryFunction concept.
         typedef NT  argument_type;
         // type for the \c AdaptableUnaryFunction concept.
         typedef bool  result_type;
         
         bool operator() ( NT x ) const {
             return x.isZeroIn(); 
         }
     };

     class Overlap {
     public:
         // type for the \c AdaptableBinaryFunction concept.
         typedef NT  first_argument_type;
         // type for the \c AdaptableBinaryFunction concept.
         typedef NT  second_argument_type;
         // type for the \c AdaptableBinaryFunction concept.
         typedef bool  result_type;
    
         bool operator() ( NT x, NT y ) const {
             Self::In_zero in_zero;
             bool result = in_zero(x-y);
             return result;
         }
     };
   
     class Hull {
     public:
         // type for the \c AdaptableBinaryFunction concept.
         typedef NT  first_argument_type;
         // type for the \c AdaptableBinaryFunction concept.
         typedef NT  second_argument_type;
         // type for the \c AdaptableBinaryFunction concept.
         typedef NT  result_type;
    
         NT operator() ( NT x, NT y ) const {
             CORE::BigFloat result;
             
             // Unfortunately, CORE::centerize(x,y) has bugs. 
             if ((x.m() == y.m()) && (x.err() == y.err()) && 
                 (x.exp() == y.exp())) {
                 
                 //rep(x) == rep(y) 
                 return x;
             }
             Self::Lower lower_functor;
             Self::Upper upper_functor;
             CORE::BigFloat lower = std::min(lower_functor(x),
                                             lower_functor(y));
             CORE::BigFloat upper = std::max(upper_functor(x),
                                             upper_functor(y));
             CORE::BigFloat mid = (lower + upper)/2;
             CORE::BigFloat err = CGAL::round((upper - lower)/2,0);
             err = CORE::BigFloat(0,err.m().longValue()+err.err(),err.exp());
             result = mid + err;  
             
             CGAL_postcondition(lower_functor(result) <= 
                          std::min(lower_functor(x),
                                   lower_functor(y)));
             CGAL_postcondition(upper_functor(result) >= 
                          std::max(upper_functor(x),
                                   upper_functor(y)));
             return result ;
         }
     };

     class Singleton {
     public:
         // type for the \c AdaptableUnaryFunction concept.
         typedef NT  argument_type;
         // type for the \c AdaptableUnaryFunction concept.
         typedef bool  result_type;
         
         bool operator() ( NT x ) const {
             return (x.err() == 0); 
         }
     };

     class Width {
     public:
         // type for the \c AdaptableUnaryFunction concept.
         typedef NT  argument_type;
         // type for the \c AdaptableUnaryFunction concept.
         typedef BF  result_type;
         
         BF operator() ( NT x ) const {
           unsigned long err = 2*x.err();
           return BF(CORE::BigInt(err),0,x.exp());
         }
     };

     class Convert_to_bfi {
     public:

         typedef NT result_type;

             
         NT operator() (const ::CORE::Expr& x){
             return round(CORE::BigFloat(x));
         }

         NT operator() (const ::CORE::BigInt& x){
             return round(CORE::BigFloat(x));
         }
         
         NT operator() (const ::CORE::BigRat& x){
             return round(CORE::BigFloat(x));
         }
     };


   };



// } // namespace NiX
CGAL_END_NAMESPACE

#endif // CGAL_USE_CORE
#endif // CGAL_CORE_INTERVAL_SUPPORT_H
