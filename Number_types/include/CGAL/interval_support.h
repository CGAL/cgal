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


/*! \file CGAL/interval_support.h
  This is experimental 
*/

#ifndef CGAL_INTERVAL_SUPPORT_H
#define CGAL_INTERVAL_SUPPORT_H

#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel.h>

// #include <NiX/Float_traits.h>
#include <CGAL/Float_traits.h>

///////// covnersion tools: 
// namespace NiX{
CGAL_BEGIN_NAMESPACE

template <class NT>  class Get_arithmetic_traits;

template <typename A,typename B,typename C,typename D,typename E> 
class Algebraic_real;

// } // namespace NiX
// namespace CGAL {
template <typename A,typename B> class Sqrt_extension;
// } // namespace CGAL
// namespace NiX {
// Bigfloat_interval_functions

template<typename BigfloatInterval> class Bigfloat_interval_traits;

template<typename BFI> long get_significant_bits(BFI bfi) {
    typename Bigfloat_interval_traits<BFI>::Get_significant_bits 
      get_significant_bits;
    return get_significant_bits(bfi);
}

/*
 template<typename BFI> long set_precision(BFI bfi,long prec) {
   typename Float_traits<typename Bigfloat_interval_traits<BFI>::BF>
     ::Set_precision set_precision;
    return set_precision(prec);
}

template<typename BFI> long get_precision(BFI bfi) {
    typename Float_traits<typename Bigfloat_interval_traits<BFI>::BF>
      ::Get_precision get_precision;
    return get_precision();
}
*/

// ONLY FOR TESTING THIS IS THE WRONG FILE!!!
 template<typename BF> long set_precision(BF bfi,long prec) {
   typename Float_traits<BF>::Set_precision set_precision;
    return set_precision(prec);
}

template<typename BF> long get_precision(BF bfi) {
    typename Float_traits<BF>::Get_precision get_precision;
    return get_precision();
}


template<typename BFI> 
  typename Bigfloat_interval_traits<BFI>::BF upper(BFI bfi) {
    typename Bigfloat_interval_traits<BFI>::Upper upper;
    return upper(bfi);
}

template<typename BFI> 
  typename Bigfloat_interval_traits<BFI>::BF lower(BFI bfi) {
    typename Bigfloat_interval_traits<BFI>::Lower lower;
    return lower(bfi);
}

template<typename BFI> BFI hull(BFI bfi1, BFI bfi2) {
   typename Bigfloat_interval_traits<BFI>::Hull hull;
   return hull(bfi1, bfi2);
}

template<typename BFI> bool in_zero(BFI bfi) {
   typename Bigfloat_interval_traits<BFI>::In_zero in_zero;
   return in_zero(bfi);
}

template<typename BFI> bool overlap(BFI bfi1, BFI bfi2) {
   typename Bigfloat_interval_traits<BFI>::Overlap overlap;
   return overlap(bfi1, bfi2);
}

template<typename BFI> typename Bigfloat_interval_traits<BFI>::BF
   width(BFI bfi) {

   typename Bigfloat_interval_traits<BFI>::Width width;
   return width(bfi);
}

template<typename BFI> bool singleton(BFI bfi) {
   typename Bigfloat_interval_traits<BFI>::Singleton singleton;
   return singleton(bfi);
}

template<typename NT> typename CGALi::Get_arithmetic_kernel<NT>::Bigfloat_interval   
    convert_to_bfi(NT x) {

    typedef typename 
      CGALi::Get_arithmetic_kernel<NT>::Bigfloat_interval BFI;
    typename Bigfloat_interval_traits<BFI>::CGALi::Get_arithmetic_kernel 
      get_arithmetic_traits;
    return get_arithmetic_traits(x);
}

#if 0
//TODO
//porting Specialisation for double-intervals from EXACUS2CGAL

// (Parital) Specialisation for double-intervals

#include <NiX/Interval.h>

template<>
  class Bigfloat_interval_traits<Interval> 
  {
   public:
    typedef Interval NT;

     typedef double BF;

     /* ??
     class Get_significant_bits {
     public:
         // type for the \c AdaptableUnaryFunction concept.
         typedef NT  argument_type;
         // type for the \c AdaptableUnaryFunction concept.
         typedef long  result_type;

         long operator()( NT x) const {}
     };
     */
     
     class Upper {
     public:
         // type for the \c AdaptableUnaryFunction concept.
         typedef NT  argument_type;
         // type for the \c AdaptableUnaryFunction concept.
         typedef BF  result_type;
         
         BF operator() ( NT a ) const {
             return a.upper();
         }
     };

     class Lower {
     public:
         // type for the \c AdaptableUnaryFunction concept.
         typedef NT  argument_type;
         // type for the \c AdaptableUnaryFunction concept.
         typedef BF  result_type;
         
         BF operator() ( NT a ) const {
             return a.lower();
         }
     };

     class In_zero {
     public:
         // type for the \c AdaptableUnaryFunction concept.
         typedef NT  argument_type;
         // type for the \c AdaptableUnaryFunction concept.
         typedef bool  result_type;
         
         bool operator() ( NT x ) const {
             return ::boost::numeric::in_zero(x);
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
             return ::boost::numeric::overlap(x,y);
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
             return ::boost::numeric::hull(x,y);
         }
     };

     class Singleton {
     public:
         // type for the \c AdaptableUnaryFunction concept.
         typedef NT  argument_type;
         // type for the \c AdaptableUnaryFunction concept.
         typedef bool  result_type;
         
         bool operator() ( NT a ) const {
           return ::boost::numeric::singleton(a); 
         }
     };

     class Width {
     public:
         // type for the \c AdaptableUnaryFunction concept.
         typedef NT  argument_type;
         // type for the \c AdaptableUnaryFunction concept.
         typedef BF  result_type;
         
         BF operator() ( NT a ) const {
             return ::boost::numeric::width(a); 
         }

     };

     class Convert_to_bfi {
     public:

         typedef NT result_type;

         template<typename NTX>
           NT operator() ( NTX x) {
           return to_interval( x );
         }
         
     };

   };
#endif


#if 0
template <class COEFF, class HP >
leda_bigfloat_interval  
inline 
convert_to_bfi(const Algebraic_real< COEFF, ::leda::real, ::leda::rational, HP >& x) {
     
    //std::cout << "convert_to_bfi(const Algebraic_real<COEFF, FWS, RAT, HP >& x)"<<std::endl;
    typedef ::leda::real FWS;
    typedef ::leda::rational RAT; 
    typedef Algebraic_real< COEFF, FWS, RAT, HP > ALG;
    
    // if zero done
    if(NiX::sign(x) == CGAL::ZERO)
        return (leda_bigfloat_interval(0));    
    CGAL_postcondition(NiX::sign(x.low()) == NiX::sign(x.high()));
  
    
    typedef ::leda::bigfloat BF;
    long final_prec = BF::set_precision(BF::get_precision()+4);
  
    ::leda::bigfloat rerror(2,-final_prec);

    leda_bigfloat_interval bfi(convert_to_bfi(x.low()).lower(), convert_to_bfi(x.high()).upper());
    
    //while( NiX::width(bfi) >=  rerror){ 
    //while( NiX::width(bfi) >=  NiX::median(NiX::abs(bfi))*rerror){ 
    //std::cout << "rel error:  "<< (NiX::width(bfi) / NiX::abs(bfi)).upper() <<std::endl;
    while( (NiX::width(bfi) / NiX::abs(bfi)).upper() > rerror ){
        
        //std::cout <<" diff  " << NiX::width(bfi) - NiX::abs(bfi).upper()*rerror  << std::endl;
        x.refine();
        bfi = leda_bigfloat_interval(convert_to_bfi(x.low()).lower(), convert_to_bfi(x.high()).upper());   
        //std::cout <<"bfi: "<< bfi << std::endl;
        //std::cout << "final_prec: "<< final_prec << std::endl;
        //std::cout << "rerror    : "<< rerror << std::endl;
        //std::cout << "rel error:  "<< (NiX::width(bfi) / NiX::abs(bfi)).upper() <<std::endl;
    }
    //std::cout << "convert_to_bfi(const Algebraic_real<COEFF, FWS, RAT, HP >& x) end"<<std::endl;
    BF::set_precision(final_prec);

#ifndef NDEBUG  
    ::leda::rational lower_num(bfi.lower().get_significant());
    ::leda::rational lower_denom(NiX::ipower(RAT(2),NiX::abs(bfi.lower().get_exponent().to_long())));
    ::leda::rational lower;
    if(bfi.lower().get_exponent().to_long() < 0 )
        lower = lower_num / lower_denom ;
    else
        lower = lower_num  * lower_denom;
    CGAL_postcondition( x.compare(lower) == CGAL::LARGER );

    ::leda::rational upper_num(bfi.upper().get_significant());
    ::leda::rational upper_denom(NiX::ipower(RAT(2),NiX::abs(bfi.upper().get_exponent().to_long())));
    ::leda::rational upper;
    if(bfi.upper().get_exponent().to_long() < 0 )
        upper = upper_num / upper_denom ;
    else
        upper = upper_num  * upper_denom;
    CGAL_postcondition( x.compare(upper) == CGAL::SMALLER );
#endif 

    return bfi; 
}
#endif

template <class NTX>
typename CGALi::Get_arithmetic_kernel<NTX>::Arithmetic_kernel::Bigfloat_interval
convert_to_bfi(const NTX& x) {
    typedef typename CGALi::Get_arithmetic_kernel<NTX>::Arithmetic_kernel AT;
    typedef typename AT::Bigfloat_interval BFI;
    typename Bigfloat_interval_traits<BFI>::Convert_to_bfi convert_to_bfi;
    return convert_to_bfi(x);
}

template <typename NT, typename ROOT>
typename CGALi::Get_arithmetic_kernel<NT>::Arithmetic_kernel::Bigfloat_interval
convert_to_bfi(const CGAL::Sqrt_extension<NT,ROOT>& x) {
    typedef typename CGALi::Get_arithmetic_kernel<NT>::Arithmetic_kernel AT;
    typedef typename AT::Bigfloat_interval BFI;
    if(x.is_extended()){
        BFI a0(convert_to_bfi(x.a0()));
        BFI a1(convert_to_bfi(x.a1()));
        BFI root(CGAL::sqrt(convert_to_bfi(x.root())));
        return a0+a1*root;
    }else{
        return convert_to_bfi(x.a0());
    }
}

template <class COEFF, class FWS, class RAT, class T1, class T2>
typename CGALi::Get_arithmetic_kernel<COEFF>::Arithmetic_kernel::Bigfloat_interval
inline
convert_to_bfi(const Algebraic_real< COEFF, FWS, RAT,T1,T2>& x) {
    typedef typename CGALi::Get_arithmetic_kernel<COEFF>::Arithmetic_kernel AT;
    typedef typename AT::Bigfloat_interval BFI;
    typedef typename AT::Bigfloat          BF;
    BFI result = x.convert_to_bfi();
    
#ifndef NDEBUG
    CGAL::set_precision(BF(),CGAL::get_precision(BF())*2);
    CGAL_postcondition(CGAL::lower(result) <= CGAL::lower(CGAL::convert_to_bfi(x.low() )));
    CGAL_postcondition(CGAL::upper(result) >= CGAL::upper(CGAL::convert_to_bfi(x.high())));
    CGAL::set_precision(BF(),CGAL::get_precision(BF())/2);
#endif 
    return result; 
};  

// } // namespace NiX
CGAL_END_NAMESPACE

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_interval_support.h>
#endif //  CGAL_USE_LEDA

#ifdef CGAL_USE_CORE
#include <CGAL/core_interval_support.h>
#endif //  CGAL_USE_CORE

#endif // CGAL_INTERVAL_SUPPORT_H
