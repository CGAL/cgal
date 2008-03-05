// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id:$
// 
//
// Author(s)     : Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

// TODO: 
// - mv general functors in a new Interval_traits 
// - mv function 'round' to Interval_traits

#ifndef CGAL_CORE_INTERVAL_SUPPORT_H
#define CGAL_CORE_INTERVAL_SUPPORT_H

#include <CGAL/basic.h>
#include <CGAL/interval_support.h>

#ifndef CGAL_USE_CORE
#warning This header file needs CORE installed in order to work properly.
#else // CGAL_USE_CORE

#include <CGAL/CORE/BigFloat.h>

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
            if(x.err() == 0 ) {
                return ::CORE::bitLength(x.m()); 
            }
            else {
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
#if 0
            // this is not possible since CORE::centerize has a bug.
            NT result = CORE::centerize(x,y);
#else 

            CORE::BigFloat result;
             
            // Unfortunately, CORE::centerize(x,y) has bugs. 
            if ((x.m() == y.m()) && (x.err() == y.err()) && (x.exp() == y.exp())) {
                return x;
            }
                         
            CORE::BigFloat lower = std::min(CGAL::lower(x),
                    CGAL::lower(y));
            CORE::BigFloat upper = std::max(CGAL::upper(x),
                    CGAL::upper(y));

            CORE::BigFloat mid = (lower + upper)/2;
             
            // Now we have to compute the error. The problem is that .err() is just a long
            CORE::BigFloat err = (upper - lower)/CORE::BigFloat(2);
                   
            //std::cout << "lower    " << lower << std::endl;
            //std::cout << "upper    " << upper << std::endl;
            //std::cout << "mid      " << mid << std::endl;
            //std::cout << "err I    " << err << std::endl;
            
            // shift such that err.m()+err.err() fits into long 
            int digits_long = std::numeric_limits<long>::digits;
            if(::CORE::bitLength(err.m()) >= digits_long){ 
                long shift = ::CORE::bitLength(err.m()) - digits_long + 1 ; 
                //std::cout << "shift " << shift<< std::endl;
                long new_err = (err.m()+err.err() >> shift).longValue()+1; 
                err = CORE::BigFloat(0,new_err,0) * CORE::BigFloat::exp2(err.exp()*14+shift);
            }else{
                err = CORE::BigFloat(0,err.m().longValue()+err.err(),err.exp());
            }

            result = mid + err;  
             
#endif 
            CGAL_postcondition( 
                    CGAL::lower(result) 
                    <=  std::min(CGAL::lower(x), CGAL::lower(y)));
            CGAL_postcondition( 
                    CGAL::upper(result) 
                    >= std::max(CGAL::upper(x), CGAL::upper(y)));

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

CGAL_END_NAMESPACE

#endif // CGAL_USE_CORE
#endif // CGAL_CORE_INTERVAL_SUPPORT_H
