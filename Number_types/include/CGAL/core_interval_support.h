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

template<> 
class Interval_traits<CORE::BigFloat> 
    : public CGALi::Interval_traits_base<CORE::BigFloat>{
public: 
    
    typedef Interval_traits<CORE::BigFloat> Self; 
    typedef CORE::BigFloat Interval;
    typedef CORE::BigFloat Boundary;
    typedef CGAL::Tag_true Is_interval; 
 
    struct Lower :public Unary_function<Interval,Boundary>{
        Boundary operator() ( Interval x ) const {   
            CORE::BigFloat result = ::CORE::BigFloat(x.m()-x.err(),0,x.exp());
            CGAL_postcondition(result <= x);
            return result; 
        }
    };
    
    struct Upper :public Unary_function<Interval,Boundary>{
        Boundary operator() ( Interval x ) const {     
            CORE::BigFloat result = ::CORE::BigFloat(x.m()+x.err(),0,x.exp());
            CGAL_postcondition(result >= x);
            return result; 
        }
    };

    struct Width :public Unary_function<Interval,Boundary>{
         
        Boundary operator() ( Interval x ) const {    
            unsigned long err = 2*x.err();
            return Boundary(CORE::BigInt(err),0,x.exp());
        }
    };

    struct Median :public Unary_function<Interval,Boundary>{
         
        Boundary operator() ( Interval x ) const {   
            return Boundary(x.m(),0,x.exp());
        }
    };

    struct Norm :public Unary_function<Interval,Boundary>{
        Boundary operator() ( Interval x ) const {      
            return std::max(Upper()(x).abs(),Lower()(x).abs());
        }
    };
    
    struct Zero_in :public Unary_function<Interval,bool>{
        bool operator() ( Interval x ) const {      
            return x.isZeroIn(); 
        }
    };

    struct In :public Binary_function<Boundary,Interval,bool>{  
        bool operator()( Boundary x, const Interval& a ) const {    
            CGAL_precondition(CGAL::singleton(x));
            return (Lower()(a) <= x && x <= Upper()(a));
        }
    };

    struct Equal :public Binary_function<Interval,Interval,bool>{  
        bool operator()( const Interval& a, const Interval& b ) const { 
            return (Upper()(a) == Upper()(b) &&  Lower()(a) == Lower()(b));
        }
    };
    
    struct Subset :public Binary_function<Interval,Interval,bool>{  
        bool operator()( const Interval& a, const Interval& b ) const {   
            return Lower()(b) <= Lower()(a) && Upper()(a) <= Upper()(b);
        }
    };
    
    struct Proper_subset :public Binary_function<Interval,Interval,bool>{ 
        bool operator()( const Interval& a, const Interval& b ) const { 
            return Subset()(a,b) && (!Equal()(a,b));
        }
    };
    
    struct Intersection :public Binary_function<Interval,Interval,Interval>{ 
        Interval operator()( const Interval& a, const Interval& b ) const {
            // std::cout <<"a= (" << a.m() << "+-" << a.err() << ")*2^" << a.exp() << std::endl;
            Boundary l(CGAL::max(Lower()(a),Lower()(b)));
            Boundary u(CGAL::min(Upper()(a),Upper()(b)));

            if(u < l ) throw Exception_intersection_is_empty();
            return Construct()(l,u);
        }
    };
 

    struct Overlap :public Binary_function<Interval,Interval,bool>{
        bool operator() ( Interval x, Interval y ) const {       
            Self::Zero_in Zero_in;
            bool result = Zero_in(x-y);
            return result;
        }
    };
   
    struct Hull :public Binary_function<Interval,Interval,Interval>{
        Interval operator() ( Interval x, Interval y ) const {
#if 0
            // this is not possible since CORE::centerize has a bug.
            Interval result = CORE::centerize(x,y);
#else 

            CORE::BigFloat result;
             
            // Unfortunately, CORE::centerize(x,y) has bugs. 
            if ((x.m() == y.m()) && (x.err() == y.err()) && (x.exp() == y.exp())) { 
                return x;
            }
                         
            CORE::BigFloat lower = std::min(CGAL::lower(x), CGAL::lower(y));
            CORE::BigFloat upper = std::max(CGAL::upper(x), CGAL::upper(y));

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
                long new_err = ((err.m()+err.err()) >> shift).longValue()+1; 
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

    struct Singleton :public Unary_function<Interval,bool> {
        bool operator() ( Interval x ) const {       
            return (x.err() == 0); 
        }
    };

    struct Construct :public Binary_function<Boundary,Boundary,Interval>{
        Interval operator()( const Boundary& l,const Boundary& r) const {
            CGAL_precondition( l < r ); 
            return Hull()(l,r);
        }
    };
};




template<typename BFI> long get_significant_bits(BFI bfi);

CORE::BigFloat 
inline 
round(const CORE::BigFloat& x, long rel_prec = CORE::defRelPrec.toLong() ){CGAL_postcondition(rel_prec >= 0);   
    


    // since there is not rel prec defined if Zero_in(x)
    if (x.isZeroIn()) return x; 
    if (CGAL::get_significant_bits(x) <= rel_prec) return x;
   
// if 1 
    BigFloat xr;
    xr.approx(x,rel_prec,1024);
    typedef CORE::BigFloat BF; 
// else       
//     typedef CORE::BigFloat BF; 
//     typedef CORE::BigFloat BFI; 
//     typedef CORE::BigInt Integer;
//     BF xr;
   
//     CORE::BigInt m = x.m();
//     long         err = x.err();
//     long         exp = x.exp(); 
   
//     long shift = ::CORE::bitLength(m) - rel_prec - 1;
//     if( shift > 0 ){    Integer new_m   = m >> shift ; 
//         if(err == 0){        xr = BF(new_m,1,0)*BF::exp2(exp*14+shift);
//         }else{        xr = BF(new_m,2,0)*BF::exp2(exp*14+shift);
//         }
//     }else{    // noting to do
//         xr = x; 
//     }
// endif     

    CGAL_postcondition(CGAL::abs(CGAL::get_significant_bits(xr) - rel_prec) <= 1);   
    CGAL_postcondition(BF(xr.m()-xr.err(),0,xr.exp()) <= BF(x.m()-x.err(),0,x.exp()));
    CGAL_postcondition(BF(xr.m()+xr.err(),0,xr.exp()) >= BF(x.m()+x.err(),0,x.exp()));
    return xr;     
}

template<> class Bigfloat_interval_traits<CORE::BigFloat> 
:public Interval_traits<CORE::BigFloat>
{
public:
    typedef CORE::BigFloat NT;
    typedef CORE::BigFloat BF;

    typedef Bigfloat_interval_traits<NT> Self;

    // How about retuning 
    struct Get_significant_bits {
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
       
    struct Set_precision {
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
     
    struct Get_precision {
        // type for the \c AdaptableGenerator concept.
        typedef long  result_type;  
     
        long operator() () const {
            return  ::CORE::defRelPrec.toLong(); 
        }
    };

    struct Convert_to_bfi {
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
