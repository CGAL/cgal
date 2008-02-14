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

/*! \file NiX/LEDA/interval_support.h
*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_LEDA_INTERVAL_SUPPORT_H
#define CGAL_ALGEBRAIC_KERNEL_D_LEDA_INTERVAL_SUPPORT_H

#if 0

#include <CGAL/basic.h>
#include <CGAL/leda_bigfloat.h>
#include <boost/numeric/interval.hpp>

//#include <NiX/LEDA/integer.h>
//#include <NiX/LEDA/rational.h>
//#include <NiX/LEDA/real.h>

namespace boost{
namespace numeric{
namespace interval_lib{

    template<>
    struct rounded_math<leda::bigfloat>{
    private:    typedef leda::bigfloat T;
    public:
        rounded_math(){};
        ~rounded_math(){};
        
        T conv_down(const T& a){
            return round(a,leda::bigfloat::get_precision(),leda::TO_N_INF);
        };
        T conv_up  (const T& a){
            return round(a,leda::bigfloat::get_precision(),leda::TO_P_INF);
        };  
        // mathematical operations
        T add_down(const T& a, const T& b){
            return add(a,b,leda::bigfloat::get_precision(),leda::TO_N_INF);
        };
        T add_up  (const T& a, const T& b){
            return add(a,b,leda::bigfloat::get_precision(),leda::TO_P_INF);
        };  
        T sub_down(const T& a, const T& b){
            return sub(a, b, leda::bigfloat::get_precision(),leda::TO_N_INF);
        };
        T sub_up  (const T& a, const T& b){
            return sub(a, b, leda::bigfloat::get_precision(),leda::TO_P_INF);
        }; 
        T mul_down(const T& a, const T& b){
            return mul(a, b, leda::bigfloat::get_precision(),leda::TO_N_INF);
        };
        T mul_up  (const T& a, const T& b){
            return mul(a, b, leda::bigfloat::get_precision(),leda::TO_P_INF);
        }; 
        T div_down(const T& a, const T& b){
            return div(a, b, leda::bigfloat::get_precision(),leda::TO_N_INF);
        };
        T div_up  (const T& a, const T& b){
            return div(a, b, leda::bigfloat::get_precision(),leda::TO_P_INF);
        };         
        T sqrt_down(const T& a){
            return sqrt(a, leda::bigfloat::get_precision(),leda::TO_N_INF);
        };
        T sqrt_up  (const T& a){
            return sqrt(a, leda::bigfloat::get_precision(),leda::TO_P_INF);
        }; 

        T median(const T& a, const T& b){ return (a+b)/2;    };   
        T int_down(const T& a)          { return T(floor(a));};   
        T int_up  (const T& a)          { return T(ceil(a)); };

        /*
          T exp_down(T);   
          T exp_up  (T); 
          T log_down(T);  
          T log_up  (T);    
          T cos_down(T);
          T cos_up  (T);  
          T tan_down(T);   
          T tan_up  (T);    
          T asin_down(T);   // [-1;1]
          T asin_up  (T);   // [-1;1]
          T acos_down(T);   // [-1;1]
          T acos_up  (T);   // [-1;1]
          T atan_down(T);   // [-?;+?]
          T atan_up  (T);   // [-?;+?]
          T sinh_down(T);   // [-?;+?]
          T sinh_up  (T);   // [-?;+?]
          T cosh_down(T);   // [-?;+?]
          T cosh_up  (T);   // [-?;+?]
          T tanh_down(T);   // [-?;+?]
          T tanh_up  (T);   // [-?;+?]
          T asinh_down(T);  // [-?;+?]
          T asinh_up  (T);  // [-?;+?]
          T acosh_down(T);  // [1;+?]
          T acosh_up  (T);  // [1;+?]
          T atanh_down(T);  // [-1;1]
          T atanh_up  (T);  // [-1;1]         
        */

        // unprotected rounding class
        //typedef ... unprotected_rounding;
    };

}//namespace interval_lib
/*
inline
std::ostream& operator << 
(std::ostream& os, const boost::numeric::interval<leda::bigfloat>& x)
{
    os << "[" 
       << x.lower().get_significant() << "*2^" << x.lower().get_exponent() 
       << " , "
       << x.upper().get_significant() << "*2^" << x.upper().get_exponent()
       << "]";
    return os;
}*/

}//namespace numeric
}//namespace boost

CGAL_BEGIN_NAMESPACE

namespace CGALi {
    
typedef boost::numeric::interval<leda::bigfloat> leda_bigfloat_interval;

CGAL::Sign inline sign(const leda::bigfloat& x){
    if (x < 0 ) return CGAL::NEGATIVE;
    if (x > 0 ) return CGAL::POSITIVE; 
    return CGAL::ZERO; 
}

long  inline get_significant_bits(const leda_bigfloat_interval& x){
    leda::bigfloat lower = x.lower();
    leda::bigfloat upper = x.upper();
   
    leda::integer lower_m = lower.get_significant();
    leda::integer upper_m = upper.get_significant();
    
    leda::integer lower_exp = lower.get_exponent();
    leda::integer upper_exp = upper.get_exponent();

    long shift = (upper_exp - lower_exp).to_long();
    if(shift >= 0 ) lower_m <<  shift;
    else            upper_m << -shift;

        
    //NiX_postcond(lower_m.length() == upper_m.length());
        
    leda::integer err = lower_m-upper_m; 
    
    return std::max(lower_m.length()-err.length(),0);
}

long  inline set_precision( ::leda::bigfloat, long prec){
    return leda::bigfloat::set_precision(prec); 
}
long  inline get_precision( ::leda::bigfloat){
    return leda::bigfloat::get_precision(); 
}

::leda::bigfloat inline upper(leda_bigfloat_interval x){
    return x.upper();
}

::leda::bigfloat inline lower(leda_bigfloat_interval x){
    return x.lower();
}

leda_bigfloat_interval inline abs(const leda_bigfloat_interval& x){
    return ::boost::numeric::abs(x);
}
leda_bigfloat_interval inline ipower(const leda_bigfloat_interval& x, int i ){
    return ::boost::numeric::pow(x,i);
}
leda_bigfloat_interval inline hull(const leda_bigfloat_interval& x, const leda_bigfloat_interval& y){
    return ::boost::numeric::hull(x,y);
}
leda_bigfloat_interval inline sqrt(const leda_bigfloat_interval& x){
    return ::boost::numeric::sqrt(x);
}
::leda::bigfloat inline width(const leda_bigfloat_interval& x){
    return ::boost::numeric::width(x);
}
::leda::bigfloat inline median(const leda_bigfloat_interval& x){
    return ::boost::numeric::median(x);
}
double inline to_double(const leda_bigfloat_interval& x){
    return CGALi::to_double(median(x));
}
bool inline in_zero(const leda_bigfloat_interval& x){
    return ::boost::numeric::in_zero(x);
}
bool inline overlap(const leda_bigfloat_interval& x, const leda_bigfloat_interval& y){
    return ::boost::numeric::overlap(x,y);
}

::leda::bigfloat inline relative_error(const leda_bigfloat_interval& x){
    if(in_zero(x)){
        return CGALi::abs(x).upper();
    }else{
        return  (CGALi::width(x) / CGALi::abs(x)).upper();
    }
}

bool inline singleton( const leda_bigfloat_interval& a ) {
    return ::boost::numeric::singleton(a); 
}


leda_bigfloat_interval  inline convert_to_bfi(const leda::real& x) {
    long current_prec = ::leda::bigfloat::get_precision();
    //x.improve_approximation_to(current_prec);
    x.guarantee_relative_error(current_prec);
    
    leda::bigfloat bnum = x.to_bigfloat();  
    leda::bigfloat berr = x.get_bigfloat_error();
    
    leda::bigfloat low = leda::sub(bnum,berr,current_prec,LEDA::TO_N_INF);
    leda::bigfloat upp = leda::add(bnum,berr,current_prec,LEDA::TO_P_INF);
    leda_bigfloat_interval bfi(low,upp) ;
    
    
//     std::cout <<"x: "<<  x << std::endl;
//     std::cout <<"bfi.lower(): "<<  bfi.lower() << std::endl;
//     std::cout <<"bfi.upper(): "<<  bfi.upper() << std::endl;

    CGAL_postcondition( bfi.lower() <= x );
    CGAL_postcondition( bfi.upper() >= x );
    
    return bfi; 
}


leda_bigfloat_interval  inline convert_to_bfi(const ::leda::integer& x) {
    long current_prec = ::leda::bigfloat::get_precision();
    leda_bigfloat_interval bfi;
    long length = x.length();
    
    if(length > current_prec) {
        ::leda::integer significant = CGAL::abs(x) >> (length - current_prec);
        ::leda::bigfloat lower,upper;
        if(x > 0){
            lower = ::leda::bigfloat(significant,length - current_prec);
            upper = ::leda::bigfloat(significant+1,length - current_prec);
        }else{
            lower = -::leda::bigfloat(significant+1,length - current_prec);
            upper = -::leda::bigfloat(significant,length - current_prec);
        }
        bfi = leda_bigfloat_interval(lower,upper);
    }else{
        ::leda::bigfloat bf(x);
        bfi = leda_bigfloat_interval(bf,bf);
    }
    CGAL_postcondition( bfi.lower() <= x );
    CGAL_postcondition( bfi.upper() >= x );
    return bfi; 
}

leda_bigfloat_interval  inline convert_to_bfi(const ::leda::rational& x) {
    long old_prec = ::leda::bigfloat::get_precision();
    ::leda::bigfloat::set_precision(old_prec*2);
    leda_bigfloat_interval num = convert_to_bfi(x.numerator());
    leda_bigfloat_interval den = convert_to_bfi(x.denominator());
    ::leda::bigfloat::set_precision(old_prec);
    leda_bigfloat_interval bfi = num/den;
    
    CGAL_postcondition( bfi.lower() <= x );
    CGAL_postcondition( bfi.upper() >= x );
    return bfi; 
}
 


} // namespace CGALi

inline
std::pair<double, double> to_interval( const CGALi::leda_bigfloat_interval& x ) {
    return std::pair< double, double >( CGAL::to_double( x.lower() ),
                                        CGAL::to_double( x.upper() ) );
}


CGAL_END_NAMESPACE

#endif

#endif // CGAL_ALGEBRAIC_KERNEL_D_LEDA_INTERVAL_SUPPORT_H
