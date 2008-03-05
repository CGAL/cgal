// TODO:
// - add to comming interval_traits: 
// -- CGAL::median(bf_interval)
// -- CGAL::ipower(bf_interval)
// -- CGAL::relatoive_error(bf_interval)

/*! \file CGAL/leda_interval_support.h
 */

#ifndef CGAL_LEDA_INTERVAL_SUPPORT_H
#define CGAL_LEDA_INTERVAL_SUPPORT_H

#include <CGAL/basic.h>

#ifndef CGAL_USE_LEDA
#warning This header file needs LEDA installed in order to work properly.
#else // CGAL_USE_LEDA

#include <CGAL/leda_bigfloat.h>
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_real.h>
#include <CGAL/leda_bigfloat_interval.h>


#include <CGAL/interval_support.h>

CGAL_BEGIN_NAMESPACE

template<>
class Bigfloat_interval_traits<leda_bigfloat_interval>
{
public:
    typedef leda_bigfloat_interval NT;

    typedef leda::bigfloat BF;

    class Get_significant_bits {
    public:
        // type for the \c AdaptableUnaryFunction concept.
        typedef NT  argument_type;
        // type for the \c AdaptableUnaryFunction concept.
        typedef long  result_type;

        long operator()( NT x) const {
            leda::bigfloat lower = x.lower();
            leda::bigfloat upper = x.upper();

            leda::integer lower_m = lower.get_significant();
            leda::integer upper_m = upper.get_significant();
             
            leda::integer lower_exp = lower.get_exponent();
            leda::integer upper_exp = upper.get_exponent();
             
            long shift = (upper_exp - lower_exp).to_long();
            if(shift >= 0 ) upper_m = (upper_m <<  shift);
            else            lower_m = (lower_m << -shift);
             
            //CGAL_postcondition(lower_m.length() == upper_m.length());
             
            leda::integer err = lower_m-upper_m; 
             
            return std::max(lower_m.length()-err.length(),0);
             
        }
    };
     
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

    class Set_precision {
    public:
        // type for the \c AdaptableUnaryFunction concept.
        typedef long  argument_type;
        // type for the \c AdaptableUnaryFunction concept.
        typedef long  result_type;  
     
        long operator() ( long prec ) const {
            return BF::set_precision(prec); 
        }
    };
     
    class Get_precision {
    public:
        // type for the \c AdaptableGenerator concept.
        typedef long  result_type;  
     
        long operator() () const {
            return BF::get_precision(); 
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

        NT operator() ( const leda::real& x ) {
            long current_prec = ::leda::bigfloat::get_precision();
            //x.improve_approximation_to(current_prec);
            x.guarantee_relative_error(current_prec);
             
            leda::bigfloat bnum = x.to_bigfloat();  
            leda::bigfloat berr = x.get_bigfloat_error();
             
            leda::bigfloat low 
                = leda::sub(bnum,berr,current_prec,LEDA::TO_N_INF);
            leda::bigfloat upp 
                = leda::add(bnum,berr,current_prec,LEDA::TO_P_INF);
            leda_bigfloat_interval bfi(low,upp) ;
             
            //     std::cout <<"x: "<<  x << std::endl;
            //     std::cout <<"bfi.lower(): "<<  bfi.lower() << std::endl;
            //     std::cout <<"bfi.upper(): "<<  bfi.upper() << std::endl;

            CGAL_postcondition( bfi.lower() <= x );
            CGAL_postcondition( bfi.upper() >= x );
             
            return bfi; 
        }


        NT operator() (const ::leda::integer& x) {
            long current_prec = ::leda::bigfloat::get_precision();
            leda_bigfloat_interval bfi;
            long length = x.length();
             
            if(length > current_prec) {
                ::leda::integer significant 
                    = CGAL::abs(x) >> (length - current_prec);
                ::leda::bigfloat lower,upper;
                if(x > 0){
                    lower = ::leda::bigfloat(significant,length - current_prec);
                    upper = ::leda::bigfloat(significant+1,length - current_prec);
                }else{
                    lower 
                        = -::leda::bigfloat(significant+1,length - current_prec);
                    upper 
                        = -::leda::bigfloat(significant,length - current_prec);
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


        NT operator() (const ::leda::rational& x) {
            long old_prec = ::leda::bigfloat::get_precision();
            ::leda::bigfloat::set_precision(old_prec*2);
            Bigfloat_interval_traits<NT>::Convert_to_bfi convert_to_bfi;
            leda_bigfloat_interval num = convert_to_bfi(x.numerator());
            leda_bigfloat_interval den = convert_to_bfi(x.denominator());
            ::leda::bigfloat::set_precision(old_prec);
            leda_bigfloat_interval bfi = num/den;
            CGAL_postcondition( bfi.lower() <= x );
            CGAL_postcondition( bfi.upper() >= x );
            return bfi; 
        }
         
    };
};

// left overs? 
::leda::bigfloat inline median(const leda_bigfloat_interval& x){
    return ::boost::numeric::median(x);
}

::leda::bigfloat inline relative_error(const leda_bigfloat_interval& x){
    if(in_zero(x)){
        return CGAL::abs(x).upper();
    }else{
        return (width(x) / CGAL::abs(x)).upper();
    }
}

leda_bigfloat_interval inline ipower(const leda_bigfloat_interval& x, int i ){
    return ::boost::numeric::pow(x,i);
}


CGAL_END_NAMESPACE
#endif // CGAL_USE_LEDA
#endif //  CGAL_LEDA_INTERVAL_SUPPORT_H
