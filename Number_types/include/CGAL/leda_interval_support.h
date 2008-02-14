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


/*! \file CGAL/leda_interval_support.h
*/

#ifndef CGAL_LEDA_INTERVAL_SUPPORT_H
#define CGAL_LEDA_INTERVAL_SUPPORT_H

#include <CGAL/basic.h>
#include <CGAL/leda_bigfloat.h>
#include <boost/numeric/interval.hpp>

#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_real.h>

// namespace NiX {
CGAL_BEGIN_NAMESPACE

// namespace Intern {
namespace CGALi {


struct Rounding_for_leda_bigfloat {
private:    typedef leda::bigfloat T;
public:
    Rounding_for_leda_bigfloat(){};
    ~Rounding_for_leda_bigfloat(){};
    
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

class Checking_for_leda_bigfloat {
        
        typedef leda::bigfloat T;

    public:
        
        static T pos_inf() {
            T b = T(5) / T(0);
            CGAL_assertion(leda::ispInf(b));
            return b;
        }

        static T neg_inf() {
            T b = T(-5) / T(0);
            CGAL_assertion(leda::isnInf(b));
            return b;
        }
        
        static T nan() {
            T b = T(0)*pos_inf();
            CGAL_assertion(leda::isNaN(b));
            return b;
        }

        static bool is_nan(const T& b) {
            return leda::isNaN(b);
        }

        static T empty_lower() {
            return T(1);
        }

        static T empty_upper() {
            return T(0);
        }

        static bool is_empty(const T& a, const T& b) {
            return a==T(1) && b == T(0);
        }
    };

} // namespace CGALi
// } // namespace NiX
CGAL_END_NAMESPACE

namespace boost {
namespace numeric {


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
    }


}//namespace numeric
}//namespace boost

// namespace NiX {
CGAL_BEGIN_NAMESPACE

typedef boost::numeric::interval
  < leda::bigfloat,
    boost::numeric::interval_lib::policies
      < CGALi::Rounding_for_leda_bigfloat,
        CGALi::Checking_for_leda_bigfloat > > 
    leda_bigfloat_interval;




/*
CGAL::Sign inline sign(const leda::bigfloat& x){
    if (x < 0 ) return CGAL::NEGATIVE;
    if (x > 0 ) return CGAL::POSITIVE; 
    return CGAL::ZERO; 
}
*/

leda_bigfloat_interval inline ipower(const leda_bigfloat_interval& x, int i ){
    return ::boost::numeric::pow(x,i);
}

::leda::bigfloat inline median(const leda_bigfloat_interval& x){
    return ::boost::numeric::median(x);
}

// forward
 template<typename BFI> class Bigfloat_interval_traits;
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

template <> class Algebraic_structure_traits< leda_bigfloat_interval >
  : public Algebraic_structure_traits_base< leda_bigfloat_interval,
                                            Field_with_sqrt_tag >  {
  public:
    typedef Tag_false           Is_exact;
    typedef Tag_true            Is_numerical_sensitive;

    class Sqrt
      : public Unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return ::boost::numeric::sqrt(x);
        }
    };
};

template <> class Real_embeddable_traits< leda_bigfloat_interval >
  : public Real_embeddable_traits_base< leda_bigfloat_interval > {
  public:

    class Abs
      : public Unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
            return ::boost::numeric::abs(x);
        }
    };

    class To_double
      : public Unary_function< Type, double > {
      public:
        double operator()( const Type& x ) const {
          return CGAL::to_double(::boost::numeric::median(x));
        }
    };

    class To_interval
      : public Unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {
// TODO
//           return ::boost::numeric::hull(CGAL::to_interval(x.lower()),CGAL::to_interval(x.upper()));

    return std::pair< double, double >( CGAL::to_double( x.lower() ),
                                        CGAL::to_double( x.upper() ) );
        }
    };
};



//  template <>
//    class NT_traits<leda_bigfloat_interval> 
//    : public NT_traits_base<leda_bigfloat_interval, 
//                           Field_with_sqrt_tag>,
//      public NT_traits_comparable_base<leda_bigfloat_interval>
// 
//     {
//     public:
//         typedef leda_bigfloat_interval NT;
//         
//         class Abs {
//         public:
//             // type for the \c AdaptableUnaryFunction concept.
//             typedef NT  argument_type;
//             // type for the \c AdaptableUnaryFunction concept.
//             typedef NT  result_type;
//             NT operator()( NT a) const {
//                 return ::boost::numeric::abs(a);
//             };
//         };
//         
//         class Sqrt {
//         public:
//             // type for the \c AdaptableUnaryFunction concept.
//             typedef NT  argument_type;
//             // type for the \c AdaptableUnaryFunction concept.
//             typedef NT  result_type;
//             NT operator()( NT a) const {
//                 return ::boost::numeric::sqrt(a);
//             };
//         };
// 
//         class To_double {
//         public:
//             // the result type.
//             typedef double result_type;
//             // the argument type.
//             typedef NT argument_type;
//             
//             double operator()(const NT& a) const {
//               return NiX::to_double(::boost::numeric::median(a));
//             } 
//         };
//         class To_Interval {
//         public:
//             // the result type.
//             typedef Interval result_type;
//             // the argument type.
//             typedef NT argument_type;
// 
//             Interval operator()(const NT& a) const {
//               return ::boost::numeric::hull(NiX::to_Interval(a.lower()),
//                                             NiX::to_Interval(a.upper()));
//             } 
//         };

//TODO
//porting Floor_log2_abs and Ceil_log2_abs from EXACUS2CGAL
//         class Floor_log2_abs {
//         public:
//           typedef NT argument_type;
//           typedef long result_type;
//           result_type operator() (argument_type x) {
//             CGAL_precondition(! ::boost::numeric::in_zero(x));
//             return floor_log2_abs(::boost::numeric::abs(x).lower());
//           }
//         };
//         class Ceil_log2_abs {
//         public:
//           typedef NT argument_type;
//           typedef long result_type;
//           result_type operator() (argument_type x) {
//             CGAL_precondition(!(::boost::numeric::in_zero(x) && 
//                           ::boost::numeric::singleton(x)));
//             return ceil_log2_abs(::boost::numeric::abs(x).upper());
//           }
//         };
//     };

::leda::bigfloat inline relative_error(const leda_bigfloat_interval& x){
    if(in_zero(x)){
        return CGAL::abs(x).upper();
    }else{
        return (width(x) / CGAL::abs(x)).upper();
    }
}

// } // namespace NiX
CGAL_END_NAMESPACE


#endif //  CGAL_LEDA_INTERVAL_SUPPORT_H
