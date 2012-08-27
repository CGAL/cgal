// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Hemmer
// TODO: add sign to RET


#ifndef CGAL_LEDA_BIGFLOAT_INTERVAL_H
#define CGAL_LEDA_BIGFLOAT_INTERVAL_H

#include <CGAL/basic.h>

#ifndef CGAL_USE_LEDA
#warning This header file needs LEDA installed in order to work properly.
#else // CGAL_USE_LEDA

#include <CGAL/LEDA_basic.h>
#if CGAL_LEDA_VERSION < 500
#include <LEDA/bigfloat.h>
#else
#include <LEDA/numbers/bigfloat.h>
#endif

#include <boost/numeric/interval.hpp>

#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_real.h>
#include <CGAL/leda_bigfloat.h>
#include <CGAL/Interval_traits.h>
#include <CGAL/Bigfloat_interval_traits.h>
#include <CGAL/ipower.h>

namespace CGAL {
namespace internal {

struct Rounding_for_leda_bigfloat {
private:  typedef leda::bigfloat T;
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
      //return a==T(1) && b == T(0);
      return a > b; 
    }
};

} // namespace internal
} //namespace CGAL

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

namespace CGAL {

typedef boost::numeric::interval
< leda::bigfloat,
    boost::numeric::interval_lib::policies
      < internal::Rounding_for_leda_bigfloat,
        internal::Checking_for_leda_bigfloat > > 
leda_bigfloat_interval;

template <> class Algebraic_structure_traits< leda_bigfloat_interval >
    : public Algebraic_structure_traits_base< leda_bigfloat_interval,
                                            Field_with_sqrt_tag >  {
public:
    typedef Tag_false           Is_exact;
    typedef Tag_true            Is_numerical_sensitive;

    class Sqrt
        : public std::unary_function< Type, Type > {
    public:
        Type operator()( const Type& x ) const {
            return ::boost::numeric::sqrt(x);
        }
    };
};

template <> class Real_embeddable_traits< leda_bigfloat_interval >
    : public INTERN_RET::Real_embeddable_traits_base< leda_bigfloat_interval , CGAL::Tag_true > {
public:
  
    class Abs
        : public std::unary_function< Type, Type > {
    public:
        Type operator()( const Type& x ) const {
            return ::boost::numeric::abs(x);
        }
    };

    class To_double
        : public std::unary_function< Type, double > {
    public:
        double operator()( const Type& x ) const {
            return CGAL::to_double(::boost::numeric::median(x));
        }
    };

    class To_interval
        : public std::unary_function< Type, std::pair< double, double > > {
    public:
        std::pair<double, double> operator()( const Type& x ) const {            
            std::pair<double, double> lower_I(CGAL::to_interval(x.lower()));
            std::pair<double, double> upper_I(CGAL::to_interval(x.upper()));
            return std::pair< double, double >(
                    CGAL::min(lower_I.first , upper_I.first ),
                    CGAL::max(lower_I.second, upper_I.second));
        }
    };
};


// Coercion traits:
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short      ,leda_bigfloat_interval)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int        ,leda_bigfloat_interval)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long       ,leda_bigfloat_interval)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float      ,leda_bigfloat_interval)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double     ,leda_bigfloat_interval)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::leda::bigfloat   ,leda_bigfloat_interval)

template <>
struct Coercion_traits< leda_bigfloat_interval , ::leda::integer >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_false Are_implicit_interoperable;
    
    typedef leda_bigfloat_interval Type;
    
    struct Cast{
        typedef Type result_type;
        Type operator()(const leda_bigfloat_interval& x)  const { return x;}
        Type operator()(const ::leda::integer& x) const {
            leda::bigfloat tmp(x);
            leda_bigfloat_interval result(
                    round(tmp,leda::bigfloat::get_precision(),leda::TO_N_INF),
                    round(tmp,leda::bigfloat::get_precision(),leda::TO_P_INF));
            CGAL_postcondition( result.lower() <= x );
            CGAL_postcondition( result.upper() >= x );
            return result; 
        }
    };
};

template <>
struct Coercion_traits< leda_bigfloat_interval , ::leda::rational >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_false Are_implicit_interoperable;
    
    typedef leda_bigfloat_interval Type;
    
    struct Cast{
        typedef Type result_type;
        Type operator()(const leda_bigfloat_interval& x)  const { return x;}
        Type operator()(const ::leda::rational& x) const {
            long prec = ::leda::bigfloat::get_precision();
            leda_bigfloat_interval result (
                    leda_bigfloat::from_rational(x,prec,leda::TO_N_INF),
                    leda_bigfloat::from_rational(x,prec,leda::TO_P_INF));
            CGAL_postcondition( result.lower() <= x );
            CGAL_postcondition( result.upper() >= x );
            return result; 
        }
    };
};

template <>
struct Coercion_traits< leda_bigfloat_interval , ::leda::real >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_false Are_implicit_interoperable;
    
    typedef leda_bigfloat_interval Type;
    
    struct Cast{
        typedef Type result_type;
        Type operator()(const leda_bigfloat_interval& x)  const { return x;}
        Type operator()(const ::leda::real& x) const {
            long current_prec = ::leda::bigfloat::get_precision();
            x.guarantee_relative_error(current_prec);
            leda_bigfloat_interval 
                result(x.get_lower_bound(), x.get_upper_bound());
            CGAL_postcondition( result.lower() <= x );
            CGAL_postcondition( result.upper() >= x );
            return result;
        }
    };
};

template <> struct Coercion_traits< ::leda::integer, leda_bigfloat_interval >
    :public Coercion_traits< leda_bigfloat_interval , ::leda::integer >{}; 

template <> struct Coercion_traits< ::leda::rational, leda_bigfloat_interval >
    :public Coercion_traits< leda_bigfloat_interval , ::leda::rational >{}; 

template <> struct Coercion_traits< ::leda::real, leda_bigfloat_interval >
    :public Coercion_traits< leda_bigfloat_interval , ::leda::real>{};        
        


template<>
class Interval_traits<leda_bigfloat_interval>
  :public internal::Interval_traits_base<leda_bigfloat_interval>
{
public: 
    typedef Interval_traits<leda_bigfloat_interval> Self; 
    typedef leda_bigfloat_interval Interval; 
    typedef leda::bigfloat Bound; 
    typedef CGAL::Tag_true Is_interval; 
    typedef CGAL::Tag_true With_empty_interval; 

    struct Construct :public std::binary_function<Bound,Bound,Interval>{
        Interval operator()( const Bound& l,const Bound& r) const {
            CGAL_precondition( l < r ); 
            return Interval(l,r);
        }
    };

    struct Lower :public std::unary_function<Interval,Bound>{
        Bound operator()( const Interval& a ) const {
            return a.lower();
        }
    };

    struct Upper :public std::unary_function<Interval,Bound>{
        Bound operator()( const Interval& a ) const {
            return a.upper();
        }
    };

    struct Width :public std::unary_function<Interval,Bound>{
        Bound operator()( const Interval& a ) const {
            return ::boost::numeric::width(a);
        }
    };

    struct Median :public std::unary_function<Interval,Bound>{
        Bound operator()( const Interval& a ) const {
            return ::boost::numeric::median(a);
        }
    };
    
    struct Norm :public std::unary_function<Interval,Bound>{
        Bound operator()( const Interval& a ) const {
            return ::boost::numeric::norm(a);
        }
    };

    struct Empty :public std::unary_function<Interval,bool>{
        bool operator()( const Interval& a ) const {
            return ::boost::numeric::empty(a);
        }
    };

    struct Singleton :public std::unary_function<Interval,bool>{
        bool operator()( const Interval& a ) const {
            return ::boost::numeric::singleton(a);
        }
    };

    struct Zero_in :public std::unary_function<Interval,bool>{
        bool operator()( const Interval& a ) const {
            return ::boost::numeric::in_zero(a);
        }
    };

    struct In :public std::binary_function<Bound,Interval,bool>{
        bool operator()( Bound x, const Interval& a ) const {
            return ::boost::numeric::in(x,a);
        }
    };

    struct Equal :public std::binary_function<Interval,Interval,bool>{
        bool operator()( const Interval& a, const Interval& b ) const {
            return ::boost::numeric::equal(a,b);
        }
    };
    
    struct Overlap :public std::binary_function<Interval,Interval,bool>{
        bool operator()( const Interval& a, const Interval& b ) const {
            return ::boost::numeric::overlap(a,b);
        }
    };
    
    struct Subset :public std::binary_function<Interval,Interval,bool>{
        bool operator()( const Interval& a, const Interval& b ) const {
            return ::boost::numeric::subset(a,b);
        }
    };
    
    struct Proper_subset :public std::binary_function<Interval,Interval,bool>{
        bool operator()( const Interval& a, const Interval& b ) const {
            return ::boost::numeric::proper_subset(a,b);
        }
    };
    
    struct Hull :public std::binary_function<Interval,Interval,Interval>{
        Interval operator()( const Interval& a, const Interval& b ) const {
            return ::boost::numeric::hull(a,b);
        }
    };
    
    struct Intersection :public std::binary_function<Interval,Interval,Interval>{
        Interval operator()( const Interval& a, const Interval& b ) const {
            Interval r = ::boost::numeric::intersect(a,b);      
            return r;
        }
    };
};

template<>
class Bigfloat_interval_traits<leda_bigfloat_interval>
    :public Interval_traits<leda_bigfloat_interval> 
{
  typedef leda_bigfloat_interval NT;
  typedef leda::bigfloat BF; 
public:
  typedef Bigfloat_interval_traits<leda_bigfloat_interval> Self;
  typedef CGAL::Tag_true Is_bigfloat_interval; 
  

//   struct Get_significant_bits : public std::unary_function<NT,long>{
//         long operator()( NT x) const {
//             CGAL_precondition(!Singleton()(x));
//             leda::bigfloat lower = x.lower();
//             leda::bigfloat upper = x.upper();
//             leda::integer lower_m = lower.get_significant();
//             leda::integer upper_m = upper.get_significant();
//             leda::integer lower_exp = lower.get_exponent();
//             leda::integer upper_exp = upper.get_exponent();
//             long shift = (upper_exp - lower_exp).to_long();
//             if(shift >= 0 ) upper_m = (upper_m <<  shift);
//             else            lower_m = (lower_m << -shift);
//             //CGAL_postcondition(lower_m.length() == upper_m.length());
//             leda::integer err = upper_m - lower_m; 
//             std::cout <<"LEDA: " << lower_m << " " << err << " " << std::endl; 
//             return CGAL::abs(lower_m.length()-err.length());
//         }
//     };

    
  struct Relative_precision: public std::unary_function<NT,long>{
    long operator()(const NT& x) const {
      CGAL_precondition(!Singleton()(x));
      CGAL_precondition(!CGAL::zero_in(x));

      leda::bigfloat w = Width()(x);
      w = leda::div(w,Lower()(x),Get_precision()(),leda::TO_P_INF); 
      return -leda::ilog2(w).to_long();
    }
  };
  
  struct Set_precision : public std::unary_function<long,long> {
    long operator()( long prec ) const {
      return BF::set_precision(prec); 
    }
  };
     
    struct Get_precision {
        // type for the \c AdaptableGenerator concept.
        typedef long  result_type;  
        long operator()() const {
            return BF::get_precision(); 
        }
    };
};


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


} //namespace CGAL
#endif // CGAL_USE_LEDA
#endif //  CGAL_LEDA_BIGFLOAT_INTERVAL_H
