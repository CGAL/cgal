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

#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_real.h>
#include <CGAL/leda_bigfloat.h>
#include <boost/numeric/interval.hpp>

CGAL_BEGIN_NAMESPACE
namespace CGALi {

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
        return a==T(1) && b == T(0);
    }
};

} // namespace CGALi
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

CGAL_BEGIN_NAMESPACE

typedef boost::numeric::interval
< leda::bigfloat,
    boost::numeric::interval_lib::policies
      < CGALi::Rounding_for_leda_bigfloat,
        CGALi::Checking_for_leda_bigfloat > > 
leda_bigfloat_interval;

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
            std::pair<double, double> lower_I(CGAL::to_interval(x.lower()));
            std::pair<double, double> upper_I(CGAL::to_interval(x.upper()));
            return std::pair< double, double >(
                    CGAL::min(lower_I.first , upper_I.first ),
                    CGAL::max(lower_I.second, upper_I.second));
        }
    };
};


// Coercion traits:
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short      ,leda_bigfloat_interval);
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int        ,leda_bigfloat_interval);
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long       ,leda_bigfloat_interval);
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float      ,leda_bigfloat_interval);
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double     ,leda_bigfloat_interval);
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::leda::bigfloat   ,leda_bigfloat_interval);

template <>
struct Coercion_traits< leda_bigfloat_interval , ::leda::integer >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_false Are_implicit_interoperable;
    
    typedef leda_bigfloat_interval Type;
    
    struct Cast{
        typedef Type result_type;
        Type operator()(const leda_bigfloat_interval& x)  const { return x;}
        Type operator()(const ::leda::integer x) const {
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
        Type operator()(const ::leda::rational x) const {
            leda_bigfloat_interval result (
                    leda_bigfloat::from_rational(x, ::leda::bigfloat::get_precision(), leda::TO_N_INF),
                    leda_bigfloat::from_rational(x, ::leda::bigfloat::get_precision(), leda::TO_P_INF));
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
        Type operator()(const ::leda::real x) const {
            long current_prec = ::leda::bigfloat::get_precision();
            x.guarantee_relative_error(current_prec);
            leda_bigfloat_interval result(x.get_lower_bound(), x.get_upper_bound());
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
        



CGAL_END_NAMESPACE
#endif // CGAL_USE_LEDA
#endif //  CGAL_LEDA_BIGFLOAT_INTERVAL_H
