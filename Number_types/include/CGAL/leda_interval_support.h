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
class Interval_traits<leda_bigfloat_interval>
{
public: 
    typedef Interval_traits<leda_bigfloat_interval> Self; 
    typedef leda_bigfloat_interval Interval; 
    typedef leda::bigfloat Boundary; 
    typedef CGAL::Tag_true Is_interval; 
    typedef CGAL::Tag_true With_empty_interval; 

    struct Construct :public Binary_function<Boundary,Boundary,Interval>{
        Interval operator()( const Boundary& l,const Boundary& r) const {
            CGAL_precondition( l < r ); 
            return Interval(l,r);
        }
    };

    struct Lower :public Unary_function<Interval,Boundary>{
        Boundary operator()( const Interval& a ) const {
            return a.lower();
        }
    };

    struct Upper :public Unary_function<Interval,Boundary>{
        Boundary operator()( const Interval& a ) const {
            return a.upper();
        }
    };

    struct Width :public Unary_function<Interval,Boundary>{
        Boundary operator()( const Interval& a ) const {
            return ::boost::numeric::width(a);
        }
    };

    struct Median :public Unary_function<Interval,Boundary>{
        Boundary operator()( const Interval& a ) const {
            return ::boost::numeric::median(a);
        }
    };
    
    struct Norm :public Unary_function<Interval,Boundary>{
        Boundary operator()( const Interval& a ) const {
            return ::boost::numeric::norm(a);
        }
    };

    struct Empty :public Unary_function<Interval,bool>{
        bool operator()( const Interval& a ) const {
            return ::boost::numeric::empty(a);
        }
    };

    struct Singleton :public Unary_function<Interval,bool>{
        bool operator()( const Interval& a ) const {
            return ::boost::numeric::singleton(a);
        }
    };

    struct Zero_in :public Unary_function<Interval,bool>{
        bool operator()( const Interval& a ) const {
            return ::boost::numeric::in_zero(a);
        }
    };

    struct In :public Binary_function<Boundary,Interval,bool>{
        bool operator()( Boundary x, const Interval& a ) const {
            return ::boost::numeric::in(x,a);
        }
    };

    struct Equal :public Binary_function<Interval,Interval,bool>{
        bool operator()( const Interval& a, const Interval& b ) const {
            return ::boost::numeric::equal(a,b);
        }
    };
    
    struct Overlap :public Binary_function<Interval,Interval,bool>{
        bool operator()( const Interval& a, const Interval& b ) const {
            return ::boost::numeric::overlap(a,b);
        }
    };
    
    struct Subset :public Binary_function<Interval,Interval,bool>{
        bool operator()( const Interval& a, const Interval& b ) const {
            return ::boost::numeric::subset(a,b);
        }
    };
    
    struct Proper_subset :public Binary_function<Interval,Interval,bool>{
        bool operator()( const Interval& a, const Interval& b ) const {
            return ::boost::numeric::proper_subset(a,b);
        }
    };
    
    struct Hull :public Binary_function<Interval,Interval,Interval>{
        Interval operator()( const Interval& a, const Interval& b ) const {
            return ::boost::numeric::hull(a,b);
        }
    };
    
    struct Intersection :public Binary_function<Interval,Interval,Interval>{
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
    
public:
    typedef Bigfloat_interval_traits<leda_bigfloat_interval> Self;
    typedef leda_bigfloat_interval NT;
    typedef leda::bigfloat BF;

    struct Get_significant_bits : public Unary_function<NT,long>{

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
  
    struct Set_precision : public Unary_function<long,long> {
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

/*
  struct Convert_to_bfi {
  
  typedef NT result_type;

  NT operator()( const leda::real& x ) {
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


  NT operator()(const ::leda::integer& x) {
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


  NT operator()(const ::leda::rational& x) {
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
*/
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


CGAL_END_NAMESPACE
#endif // CGAL_USE_LEDA
#endif //  CGAL_LEDA_INTERVAL_SUPPORT_H
