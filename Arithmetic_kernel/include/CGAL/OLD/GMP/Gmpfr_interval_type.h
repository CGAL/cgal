// TODO: add sign to RET


#ifndef CGAL_GMPFR_INTERVAL_TYPE_H
#define CGAL_GMPFR_INTERVAL_TYPE_H

#include <CGAL/basic.h>
#include <CGAL/Gmpfr.h>
#include <boost/numeric/interval.hpp>

#define OP(o,d) Gmpfr r;                                \
  mpfr_ ## o (r.fr(), a.fr(), b.fr(), GMP_RND ## d);    \
  return r;

namespace CGAL {

namespace internal {

struct Rounding_for_gmpfr {
private:  typedef CGAL::Gmpfr T;
public:
  Rounding_for_gmpfr(){};
  ~Rounding_for_gmpfr(){};
    
  T conv_down(const T& a){ // TODO: fix this
    return add_down(T(0),a); 
  };
  T conv_up  (const T& a){ // TODO: fix this
    return add_up(T(0),a); 
  };  
  // mathematical operations
  T add_down(const T& a, const T& b) {
    OP(add,D);
  };
  T add_up  (const T& a, const T& b){
    OP(add,U);
  };  
  T sub_down(const T& a, const T& b){
    OP(sub,D);
  };
  T sub_up  (const T& a, const T& b){
    OP(sub,U);
  }; 
  T mul_down(const T& a, const T& b){
    OP(mul,D);
  };
  T mul_up  (const T& a, const T& b){
    OP(mul,U);
  }; 
  T div_down(const T& a, const T& b){
    OP(div,D);
  };
  T div_up  (const T& a, const T& b){
    OP(div,U);
  };

  T sqrt_down(const T& a){
    T b;
    mpfr_sqrt(b.fr(), a.fr(), GMP_RNDD);
    return b;
  };
  T sqrt_up  (const T& a){
    T b;
    mpfr_sqrt(b.fr(), a.fr(), GMP_RNDU);
    return b;
  }; 

  T median(const T& a, const T& b) { 
    T result(a + (b-a)/2);
    CGAL_postcondition(a <= result);
    CGAL_postcondition(result <= b);
    return result;
  };   
  T int_down(const T& a) { 
    T r;
    mpfr_floor(r.fr(), a.fr());
    return r;
  };   
  T int_up  (const T& a) {
    T r;
    mpfr_ceil(r.fr(), a.fr());
    return r;
  };
};

class Checking_for_gmpfr {
        
  typedef CGAL::Gmpfr T;

public:
        
  static T pos_inf() {
    T b;
    mpfr_set_inf(b.fr(), 1);
    return b;
  }

  static T neg_inf() {
    T b;
    mpfr_set_inf(b.fr(), -1);
    return b;
  }
        
  static T nan() {
    T b;
    mpfr_set_nan(b.fr());
    return b;
  }

  static bool is_nan(const T& b) {
    return mpfr_nan_p(b.fr());
  }

  static T empty_lower() {
    return T(1);
  }

  static T empty_upper() {
    return T(0);
  }

  static bool is_empty(const T& a, const T& b) {
//        return a==T(1) && b == T(0);
//        return false;
    return a > b; // TODO: optimize this
  }
};

} // namespace internal
} //namespace CGAL

namespace boost {
namespace numeric {
inline
std::ostream& operator << 
    (std::ostream& os, const boost::numeric::interval<CGAL::Gmpfr>& x)
{
  os << "[" 
     << x.lower() << ", " << x.upper() << "]";
  return os;
}


}//namespace numeric
}//namespace boost

namespace CGAL {

typedef boost::numeric::interval< 
        Gmpfr, 
        boost::numeric::interval_lib::policies
      < internal::Rounding_for_gmpfr, internal::Checking_for_gmpfr > > 
Gmpfr_interval;

// I have to redfine these operators, since the test reports an 
// ambiguity with operators in CGAL::Polynomial 
// However, it seems that this is due to the use of struct interval_holder 
// in the definition of the operators in boost. 

inline bool operator == (const Gmpfr_interval& x , const Gmpfr_interval& y)
{return x.operator==(y) ;}

} //namespace CGAL
//#endif // CGAL_USE_LEDA
#undef OP
#endif //  CGAL_GMPFR_INTERVAL_TYPE_H
