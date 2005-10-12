#ifndef CGAL_POLYNOMIAL_SIMPLE_INTERVAL_ROOT_H
#define CGAL_POLYNOMIAL_SIMPLE_INTERVAL_ROOT_H
#include <CGAL/Polynomial/basic.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE;

//! A root represented as a bounding interval and a polynomial.
/*!
  Representing an interval which contains one root of a function. 
  \todo cache sturm sequence
*/
template <class Traits>
class Simple_interval_root {
  typedef Simple_interval_root<Traits> This;
  //! The bit field for storing the type
  /*!
    - UP means the first non-zero derivative is positive
    - EVEN the multiplicity of the root is even
    - INF the root is infinite
    - CONST the root is a rational number
    
    The valid combinations are UP, UP|EVEN, UP | INF,
    EVEN|CONST. EVEN|INF is the unitialized value.
  */
  //typedef enum Fields {UP=1, EVEN=2, INF=4, CONST=8};
  /*typedef enum Type {UP, DOWN, 
		     EVEN_DOWN, EVEN_UP, 
		     POS_INF, NEG_INF, 
		     EVEN_CONST, CONST, UNINITIALIZED} Type;*/

  typedef enum Fields {UP=1, EVEN=2, INF=4, CONST=8} Type_fields;
  typedef unsigned char Type;
  typedef typename Traits::Function Polynomial;
  typedef typename Traits::NT NT;
  typedef typename Traits::Isolating_interval Interval;
  //typedef internal::Isolating_interval_tools<Polynomial, NT, Interval> IIT;
public:
  Simple_interval_root(){set_type(INF|EVEN); assert(is_null());}
  /*template <class RNT>
  Simple_interval_root(const RNT &nt): ii_(nt){
    bool is_this_used;
    set_type(CONST);
    audit();
    compute_approximation();
    }*/

  template <class CNT>
  Simple_interval_root(CNT nt) {
    if (std::numeric_limits<CNT>::has_infinity && (nt == std::numeric_limits<CNT>::infinity() 
						    || -nt == std::numeric_limits<CNT>::infinity())){
      if (nt == std::numeric_limits<CNT>::infinity() ){
	set_type(INF|UP);
      } else {
	set_type(INF);
      }
    }else {
      set_type(CONST);
      ii_= Interval(nt); 
    } 
    audit();
    compute_approximation();
  }
  
  Simple_interval_root(Type t): type_(t){
    audit();
    compute_approximation();
  };

  //! represent a rational root
  Simple_interval_root(const NT &nt, bool is_odd, Traits k): ii_(nt),  kernel_(k){
    if (is_odd) set_type(CONST);
    else set_type(CONST | EVEN);
    audit();
    compute_approximation();
  }

  //! Represent a rational root, another way. This needs to be public since  intervals can be opaque
  Simple_interval_root(const Interval &nt, bool is_odd=true, Traits k= Traits()): ii_(nt), 
										  kernel_(k){
    if (is_odd) set_type(CONST);
    else set_type(CONST | EVEN);
    compute_approximation();
    audit();
  }

  //! Represent a root by an interval and a polynomial
  Simple_interval_root(const Interval &ii, const Polynomial &sa, 
		       Sign slb, Sign sub, 
		       Traits k): ii_(ii), function_(sa),
				  kernel_(k){
    CGAL_Polynomial_precondition(!interval().is_singular());
    //Sign slb= sign_(ii_.lower_bound());
    if (slb == sub){
      if (slb== POSITIVE){
	set_type(EVEN|UP);
      } else {
	set_type(EVEN);
      }
    } else if (slb == POSITIVE){
      set_type(0);
    } else if (slb== NEGATIVE){
      set_type(UP);
    } else {
      set_type(INF|EVEN);
      CGAL_Polynomial_assertion(0);
    }
    compute_approximation();
    audit();
  }

  static This infinity(){
    This ret(Type(UP|INF));
    CGAL_Polynomial_postcondition(ret.is_infinite());
    //ret.compute_approximation();
    return ret;
  }

  bool operator<(const This &o) const {
    CGAL_Polynomial_expensive_precondition(!is_null() && !o.is_null());
    Comparison_result r= compare(o);
    audit(); o.audit();
    return r==SMALLER;
  }
  
  bool operator>(const This &o) const {
    CGAL_Polynomial_expensive_precondition(!is_null() && !o.is_null());
    Comparison_result r= compare(o);
    audit(); o.audit();
    return r==LARGER;
  }
  
  bool operator!=(const This &o) const {
    /*if (is_null()) return !o.is_null();
      else if (o.is_null()) return true;*/
    CGAL_Polynomial_expensive_precondition(!is_null() && !o.is_null());
    Comparison_result r= compare(o);
    audit(); o.audit();
    return r != EQUAL;
  }
  bool operator<=(const This &o) const {
    CGAL_Polynomial_expensive_precondition(!is_null() && !o.is_null());
    Comparison_result r= compare(o);
    audit(); o.audit();
    return r != LARGER;
  }
  bool operator>=(const This &o) const {
    CGAL_Polynomial_expensive_precondition(!is_null() && !o.is_null());
    Comparison_result r= compare(o);
    audit(); o.audit();
    return r != SMALLER;
  }
  
  bool operator==(const This &o) const {
    /*if (is_null()) return o.is_null();
      else if (o.is_null()) return false;*/
    CGAL_Polynomial_expensive_precondition(!is_null() && !o.is_null());
    Comparison_result r= compare(o);
    audit(); o.audit();
    return r==EQUAL;
  }

  bool is_even_multiplicity() const {
    return type_&EVEN;
  }

  //! \todo implement multiplicity
  int multiplicity() const {
    Polynomial_expensive_precondition(!is_null());
    bool I_have_not_implemented_this;
    //if (is_odd_) return 1;
    //else return 2;
    return 0;
  }
  
  //! Compute a value as a double
  /*!
    This currently does not compute the closest double. I should figure out what value to use
    for accuracy to make the result be the closest double.
    \todo compute closest double rather than stupid approximation
  */
  double double_approximation(double accuracy=.000001) const {
    CGAL_Polynomial_expensive_precondition(!is_null());
    return compute_double(accuracy);
  }

  //! Represent an interval by an exact number type.
  /*!
    I forget why I use this. 
  */
  std::pair<NT, NT> isolating_interval() const {
    bool do_not_use;
    CGAL_Polynomial_precondition(!is_infinite()); 
    CGAL_Polynomial_expensive_precondition(!is_null());
    return interval().to_exact_interval();
  }

  //! To use for interval arithmetic
  /*!
    This refines the current interval too.
  */
  std::pair<double, double> double_interval(double accuracy=.001) const {
    CGAL_Polynomial_expensive_precondition(!is_null());
    return compute_interval(accuracy);
  }
  
  //! Write lots of info about the interval
  void write(std::ostream &o) const {
#ifndef NDEBUG
    This t=*this;
    t.write_internal(o);
#else
    write_internal(o);
#endif
  }

  void print() const {
    write(std::cout);
  }
  
  //! Negate the interval. 
  This operator-() const {
    CGAL_Polynomial_expensive_precondition(!is_null());
    if (is_pos_inf()) return This(Type(INF));
    else if (is_neg_inf()) return infinity();
    else if (is_rational()){
      return This(-interval(), !is_even_multiplicity(), kernel_);
    } else {
      This copy= *this;
      copy.ii_= -ii_;
      typename Traits::Negate_variable nf= kernel_.negate_variable_object();
      copy.function_= nf(function_);
      return copy;
    }
  }

  //! Return true if the root is known to be a rational number
  bool is_rational() const {
    return type_&CONST;
  }
  NT to_rational() const {
    CGAL_Polynomial_precondition(is_rational());
    return ii_.to_nt();
  }
  //! Return true if the root is +/- infinity.
  bool is_infinite() const {
    return type_&INF;
  }
  //! Return true if the root is a real root
  bool is_normal() const {
    return (is_up() || is_down()) && (!is_infinite());
  }
  //! This is needed by the solvers
  Interval interval() const {
    return ii_;
  }
protected:

  void write_internal(std::ostream &o) const {
   if (is_pos_inf()){
	o << "inf";
    } else if (is_neg_inf()){
      o << "-inf";
    } else {
      o<< interval();
      if (is_rational()){
	if (is_even_multiplicity()){
	  o<<"(Even)";
	}
      } else {
	if (type_ & EVEN) {
	  if (type_& UP) o <<"[--]";
	  else o << "[++]";
	} else {
	  if (type_& UP) o << "[-+]";
	  else o << "[+-]";
	}
	o<< " = " << immutable_double_approximation();
	o << "(" << function_ << ")";
      }
    }
  }

  void set_type(Type t) const {
    type_=t;
  }

  Type type() const {
    return type_;
  }

  bool is_pos_inf() const {
    return type_&INF && is_up();
  }
  bool is_neg_inf() const {
    return type_&INF && !(is_up());
  }

  void set_is_rational(const Interval &i) const {
    set_interval(i);
    set_type(CONST);
    function_=Polynomial();
  }

  std::pair<double, double> compute_interval(double accuracy) const {
    if (type_&INF) {
      if (is_up()){
	return std::pair<double, double>(double_inf_rep(), double_inf_rep());
      } else {
	return std::pair<double, double>(-double_inf_rep(), -double_inf_rep());
      }
    }
    
    double oaw;//= interval().approximate_width();
    while (interval().approximate_relative_width() > accuracy){
      oaw= interval().approximate_width();
      std::pair<double,double> before= to_interval(interval());
      refine();
      std::pair<double,double> after= to_interval(interval());
      CGAL_assertion(oaw != interval().approximate_width());
      /*if (oaw == interval().approximate_width()){
	break;
	}*/
    }
    return to_interval(interval());
  }

  double compute_double(double accuracy) const {
    if (is_infinite()) {
      if (is_up()){
	return double_inf_rep();
      } else {
	return -double_inf_rep();
      }
    }
    //This t= *this;
    std::pair<double, double> i= double_interval(accuracy);
    return (i.first+i.second)/2.0;
  }

  bool contains_root(const Interval &i, Sign ls, Sign us) const {
     if (!is_even_multiplicity()){
       return (is_up() && ls==NEGATIVE && us==POSITIVE) || (is_down() && ls==POSITIVE && us==NEGATIVE);
    } else {
      return (i.apply_to_interval(kernel_.Sturm_root_count_object(function_)) != 0);
    }
  }

  void refine() const {
    CGAL_Polynomial_precondition(is_normal());
    Sign sn= interval().apply_to_midpoint(sign_at());
    if (sn == ZERO){
      set_is_rational(interval().midpoint_interval());
    } else if (contains_root(interval().first_half(), lower_sign(), sn)){
      set_interval(interval().first_half());
    } else {
      set_interval(interval().second_half());
    }
  }

  void refine_using(const Interval &o) const {
    int noi= interval().number_overlap_intervals(o);
    if (noi == 1) return;

    Interval fi= interval().first_overlap_interval(o);
    Sign fsn= fi.apply_to_endpoint(sign_at(), Interval::UPPER);

    if (fsn == ZERO){
      set_is_rational(fi.upper_endpoint_interval());
    } else if (contains_root(fi, lower_sign(), fsn)){
      set_interval(fi);
    } else {
      Interval si= interval().second_overlap_interval(o);
      if (noi==2){
	set_interval(si);
      } else {
	Sign ssn= si.apply_to_endpoint(sign_at(), Interval::UPPER);
	if (ssn== ZERO){
	  set_is_rational(si.upper_endpoint_interval());
	} else if (contains_root(si, fsn, ssn)) {
	  set_interval(si);
	} else {
	  set_interval(interval().third_overlap_interval(o));
	}
      }
    }
  }

  Comparison_result compare(const This &o) const {
    audit();
    o.audit();
    CGAL_Polynomial_precondition(-SMALLER == LARGER);
    // Elimate all cases where the functions are equal or negations
    if (is_normal() && o.is_normal()) {
      Order rel= interval().order(o.interval());
      if (rel == STRICTLY_BELOW) return SMALLER;
      else if (rel == STRICTLY_ABOVE) return LARGER;
    } else if (is_infinite() || o.is_infinite()){
      if (type() == o.type()) return EQUAL;
      else if (is_pos_inf() || o.is_neg_inf()) return LARGER;
      else if (is_neg_inf() || o.is_pos_inf()) return SMALLER;
      else {
	CGAL_Polynomial_assertion(0); return EQUAL;
      }
    }

    // The functions are now not equal and both roots are finite
    return compare_finite(o);
  }

  Comparison_result compare_finite(const This &o) const {
    audit();
    o.audit();
    if (is_rational() || o.is_rational()) {
      return compare_rational(o);
    } else {
      //Polynomial_assertion(function_==function_ && o.function_==o.function_);


      return compare_normal(o);
    }
  }
  
 
  Comparison_result compare_rational(const This &o) const {
    if (is_rational() && o.is_rational()){
      if (interval()== o.interval()) return EQUAL;
      else if (interval() < o.interval()) return SMALLER;
      else return LARGER;
    } else if (o.is_rational()) {
      return compare_with_rational(o);
    } else {
      return Comparison_result(-o.compare_with_rational(*this));
    }
  }

  
  //! Compare where this has type CONST
  Comparison_result compare_with_rational(const This &o) const {
    //std::cout << "Comparing " << *this << " and " << o << std::endl;
    CGAL_Polynomial_assertion(!is_rational());
    CGAL_Polynomial_assertion(o.is_rational());

    Order rel= interval().order(o.interval());
    if (rel == STRICTLY_BELOW) return SMALLER;
    else if (rel== STRICTLY_ABOVE) return LARGER;
    else {
      typename Traits::Sign_at osa= sign_at();
      Sign sn= o.interval().apply_to_endpoint(osa, Interval::LOWER);
      //std::cout << "Sign is " << sn << " and type of o is " << o.type_ << std::endl;
      if (sn==ZERO) {
	set_is_rational(o.interval());
	audit();
	return EQUAL;
      } else if (sn==POSITIVE){
	if (is_up()) return SMALLER;
	else return LARGER;
      } else {
	if (is_up()) return LARGER;
	else return SMALLER;
      }
    }
  }

  //! Compare two different intervals (where the endpoints are not the same, but they don't overlap).
  /*!
    This compares two when neither is const or inf, when this is to the left of o.
  */
  Comparison_result compare_normal(const This &o) const {
    CGAL_Polynomial_assertion(is_normal() && o.is_normal());

    refine_using(o.interval());
    o.refine_using(interval());

    do {
      audit(); o.audit();
     
       // See if that refinement changed anything
      if (is_rational() || o.is_rational()){
	return compare_rational(o);
      }
      Order ord= interval().order(o.interval());
      if (ord== STRICTLY_BELOW) return SMALLER;
      else if (ord == STRICTLY_ABOVE) return LARGER;
     
      // They overlap so if the functions are the same they must be equal
      typename Traits::Are_negations an= kernel_.are_negations_object();
      if (function_== o.function_ || an(function_, o.function_)){
	//std::cout << "Comparing the same function with " << *this << " and " << o << std::endl;
	return EQUAL;
      }

      CGAL_Polynomial_assertion(interval()==o.interval());
      
      refine();
      o.refine();
    } while (interval().approximate_width() > .0000001);

    std::cout << "Using sturm to compare " << *this << " and " << o << std::endl;
    /*Polynomial_assertion_code(typename Traits::Compare_isolated_roots_in_interval apred
			      = kernel_.compare_isolated_roots_in_interval_object(function_, function_));
    Polynomial_assertion_code(Comparison_result cc= interval().apply_to_interval(apred));
    Polynomial_assertion(cc == EQUAL);*/
    typename Traits::Compare_isolated_roots_in_interval pred=kernel_.compare_isolated_roots_in_interval_object(function_,o.function_);
    Comparison_result co= interval().apply_to_interval(pred);
    std::cout << "The result is " << co << std::endl;
    return co;
  }
  

  //! Check that everything is correct
  void audit() const {
    CGAL_Polynomial_assertion(!is_null());
#ifndef NDEBUG
    bool problem=false;
    if (is_infinite()){
      problem= problem || is_even_multiplicity();
    } else if (is_rational()){
      //problem = problem || (interval().apply_to_endpoint(sign_at(), Interval::LOWER) != ZERO);
    } else if (is_even_multiplicity()){
      if (is_up()){
	problem = problem || interval().is_singular();
	problem = problem || (interval().apply_to_endpoint(sign_at(), Interval::LOWER) != POSITIVE);
	problem = problem || (interval().apply_to_endpoint(sign_at(), Interval::UPPER) != POSITIVE);
      } else {
	problem = problem || interval().is_singular();
	problem = problem || (interval().apply_to_endpoint(sign_at(), Interval::LOWER) != NEGATIVE);
	problem = problem || (interval().apply_to_endpoint(sign_at(), Interval::UPPER) != NEGATIVE);
      }
    } else {
      if (is_up()){
	problem = problem || interval().is_singular();
	problem = problem || (interval().apply_to_endpoint(sign_at(), Interval::LOWER) != NEGATIVE);
	problem = problem || (interval().apply_to_endpoint(sign_at(), Interval::UPPER) != POSITIVE);
      } else {
	problem = problem || interval().is_singular();
	problem = problem || (interval().apply_to_endpoint(sign_at(), Interval::LOWER) != POSITIVE);
	problem = problem || (interval().apply_to_endpoint(sign_at(), Interval::UPPER) != NEGATIVE);
      }
    }

    if (problem){
      std::cerr << "Problem with interval.\n";
      std::cerr << "Type is " << int(type_ &1) << int(type_&2)<<int(type_&4) << int(type_&8) << std::endl;
      std::cerr << interval()<< ", " << interval().apply_to_endpoint(sign_at(), Interval::LOWER)
		<< interval().apply_to_endpoint(sign_at(), Interval::UPPER) << std::endl;
      CGAL_Polynomial_exactness_assertion(0);
    }
#endif
  }



  //! Is this used? 
  /*Simple_interval_root(Interval ii, Type type, Polynomial sign, bool mult): interval()(ii), type_(type),
									    function_(sign), is_odd_(mult){
    bool is_this_used;
    negated_=false;
    audit();
    compute_approximation();
    bool is_this_used;
    }*/


  typename Traits::Sign_at sign_at() const {
    // we loose the function when we find a rational value
    CGAL_Polynomial_precondition(!is_rational());
    return kernel_.sign_at_object(function_);
  }


  //! The representation of negative infinity
  /*static This negative_infinity(){
    This ret(INF);
    //ret.set_type(INF);
    //ret.compute_approximation();
    return ret;
    }*/
  
  bool is_up() const {
    return type_&UP;
  }
  bool is_down() const {
    return !(type_&UP);
  }
  
  

  void set_interval(const Interval &ii) const {
    ii_=ii;
  }

  static double double_inf_rep() {
    if (std::numeric_limits<double>::has_infinity){
      return (std::numeric_limits<double>::infinity());
    } else return (std::numeric_limits<double>::max());
  }

  Sign lower_sign() const {
    CGAL_Polynomial_precondition(!is_infinite());
    if (is_even_multiplicity()){
      if (is_up()) return POSITIVE;
      else return NEGATIVE;
    } else {
      if (is_up()) return NEGATIVE;
      else return POSITIVE;
    }
  }

  /*Type negate(const Type &t) const {
    Polynomial_assertion(t != POS_INF && t != NEG_INF);
    switch(t){
    case UP:
    return DOWN;
    case DOWN:
    return UP;
    case POS_INF:
    return NEG_INF;
    case NEG_INF:
    return POS_INF;
    case CONST:
    return CONST;
    case EVEN_UP:
    return EVEN_DOWN;
    case EVEN_DOWN:
    return EVEN_UP;
    default:
    Polynomial_assertion(0);
    return CONST;
    }
    }*/


  //! comptute a value that can be inspected in the compiler
  void compute_approximation(){
#ifndef NDEBUG
    approximation_=immutable_double_approximation(0.0001);
#endif
  }

  double immutable_double_approximation(double accuracy=0.00001) const {
    CGAL_Polynomial_expensive_precondition(!is_null());
    This temp = *this;
    return temp.compute_double(accuracy);
  }

  //! return true if the this is uninitialized
  bool is_null() const {
    return (type_&EVEN && type_&INF);
  }

  mutable Interval ii_;
  //Function f_;
  mutable Type type_;
  mutable Polynomial function_;
  Traits kernel_;
#ifndef NDEBUG
  double approximation_;
#endif
};

/*
template <class F>
double to_double(const Simple_interval_root<F> &f){
  return f.to_double();
}

template <class F>
std::pair<double,double> to_interval(const Simple_interval_root<F> &f){
  return f.to_interval();
  }*/

template <class F>
std::ostream &operator<<(std::ostream &out, const Simple_interval_root<F> &f){
  f.write(out);
  return out;
}

template <class F>
bool root_is_even_multiplicity(const Simple_interval_root<F> &f){
  return f.is_even_multiplicity();
}

template <class F>
typename F::NT to_rational(const Simple_interval_root<F> &f){
  return f.to_rational();
}

template <class F>
bool is_rational(const Simple_interval_root<F> &f){
  return f.is_rational();
}


CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE


//CGAL_BEGIN_NAMESPACE
namespace CGAL {
template <class F>
double to_double(const CGAL_POLYNOMIAL_NS::internal::Simple_interval_root<F> &f){
  //bool to_double_in_sir;
  return f.double_approximation();
}

template <class F>
std::pair<double, double> to_interval(const CGAL_POLYNOMIAL_NS::internal::Simple_interval_root<F> &f){
  //bool to_interval_in_sir;
  return f.double_interval();
}

}
//CGAL_END_NAMESPACE

namespace std {
  template <class Tr>
  struct numeric_limits<CGAL_POLYNOMIAL_NS::internal::Simple_interval_root<Tr> > {
    typedef CGAL_POLYNOMIAL_NS::internal::Simple_interval_root<Tr> T;
    static const bool is_specialized = true;
    static T min() throw () {return -T::infinity();}
    static T max() throw () {return T::infinity();}
    static const int digits =0;
    static const int digits10 =0;
    static const bool is_signed = true;
    static const bool is_integer = false;
    static const bool is_exact = true;
    static const int radix =0;
    static T epsilon() throw(){return T(0);}
    static T round_error() throw(){return T(0);}
    static const int min_exponent=0;
    static const int min_exponent10=0;
    static const int max_exponent=0;
    static const int max_exponent10=0;
    static const bool has_infinity=true;
    static const bool has_quiet_NaN = false;
    static const bool has_signaling_NaN= false;
    static const float_denorm_style has_denorm= denorm_absent;
    static const bool has_denorm_loss = false;
    static T infinity() throw() {return T::infinity();}
    static T quiet_NaN() throw(){return T(0);}
    static T denorm_min() throw() {return T(0);}
    static const bool is_iec559=false;
    static const bool is_bounded =false;
    static const bool is_modulo= false;
    static const bool traps = false;
    static const bool tinyness_before =false;
    static const float_round_style round_stype = round_toward_zero;
  };
};

#endif
