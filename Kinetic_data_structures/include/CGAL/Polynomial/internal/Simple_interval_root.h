// Copyright (c) 2005  Stanford University (USA).
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_SIMPLE_INTERVAL_ROOT_H
#define CGAL_POLYNOMIAL_SIMPLE_INTERVAL_ROOT_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Real_embeddable_traits.h>
#include <vector>
namespace CGAL { namespace POLYNOMIAL { namespace internal {

//! A root represented as a bounding interval and a polynomial.
/*!
  Representing an interval which contains one root of a function.
  \todo cache sturm sequence
*/
template <class Traits>
class Simple_interval_root
{
  typedef Simple_interval_root<Traits> This;


  typedef enum Fields {INVALID= 4,UP=1, INF=2} Type_fields;
  typedef unsigned char Type;
  typedef typename Traits::Function Polynomial;
  typedef typename Traits::FT NT;
public:
  Simple_interval_root(){
    set_type(INVALID);
    CGAL_Polynomial_assertion(is_null());
  }

  Simple_interval_root(double nt): function_(0), interval_(nt,nt) {
    if (std::numeric_limits<double>::has_infinity && (nt == std::numeric_limits<double>::infinity()
						   || -nt == std::numeric_limits<double>::infinity())) {
      if (nt == std::numeric_limits<double>::infinity() ) {
	set_type(INF|UP);
      }
      else {
	set_type(INF);
      }
    } else {
      set_type(0);
      ii_= std::make_pair(NT(nt),NT(nt));
    }
    audit();
    compute_approximation();
  }

  Simple_interval_root(const Simple_interval_root<Traits> &o):ii_(o.ii_), type_(o.type_), 
							      function_(o.function_),
							      kernel_(o.kernel_), interval_(o.interval_)
#ifndef NDEBUG
							     , approximation_(o.approximation_)
#endif
  {
    
  }

  /*Simple_interval_root(Type t): type_(t) {
    audit();
    compute_approximation();
    };*/

  //! represent a rational root
  Simple_interval_root(const NT &nt): ii_(std::make_pair(nt,nt)), function_(0), interval_(0,-1) {
    set_type(0);
    audit();
    compute_approximation();
  }


  //! Represent a root by an interval and a polynomial
  Simple_interval_root(const std::pair<NT,NT> &ii,
		       const Polynomial &sa,
		       Sign slb, Sign,
		       Traits k): ii_(ii), function_(sa),
				  kernel_(k), interval_(0,-1){
    //CGAL_Polynomial_precondition(!ii_.is_singular());
    //Sign slb= sign_(ii_.lower_bound());
    if (slb == CGAL::NEGATIVE) {
      set_type(UP);
    } else {
      set_type(0);
    }
    compute_approximation();
    audit();
  }

  static This infinity() {
    static This ret(std::numeric_limits<double>::infinity());
    CGAL_Polynomial_postcondition(ret.is_infinite());
    //ret.compute_approximation();
    return ret;
  }

  bool operator<(const This &o) const
  {
    CGAL_Polynomial_expensive_precondition(!is_null() && !o.is_null());
    Comparison_result r= compare(o);
    audit(); o.audit();
    return r==SMALLER;
  }

  bool operator>(const This &o) const
  {
    CGAL_Polynomial_expensive_precondition(!is_null() && !o.is_null());
    Comparison_result r= compare(o);
    audit(); o.audit();
    return r==LARGER;
  }

  bool operator!=(const This &o) const
  {
    if (is_null()) return !o.is_null();
    else if (o.is_null()) return true;
    CGAL_Polynomial_expensive_precondition(!is_null() && !o.is_null());
    Comparison_result r= compare(o);
    audit(); o.audit();
    return r != EQUAL;
  }
  bool operator<=(const This &o) const
  {
    CGAL_Polynomial_expensive_precondition(!is_null() && !o.is_null());
    Comparison_result r= compare(o);
    audit(); o.audit();
    return r != LARGER;
  }
  bool operator>=(const This &o) const
  {
    CGAL_Polynomial_expensive_precondition(!is_null() && !o.is_null());
    Comparison_result r= compare(o);
    audit(); o.audit();
    return r != SMALLER;
  }

  bool operator==(const This &o) const
  {
    if (is_null()) return o.is_null();
    else if (o.is_null()) return false;
    CGAL_Polynomial_expensive_precondition(!is_null() && !o.is_null());
    Comparison_result r= compare(o);
    audit(); o.audit();
    return r==EQUAL;
  }

 
  const std::pair<NT, NT>& isolating_interval() const
  {
    CGAL_Polynomial_precondition(!is_infinite());
    CGAL_Polynomial_expensive_precondition(!is_null());
    return ii_;
  }

  //! Write lots of info about the interval
  void write(std::ostream &o) const
  {
#ifndef NDEBUG
    This t=*this;
    t.write_internal(o);
#else
    write_internal(o);
#endif
  }

  void print() const
  {
    write(std::cout);
  }

  //! Negate the interval.
  This operator-() const
  {
    CGAL_Polynomial_expensive_precondition(!is_null());
    if (is_pos_inf()) {
      This ret;
      ret.type_=INF;
      return ret;
    } else if (is_neg_inf()) return infinity();
    else {
      This copy= *this;
      copy.ii_= std::make_pair(-ii_.second, -ii_.first);
      typename Traits::Negate_variable nf= kernel_.negate_variable_object();
      copy.function_= nf(function_);
      if (type_&UP) {
	copy.type_= 0;
      } else {
	copy.type_=UP;
      }
      return copy;
    }
  }

  //! Return true if the root is +/- infinity.
  bool is_infinite() const
  {
    return ((type_&INF)!= 0);
  }
 
  
  Comparison_result compare(const This &o) const
  {
    audit();
    o.audit();
    CGAL_Polynomial_precondition(-SMALLER == LARGER);
    if (is_infinite() || o.is_infinite()) {
      if (is_pos_inf() && o.is_pos_inf()) return EQUAL;
      if (is_neg_inf() && o.is_neg_inf()) return EQUAL;
      else if (is_pos_inf() || o.is_neg_inf()) return LARGER;
      else if (is_neg_inf() || o.is_pos_inf()) return SMALLER;
      else {
	CGAL_Polynomial_assertion(0); return EQUAL;
      }
    } else {
      CGAL::Comparison_result cmp;
      if (try_compare(o, cmp)) return cmp;
      else return compare_finite(o);
    }
  }
 std::pair<double, double> compute_interval(double accuracy=.0001) const
  {
    if (interval_.first > interval_.second) {
      //std::cout << "Computing interval for ";
      //write_raw(std::cout) << std::endl;
      internal_compute_interval(accuracy);
      /*std::cout << "Got " << CGAL::to_interval(ii_.first).first << "..." 
	<< CGAL::to_interval(ii_.second).second << std::endl;*/
      interval_=std::make_pair(CGAL::to_interval(ii_.first).first, CGAL::to_interval(ii_.second).second);
    } 
    return interval_;
  }

  double compute_double(double accuracy) const
  {
    if (is_infinite()) {
      if (is_up()) {
	return double_inf_rep();
      }
      else {
	return -double_inf_rep();
      }
    }
    //This t= *this;
    std::pair<double, double> i= compute_interval(accuracy);
    return (i.first+i.second)/2.0;
  }


  NT rational_between(const This &o) const {
    CGAL_precondition(CGAL::compare(*this,o) == CGAL::SMALLER);
    CGAL_precondition(ii_.first != ii_.second || ii_.second != o.ii_.first || o.ii_.first != o.ii_.second);
    compare(o);
    CGAL_precondition(ii_.first != ii_.second || ii_.second != o.ii_.first || o.ii_.first != o.ii_.second);
    //write_internal(std::cout) << std::endl;
    //write_internal(std::cout) << std::endl;
    if (is_neg_inf()) return o.ii_.first -1;
        
    // hopefully disjoint
    NT ret = ii_.second;
    NT step= NT(.0000000596046447753906250000000);
    
    do {
      while (This(ret) <= *this) {
	ret+= step;
      }
      while (This(ret) >= o) {
	ret -= step;
      }
      step= step/NT(2.0);
      /*std::cout << ret << "  (" << step << ")" 
		<< o.ii_.first - ii_.second 
		<< " " << o.ii_.second- ii_.first << std::endl;*/
    } while (This(ret) >= o || This(ret) <= *this);
    return ret;
  }

protected:
  std::pair<double, double> internal_compute_interval(double accuracy) const
  {
 
    if (type_&INF) {
      if (is_up()) {
	return std::pair<double, double>(double_inf_rep(), double_inf_rep());
      }
      else {
	return std::pair<double, double>(-double_inf_rep(), -double_inf_rep());
      }
    }

    //double oaw;                           //= ii_.approximate_width();
    // int ct=0;
    while (ii_.second != ii_.first && ii_.second-ii_.first > NT(accuracy)) {
      refine();
      /*++ct;
      if (ct== 30) {
	std::cerr << "Error subdividing ";
	write_internal(std::cerr) << std::endl;
	break;
	}*/
    }
    return std::make_pair(CGAL::to_interval(ii_.first).first, CGAL::to_interval(ii_.second).second);
  }
  std::ostream& write_internal(std::ostream &o) const
  {
    if (is_pos_inf()) {
      o << "inf";
    }
    else if (is_neg_inf()) {
      o << "-inf";
    }
    else {
      if (ii_.first == ii_.second) {
	o << ii_.first;
      } else {
	o << function_ << " in [" << ii_.first << "," << ii_.second << "]";
	o << " = " << internal_compute_interval(.00001).first 
	  << "..." << internal_compute_interval(.00001).second;
      }
    }
    return o;
  }

  std::ostream & write_raw(std::ostream &o) const
  {
    if (is_pos_inf()) {
      o << "inf";
    }
    else if (is_neg_inf()) {
      o << "-inf";
    }
    else {
      if (ii_.first == ii_.second) {
	o << ii_.first;
      } else {
	o << function_ << " in [" << ii_.first << "," << ii_.second << "]";
      }
    }
    return o;
  }

  void set_type(Type t)
  {
    type_=t;
  }

  bool is_pos_inf() const
  {
    return (type_&INF) && (type_&UP);
  }

  bool is_neg_inf() const
  {
    return (type_&INF) && !(type_&UP);
  }

  bool contains_root(Sign ls, Sign us) const
  {
    return (is_up() && ls==NEGATIVE && us==POSITIVE) || (is_down() && ls==POSITIVE && us==NEGATIVE);
  }

  void refine() const
  {
    if (ii_.first == ii_.second) return;
    //CGAL_Polynomial_precondition(!is_rational());
    //CGAL_Polynomial_precondition(is_normal());
    NT mp = (ii_.first + ii_.second)*NT(0.5);
    Sign sn= sign_at(mp);
    if (sn == ZERO) {
      ii_= std::make_pair(mp,mp);
    } else if (contains_root(lower_sign(), sn)) {
      ii_.second= mp;
    } else {
      ii_.first= mp;
    }
  }

  void refine_using(const std::pair<NT, NT> &o) const
  {
    if (ii_.first== ii_.second) return;
    std::vector<NT > plist;
    plist.push_back(ii_.first);
    if (o.first < ii_.second && o.first > ii_.first) {
      plist.push_back(o.first);
    }
    if (o.first != o.second && o.second < ii_.second && o.second > ii_.first) {
      plist.push_back(o.second);
    }
    plist.push_back(ii_.second);

    /*for (unsigned int i=0; i< plist.size(); ++i) {
      std::cout << plist[i] << "   ";
    }
    std::cout << std::endl;*/

    if (plist.size()==2) return;
    
    CGAL::Sign ps= lower_sign();
    //std::cout << "ps is " << ps << std::endl;
    CGAL_assertion(ps != CGAL::ZERO);
    for (unsigned int i=1; i< plist.size()-1; ++i) {
      CGAL::Sign sn= sign_at(plist[i]);
      //std::cout << "sn is " << ps << std::endl;
      if (sn==0) {
	ii_= std::make_pair(plist[i], plist[i]);
	audit();
	return;
      } else if (sn != ps) {
	ii_= std::make_pair(plist[i-1], plist[i]);
	audit();
	return;
      }
    }

    ii_= std::make_pair(plist[plist.size()-2], plist[plist.size()-1]);
    CGAL_postcondition(sign_at(plist[plist.size()-2]) == ps);
    CGAL_postcondition(sign_at(plist[plist.size()-1]) == -ps);
    audit();
  }

  Comparison_result compare_finite(const This &o) const
  {
    audit();
    o.audit();

    refine_using(o.ii_);
    o.refine_using(ii_);

    do {
      audit(); o.audit();

      CGAL::Comparison_result cmp;
      if (try_compare(o, cmp)) return cmp;
    
      refine();
      o.refine();
    } while (ii_.second-ii_.first  > NT(.0000001));

    //std::cout << "Using sturm to compare " << *this << " and " << o << std::endl;

    typename Traits::Compare_isolated_roots_in_interval pred=kernel_.compare_isolated_roots_in_interval_object(function_,o.function_);
    Comparison_result co= pred(ii_.first, ii_.second);
    //std::cout << "The result is " << co << std::endl;
    return co;
  }

  //! Check that everything is correct
  void audit() const
  {
    CGAL_Polynomial_assertion(!is_null());
#ifndef NDEBUG
    bool problem=false;
    if (type_&INF) {

    } else if (ii_.first == ii_.second) {
      
    } else {
      if (is_up()) {
	problem = problem || (sign_at(ii_.first) != NEGATIVE);
	problem = problem || (sign_at(ii_.second) != POSITIVE);
      } else {
	problem = problem || (sign_at(ii_.first) != POSITIVE);
	problem = problem || (sign_at(ii_.second) != NEGATIVE);
      }
    }
    if (problem) {
      std::cerr << "Problem with interval.\n";
      std::cerr << "Type is " << type_ << std::endl;
      std::cerr << ii_.first << "..." << ii_.second << ", " << sign_at(ii_.first) 
		<< sign_at(ii_.second) << std::endl;
      CGAL_Polynomial_exactness_assertion(0);
    }
#endif
  }



  CGAL::Sign sign_at(const NT &nt) const
  {
    return kernel_.sign_at_object()(function_, nt);
  }

  //! The representation of negative infinity
  /*static This negative_infinity(){
    This ret(INF);
    //ret.set_type(INF);
    //ret.compute_approximation();
    return ret;
    }*/

  bool try_compare(const This &o, CGAL::Comparison_result &cmp) const {
    if (ii_.first == ii_.second){
      if (o.ii_.first == o.ii_.second) {
	cmp= CGAL::compare(ii_.first, o.ii_.second);
	return true;
      } else {
	if (o.ii_.first < ii_.first && o.ii_.second > ii_.first) {
	  return false;
	} else {
	  cmp= CGAL::compare(ii_.first, o.ii_.first);
	  if (cmp == CGAL::EQUAL){
	    cmp=  CGAL::compare(ii_.first, o.ii_.second);
	  }
	  //CGAL_assertion(CGAL::compare(ii_.second, o.ii_.second) == cmp);
	  return true;
	}
      }
    } else if (o.ii_.first == o.ii_.second) {
      if (ii_.first < o.ii_.first && ii_.second > o.ii_.first) {
	  return false;
	} else {
	  cmp= CGAL::compare(ii_.first, o.ii_.first);
	  if (cmp == CGAL::EQUAL) {
	    cmp=  CGAL::compare(ii_.second, o.ii_.first);
	  }
	  //CGAL_assertion(CGAL::compare(ii_.second, o.ii_.second) == cmp);
	  return true;
	}
    } else {
      if (ii_.first >= o.ii_.second) {
	cmp= CGAL::LARGER;
	return true;
      } else if (ii_.second <= o.ii_.first) {
	cmp= CGAL::SMALLER;
	return true;
      } else return false;
    }
  }

  bool is_up() const
  {
    return type_&UP;
  }
  bool is_down() const
  {
    return !(type_&UP);
  }

  static double double_inf_rep() {
    if (std::numeric_limits<double>::has_infinity) {
      return (std::numeric_limits<double>::infinity());
    } else return ((std::numeric_limits<double>::max)());
  }

  Sign lower_sign() const
  {
    CGAL_Polynomial_precondition(!is_infinite());
    if (is_up()) return NEGATIVE;
    else return POSITIVE;
  }

  //! comptute a value that can be inspected in the compiler
  void compute_approximation() {
#ifndef NDEBUG
    approximation_=immutable_double_approximation(0.0001);
#endif
  }

  double immutable_double_approximation(double accuracy=0.00001) const
  {
    CGAL_Polynomial_expensive_precondition(!is_null());
    This temp = *this;
    return temp.internal_compute_interval(accuracy).first;
  }

  //! return true if the this is uninitialized
  bool is_null() const
  {
    return (type_&INVALID)!= 0;
  }

  mutable std::pair<NT, NT> ii_;
  //Function f_;
  Type type_;
  Polynomial function_;
  Traits kernel_;
  mutable std::pair<double,double> interval_;
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
  return f.to_ii_;
  }*/

template <class F>
std::ostream &operator<<(std::ostream &out, const Simple_interval_root<F> &f)
{
  f.write(out);
  return out;
}

/*
template <class F>
bool root_is_even_multiplicity(const Simple_interval_root<F> &f)
{
  return f.is_even_multiplicity();
}


template <class F>
typename F::NT to_rational(const Simple_interval_root<F> &f)
{
  return f.to_rational();
}


template <class F>
bool is_rational(const Simple_interval_root<F> &f)
{
  return f.is_rational();
  }*/


} } } //namespace CGAL::POLYNOMIAL::internal

namespace CGAL {


template <class T>
class Real_embeddable_traits< CGAL::POLYNOMIAL::internal::Simple_interval_root<T> > 
  : public INTERN_RET::Real_embeddable_traits_base< CGAL::POLYNOMIAL::internal::Simple_interval_root<T> , Tag_true > {
public:
  typedef CGAL::POLYNOMIAL::internal::Simple_interval_root<T>  Type;
  class Abs 
    : public std::unary_function< Type, Type > {
  public:
    Type operator()( const Type& x ) const {
      if (x < Type(0)) return -x;
      else return x;
    }
  };
    
  class Sgn 
    : public std::unary_function< Type, ::CGAL::Sign > {
  public:
    ::CGAL::Sign operator()( const Type& x ) const {
      return static_cast<CGAL::Sign>(x.compare(0));
    }        
  };
    
  class Compare 
    : public std::binary_function< Type, Type,
			      Comparison_result > {
  public:
    Comparison_result operator()( const Type& x, 
				  const Type& y ) const {
      return x.compare(y);
    }
        
    CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( Type,
							 Comparison_result )
        
      };
    
  class To_double 
    : public std::unary_function< Type, double > {
  public:
    double operator()( const Type& x ) const {
      // this call is required to get reasonable values for the double
      // approximation
      return x.compute_double(.00000001);
    }
  };
    
  class To_interval 
    : public std::unary_function< Type, std::pair< double, double > > {
  public:
    std::pair<double, double> operator()( const Type& x ) const {

      return x.compute_interval(.00001);
    }          
  };
};



} //namespace CGAL

namespace std
{
  template <class Tr>
  class numeric_limits<CGAL_POLYNOMIAL_NS::internal::Simple_interval_root<Tr> >
  {
  public:
    typedef CGAL_POLYNOMIAL_NS::internal::Simple_interval_root<Tr> T;
    static const bool is_specialized = true;
    static T min BOOST_PREVENT_MACRO_SUBSTITUTION () throw () {return -T::infinity();}
    static T max BOOST_PREVENT_MACRO_SUBSTITUTION () throw () {return T::infinity();}
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
    static const bool has_quiet_NaN = true;
    static const bool has_signaling_NaN= false;
    static const float_denorm_style has_denorm= denorm_absent;
    static const bool has_denorm_loss = false;
    static T infinity() throw() {return T::infinity();}
    static T quiet_NaN() throw(){return T();}
    static T denorm_min() throw() {return T(0);}
    static const bool is_iec559=false;
    static const bool is_bounded =false;
    static const bool is_modulo= false;
    static const bool traps = false;
    static const bool tinyness_before =false;
    static const float_round_style round_stype = round_toward_zero;
  };
}
#endif
