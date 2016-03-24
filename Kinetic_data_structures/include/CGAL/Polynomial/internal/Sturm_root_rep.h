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

#ifndef CGAL_POLYNOMIAL_STURM_ROOT_REP_H
#define CGAL_POLYNOMIAL_STURM_ROOT_REP_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Rational/Sign_Sturm_sequence.h>
#include <limits>
#include <cfloat>



namespace CGAL { namespace POLYNOMIAL { namespace internal {

//==================
// the Root class
//==================
template<class Solver_t, class Interval_t>
class Sturm_root_rep
{
public:
  typedef Solver_t                            Solver;
  typedef Interval_t                          Interval;
  typedef typename Solver::Traits             Traits;

  //  typedef typename Solver::Exact_nt           Exact_nt;
  //  typedef typename Solver::Storage_function   Storage_function;
  //  typedef typename Solver::Exact_function     Exact_function;
  typedef typename Solver::Method_tag         Method_tag;

  //  typedef typename Solver::Function_handle    Function_handle;
  typedef typename Solver::Standard_sequence  Standard_sequence;

  typedef typename Solver::Polynomial         Polynomial;
  typedef typename Solver::NT                 NT;
  typedef typename Solver::Sign_at            Sign_at;
  typedef NT                                  Exact_nt;

  typedef Sturm_root_rep<Solver,Interval>  Self;

  class Sign_at_functor
  {
  public:
    typedef Sign result_type;

    Sign_at_functor(const Self* outer) : outer(outer) {}

    template<class T>
    result_type operator()(const T& x) const
    {
      Sign s1 = outer->sseq.sign_at(x, 0);
      Sign s2 = outer->sseq.sign_at_gcd(x);

      CGAL_assertion( s1 == CGAL::ZERO || s2 != CGAL::ZERO );
      return s1 * s2;
    }

  private:
    const Self* outer;
  };

  typedef typename Traits::Sign_Sturm_sequence   Sign_Sturm_sequence;

protected:
  typedef std::pair<Interval,Interval>               Interval_pair;

public:
  //==========================================
  // METHODS FOR COMPUTING THE MULTIPLICITY
  //==========================================

  bool is_rational() const
  {
    return is_exact();
  }

  NT to_rational() const
  {
    CGAL_Polynomial_precondition( is_rational() );
    return ivl.lower_bound();
  }

  bool is_exact() const
  {
    return ivl.is_singular();
  }

  // computes the multiplicity using the interval number type;
  // if computation was successful the method sets success to true;
  // if computation was not possible the method sets success to
  // false;

  void compute_multiplicity() const
  {
    CGAL_precondition( multiplicity_ == 0 );

    if ( p_.degree() == 1 ) {
      multiplicity_ = 1;
      return;
    }

    if ( is_exact() ) {
      compute_multiplicity_exact();
    }
    else {
      compute_multiplicity_interval();
    }
  }

  void compute_multiplicity_exact() const
  {
    Polynomial q = p_;

    while ( true ) {
      multiplicity_++;
      q = tr_.differentiate_object()(q);

      if ( q.is_zero() ) { break; }

      Sign_at  sign_at_q(q);
      Sign s = ivl.apply_to_endpoint(sign_at_q, Interval::LOWER);
      if ( s != CGAL::ZERO ) {
	break;
      }
    }
  }

  void compute_multiplicity_interval() const
  {
    bool is_even_multiplicity_ = is_even_multiplicity();

    typename Traits::Differentiate differentiate =
      tr_.differentiate_object();

    Polynomial q = p_;
    bool first_time = true;
    while ( true ) {
      if ( first_time ) {
	first_time = false;
	if ( !is_even_multiplicity_ ) {
	  multiplicity_++;
	  q = differentiate(q);
	}
	else {
	  multiplicity_ += 2;
	  q = differentiate( differentiate(q) );
	}
      }
      else {
	multiplicity_ += 2;
	q = differentiate( differentiate(q) );
      }

      if ( q.degree() <= 0 ) { break; }

      Sign_Sturm_sequence sign_sturm =
	tr_.sign_Sturm_sequence_object(p_, q);

      int sign = ivl.apply_to_interval(sign_sturm);

      if ( sign != 0 ) { break; }
    }
  }

  //===========================================
  // HELPER METHODS FOR UPDATING THE VALUES
  // ASSOCIATED WITH THE INTERVAL
  //===========================================

  void set_lower(const typename Interval::NT& l,
		 const Sign& s_l) const
  {
    ivl.set_lower(l);
    s_lower = s_l;
  }

  void set_upper(const typename Interval::NT& u,
		 const Sign& s_u) const
  {
    ivl.set_upper(u);
    s_upper = s_u;
  }

  void set_interval(const typename Interval::NT& x) const
  {
    ivl.set_lower(x);
    ivl.set_upper(x);
    s_lower = s_upper = CGAL::ZERO;
  }



  template <class This>
  Comparison_result
  compare_finite(const This &r, bool subdiv=true) const {
    // consider the cases that the root is known exactly;
    // this is equivalent to saying that the two endpoints for the
    // interval containing the root are the same;
    // moreover if the two interval endpoints are not the same we
    // make sure that the interval endpoints are not roots.
    // this will make life easier afterwards.

    if ( is_exact() ) {
      //	std::cout << "first is exact" << std::endl;
      if ( r.is_exact() ) {
	//	  std::cout << "second is exact" << std::endl;
	return CGAL::compare( lower_bound(), r.lower_bound() );
      }
      else {
	if ( upper_bound() <= r.lower_bound() ) {
	  return CGAL::SMALLER;
	}
	else if ( lower_bound() >= r.upper_bound() ) {
	  return CGAL::LARGER;
	} else if ( upper_bound() > r.lower_bound() &&
                    lower_bound() < r.upper_bound() ) {
	  Sign s_r_lb = r.sign_lower();
	  {
	    Sign s_r_ub = r.sign_upper();
	    CGAL_assertion( s_r_lb != CGAL::ZERO && s_r_ub != CGAL::ZERO );
	    if(0) s_r_ub= CGAL::ZERO;
	  }

	  Sign s_at_r = r.sign_at( *this, Interval::LOWER );

	  if ( s_at_r == CGAL::ZERO ) { return CGAL::EQUAL; }
	  return ( s_at_r == s_r_lb ) ? CGAL::SMALLER : CGAL::LARGER;
	}
      }
    }
    else if ( r.is_exact() ) {

      // we have already checked the case that this root is also
      // known exactly
      //	std::cout << "second is exact" << std::endl;
      if ( upper_bound() <= r.lower_bound() ) {
	//	  std::cout << "this.upper bound <= other.lower bound" << std::endl;
	return CGAL::SMALLER;
      }
      else if ( lower_bound() >= r.upper_bound() ) {
	//	  std::cout << "this.lower bound >= other.upper bound" << std::endl;
	return CGAL::LARGER;
      } else if (  upper_bound() > r.lower_bound() &&
		   lower_bound() < r.upper_bound() ) {
	//	   std::cout << "other in interval of this" << std::endl;
	Sign s_ub = sign_upper();
	{
	  Sign s_lb = sign_lower();
	  CGAL_assertion( s_lb != CGAL::ZERO && s_ub != CGAL::ZERO );
	  if(0) s_lb= CGAL::ZERO;
	}

	Sign s_at_this = sign_at( r, Interval::LOWER );

	if ( s_at_this == CGAL::ZERO ) { return CGAL::EQUAL; }
	return ( s_ub == s_at_this ) ? CGAL::SMALLER : CGAL::LARGER;
      }
    }
    else {
      // now the roots are in the interior of interval
      // check if the interiors of the intervals are disjoint
      //	std::cout << "both not exact" << std::endl;
      //	int prec = std::cout.precision();
      //	std::cout.precision(16);
      //	std::cout << "this:  " << lower_bound() << " " << upper_bound()
      //		  << std::endl;
      //	std::cout << "other: " << r.lower_bound() << " "
      //		  << r.upper_bound() << std::endl;
      //	std::cout.precision(prec);

      if ( upper_bound() <= r.lower_bound() ) {
	//	  std::cout << "this.upper bound <= other.lower bound" << std::endl;
	return CGAL::SMALLER;
      }
      else if ( lower_bound() >= r.upper_bound() ) {
	//	  std::cout << "this.lower bound >= other.upper bound" << std::endl;
	return CGAL::LARGER;
      } else {
	if (subdiv) {
	  int count=0;
	  while ( upper_bound() > r.lower_bound() &&
		  lower_bound() < r.upper_bound() ) {
	    subdivide();
	    ++count;
	    if (count==4) break;
	    //	  std::cout << "intersecting intervals" << std::endl;
	  }
	  return compare_finite(r, false);
	} else {
	  return compare_intersecting(r);
	}
      }
    }

    bool this_line_should_not_have_been_reached = false;
    CGAL_assertion( this_line_should_not_have_been_reached );
    if (0) this_line_should_not_have_been_reached= false;
    return CGAL::EQUAL;
  }
  

  //===========================================
  // HELPER METHODS FOR COMPARISON OPERATORS
  //===========================================
  template<class This>
  Comparison_result compare(const This& r) const
  {
    // check against positive and negative infinity
    if (idx == -2) {
      if ( r.idx == -2 ) { return CGAL::EQUAL; }
      return CGAL::SMALLER;
    }
    if (idx == -1) {
      if ( r.idx == -1 ) { return CGAL::EQUAL; }
      return CGAL::LARGER;
    }
    if (r.idx == -2) {
      // I have already checked if both are negative infinity
      return CGAL::LARGER;
    }
    if (r.idx == -1) {
      // I have already checked if both are positive infinity
      return CGAL::SMALLER;
    }
    return compare_finite(r);
  }

  
  //--------------------------------------------------------------


  template<class This>
  Comparison_result
  compare_intersecting(const This& r) const
  {
    // in this method we know that the intervals are not
    // trivial, are intersecting and that the endpoints are not
    // roots

    // Case 1: the intervals have common left endpoint
    if ( lower_bound() == r.lower_bound() ) {
      if ( upper_bound() > r.upper_bound() ) {
	Sign s = sign_at( r, Interval::UPPER );

	if ( s == CGAL::ZERO ) {
	  // the root of this is equal to the upper bound of r
	  //	    set_interval( r, Interval::UPPER );
	  set_interval( r.upper_bound() );
	  return CGAL::LARGER;
	}
	else {
	  if ( s == sign_lower() ) {
	    return CGAL::LARGER;
	  }
	  else {
	    //	      set_upper( r, Interval::UPPER, s );
	    set_upper( r.upper_bound(), s );
	    return compare_same_interval(r);
	  }
	}
      }
      else if ( r.upper_bound() > upper_bound() ) {
	const This* new_this = static_cast<const This*>(this);
	return opposite(r.compare_intersecting( *new_this ));
      }
      else if ( upper_bound() == r.upper_bound() ) {
	return compare_same_interval(r);
      }
    }

    // Case 2: the intervals have common right endpoint
    if ( upper_bound() == r.upper_bound() ) {
      if ( lower_bound() < r.lower_bound() ) {
	Sign s = sign_at( r, Interval::LOWER );
	if ( s == CGAL::ZERO ) {
	  // the root of this is equal to the lower bound of r
	  //	    set_interval( r, Interval::LOWER );
	  set_interval( r.lower_bound() );
	  return CGAL::SMALLER;
	}
	else {
	  if ( s == sign_upper() ) {
	    return CGAL::SMALLER;
	  }
	  else {
	    //	      set_lower( r, Interval::LOWER, s );
	    set_lower( r.lower_bound(), s );
	    return compare_same_interval(r);
	  }
	}
      }
      else if ( r.lower_bound() < lower_bound() ) {
	const This* new_this = static_cast<const This*>(this);
	return opposite(r.compare_intersecting( *new_this ));
      }
      else if ( lower_bound() == r.lower_bound() ) {
	return compare_same_interval(r);
      }
    }

    // Case 3: the upper bound of r is an interior point of the
    //         interval of this and the lower bound of r is before
    //         the lower bound of this
    if ( lower_bound() > r.lower_bound() &&
	 upper_bound() > r.upper_bound() ) {

      Sign s_this_at_r = r.sign_at( *this, Interval::LOWER );
      if ( s_this_at_r == CGAL::ZERO ) {
	//	  r.set_interval( *this, Interval::LOWER );
	r.set_interval( lower_bound() );
	return CGAL::LARGER;
      }

      Sign s_r_at_this = sign_at( r, Interval::UPPER );
      if ( s_r_at_this == CGAL::ZERO ) {
	//	  set_interval( r, Interval::UPPER );
	set_interval( r.upper_bound() );
	return CGAL::SMALLER;
      }

      if ( s_r_at_this != CGAL::ZERO && s_this_at_r != CGAL::ZERO ) {

	if ( s_this_at_r != r.sign_lower() ) {
	  //	    r.set_upper( *this, Interval::LOWER, s_this_at_r );
	  r.set_upper( lower_bound(), s_this_at_r );
	  return CGAL::LARGER;
	}

	if ( s_r_at_this != sign_upper() ) {
	  //	    set_lower( r, Interval::UPPER, s_r_at_this );
	  set_lower( r.upper_bound(), s_r_at_this );
	  return CGAL::LARGER;
	}

	//	  set_upper( r, Interval::UPPER, s_r_at_this );
	//	  r.set_lower( *this, Interval::LOWER, s_this_at_r );
	set_upper( r.upper_bound(), s_r_at_this );
	r.set_lower( lower_bound(), s_this_at_r );

	return compare_same_interval(r);
      }
    }

    // Case 4: the upper bound of this is an interior point of the
    //         interval of r and the lower bound of this is before
    //         the lower bound of r
    if ( lower_bound() < r.lower_bound() &&
	 upper_bound() < r.upper_bound() ) {
      const This* new_this = static_cast<const This*>(this);
      return opposite(r.compare_intersecting( *new_this ));
    }

    // Case 5: the interval of r is contained in the interval of
    //         this
    if ( lower_bound() < r.lower_bound() &&
	 upper_bound() > r.upper_bound() ) {

      Sign s_rl_at_this = sign_at( r, Interval::LOWER );
      if ( s_rl_at_this == CGAL::ZERO ) {
	//	  set_interval( r, Interval::LOWER );
	set_interval( r.lower_bound() );
	return CGAL::SMALLER;
      }

      Sign s_ru_at_this = sign_at( r, Interval::UPPER );
      if ( s_ru_at_this == CGAL::ZERO ) {
	//	  set_interval( r, Interval::UPPER );
	set_interval( r.upper_bound() );
	return CGAL::LARGER;
      }

      Sign s_l = sign_lower();
      if ( s_l != s_rl_at_this ) {
	//	  set_upper( r, Interval::LOWER, s_rl_at_this );
	set_upper( r.lower_bound(), s_rl_at_this );
	return CGAL::SMALLER;
      }

      Sign s_u = sign_upper();
      if ( s_u != s_ru_at_this ) {
	//	  set_lower( r, Interval::UPPER, s_ru_at_this );
	set_lower( r.upper_bound(), s_ru_at_this );
	return CGAL::LARGER;
      }

      if ( s_rl_at_this != CGAL::ZERO &&
	   s_ru_at_this != CGAL::ZERO   ) {
	//	  set_lower( r, Interval::LOWER, s_rl_at_this );
	//	  set_upper( r, Interval::UPPER, s_ru_at_this );
	set_lower( r.lower_bound(), s_rl_at_this );
	set_upper( r.upper_bound(), s_ru_at_this );

	return compare_same_interval(r);
      }
    }

    // Case 6: the interval of this is contained in
    //         the interval of r
    if ( lower_bound() > r.lower_bound() &&
	 upper_bound() < r.upper_bound() ) {
      const This* new_this = static_cast<const This*>(this);
      return opposite(r.compare_intersecting( *new_this ));
    }

    CGAL_postcondition_msg(0, "This line should not have been reached.\n");
    //bool this_line_should_not_have_been_reached = false;
    //CGAL_assertion( this_line_should_not_have_been_reached );

    return CGAL::EQUAL;
  }                                         // end-compare_intersecting

  //--------------------------------------------------------------
  void subdivide() const
  {
    if (!refined_) {
      refined_=true;
      //std::cout << "Refining " << *this << " " << std::endl;
    }


    Interval_pair ivl_pair = ivl.split();
    Interval first_half  = ivl_pair.first;
    Interval second_half = ivl_pair.second;

    Sign_at_functor sign_at_f(this);

    Sign s_mid =
      first_half.apply_to_endpoint( sign_at_f, Interval::UPPER );

    if ( s_mid == CGAL::ZERO ) {
      ivl = first_half.endpoint_interval( Interval::UPPER );
      return;
    }

    if ( s_mid == sign_upper() ) {
      ivl = first_half;
    }
    else {
      ivl = second_half;
    }
  }

  

  template<class Child>
  Polynomial compute_simple(Field_tag, const Child&) const
  {
    Polynomial gcd = sseq.exact( sseq.exact_size() - 1 );
    return p_ / gcd;
  }

  /*Polynomial compute_simple(Integral_domain_without_division_tag, const Self&) const
    {
    Polynomial gcd = sseq[sseq.size() - 1];
    return tr_.pseudo_quotient_object()(p_, gcd);
    }*/
  
  Polynomial compute_simple(Field_tag, const Self&) const
  {
    Polynomial gcd = sseq[sseq.size() - 1];
    return tr_.quotient_object()(p_, gcd);
  }

  template<class This>
  Comparison_result
  //    compare_same_interval(const Self& r) const
  compare_same_interval(const This& r) const
  {
    CGAL_precondition( lower_bound() == r.lower_bound() );
    CGAL_precondition( upper_bound() == r.upper_bound() );

    // here I have to do some subdivisions...

    Method_tag mtag;
    Sign_Sturm_sequence seq
      = tr_.sign_Sturm_sequence_object(r.compute_simple(mtag, r),
				       compute_simple(mtag, r) );

    int sign_of_r_at_this;

    sign_of_r_at_this =
      seq.sum_of_signs(lower_bound(), upper_bound());

    if ( sign_of_r_at_this == 0 ) { return CGAL::EQUAL; }

    Sign s_r;
    if ( sign_of_r_at_this < 0 ) {
      s_r = CGAL::NEGATIVE;
    }
    else {
      s_r = CGAL::POSITIVE;
    }

    Sign s_l = sign_lower();

    return ( s_r != s_l ) ? CGAL::SMALLER : CGAL::LARGER;
  }

  //===================
  // SIGNS EVALUATORS
  //===================
  Sign sign_lower() const
  {
    return s_lower;
  }

  Sign sign_upper() const
  {
    return s_upper;
  }

  template<class This>
  Sign sign_at(const This& r, typename Interval::Endpoint b) const
  {
    Sign_at_functor sign_at_f(this);

    return r.ivl.apply_to_endpoint(sign_at_f, b);
  }

public:
  //==================================
  // POSITIVE AND NEGATIVE INFINITY
  //==================================
  static Self infinity(const Traits& tr = Traits()) {
    return Self(-1, 0, tr);
  }

  //================
  // CONSTRUCTORS
  //================
protected:
  Sturm_root_rep(int i, int, const Traits& tr)
    : idx(i), ivl(), p_(), sseq(), s_lower(CGAL::ZERO),
      s_upper(CGAL::ZERO), multiplicity_(0), tr_(tr), refined_(false)  {}

public:
  Sturm_root_rep(const Traits& tr = Traits())
    : idx(-1), ivl(), p_(), sseq(), s_lower(CGAL::ZERO),
      s_upper(CGAL::ZERO), multiplicity_(0), tr_(tr), refined_(false)  {}
  //-------
  template<class T>
  Sturm_root_rep(const T& a, const Traits& tr = Traits())
    : idx(1), ivl(a), p_(), sseq(), s_lower(CGAL::ZERO),
      s_upper(CGAL::ZERO), multiplicity_(0), tr_(tr), refined_(false)  {}
  //-------
  Sturm_root_rep(const Interval& ivl,
		 const Polynomial& p,
		 const Standard_sequence& sseq, int idx,
		 const Traits& tr)
    : idx(idx), ivl(ivl), p_(p), sseq(sseq),
      s_lower(CGAL::ZERO), s_upper(CGAL::ZERO), multiplicity_(0),
      tr_(tr), refined_(false) {
    //++sturm_created__;
    Sign_at_functor sign_at_p(this);
    s_lower = apply(sign_at_p, ivl.lower_bound());
    s_upper = apply(sign_at_p, ivl.upper_bound());

    //      s_lower = sign_at( *this, Interval::LOWER );
    //      s_upper = sign_at( *this, Interval::UPPER );
  }
  //-------
  Sturm_root_rep(const Self& other)
    : idx(other.idx), ivl(), p_(other.p_), sseq(other.sseq),
      s_lower(other.s_lower),   s_upper(other.s_upper),
      multiplicity_(other.multiplicity_), tr_(other.tr_), refined_(other.refined_){
    if ( other.idx >= 1 ) {
      ivl = other.ivl;
    }
  }

  //=======================
  // ASSIGNMENT OPERATOR
  //=======================
  const Self& operator=(const Self& other) {
    if ( &other == this ) { return *this; }
    idx = other.idx;
    p_ = other.p_;
    sseq = other.sseq;
    s_lower = other.s_lower;
    s_upper = other.s_upper;
    multiplicity_ = other.multiplicity_;
    tr_ = other.tr_;

    if ( idx >= 1 ) {
      ivl = other.ivl;
    }
    return *this;
  }

  //=========================
  // COMPARISON OPERATORS
  //=========================
  bool operator!=(const Self& r) const
  {
    return compare(r) != CGAL::EQUAL;
  }

  bool operator==(const Self& r) const
  {
    return compare(r) == CGAL::EQUAL;
  }

  bool operator<=(const Self& r) const
  {
    return compare(r) != CGAL::LARGER;
  }

  bool operator>=(const Self& r) const
  {
    return compare(r) != CGAL::SMALLER;
  }

  bool operator>(const Self& r) const
  {
    return compare(r) == CGAL::LARGER;
  }

  bool operator<(const Self& r) const
  {
    return compare(r) == CGAL::SMALLER;
  }

  Self operator-() const
  {
    if ( idx == -1 ) {
      return Self(-2, 0, tr_);
    }
    else if ( idx == -2 ) {
      return Self(-1,0, tr_);
    }
    else if (idx == 1) {
      // HACK by Daniel
      return Self(-ivl, tr_);
    }
    else {

      Interval m_ivl = -ivl;
      Polynomial m_p = tr_.negate_variable_object()(p_);
      Standard_sequence m_sseq= tr_.standard_sequence_object(m_p);

      Self m_root(m_ivl, m_p, m_sseq, idx, tr_);

      m_root.multiplicity_ = multiplicity_;

      return m_root;
    }
  }

  //===========================
  // THE MULTIPLICITY METHOD
  //===========================
  unsigned int multiplicity() const
  {
    if ( idx < 0 ) { return 0; }
    if ( multiplicity_ == 0 ) {
      compute_multiplicity();
    }last_zero_=false;
    if (multiplicity_%2==0) return multiplicity_/2;
    else return multiplicity_;
  }

  //====================
  // THE PARITY METHOD
  //====================
  bool is_even_multiplicity() const
  {
    return false;
    /*if ( idx < 0 ) { return false; }

    if ( is_exact() ) {
      return (multiplicity() % 2 == 0);
    }
    else {
      Sign_at_functor sign_at_f(this);

      Sign sl = apply(sign_at_f, ivl.lower_bound());
      Sign su = apply(sign_at_f, ivl.upper_bound());
      return (sl == su);
      }*/
  }

  //==========================
  // ACCESS TO THE INTERVAL
  //==========================
  typename Interval::NT lower_bound() const { 
    CGAL_precondition(idx >= 0);
    return ivl.lower_bound(); }
  typename Interval::NT upper_bound() const { 
    CGAL_precondition(idx >= 0);
    return ivl.upper_bound(); }
  Interval    interval()    const { return ivl; }
  int         index()       const { return idx; }
  const Polynomial&     polynomial() const { return p_; }
  //    Polynomial  simple_polynomial() const { return p; }

public:
  //======================
  // CONVERTOR TO double
  //======================
  double compute_double(double acc = 1e-10) const
  {
    if (idx < 0) {
      double inf=std::numeric_limits<double>::has_infinity? std::numeric_limits<double>::infinity() : (std::numeric_limits<double>::max)();
      if ( idx == -1 ){
	return inf;
      } else return -inf;
    }
    /*if ( idx == -2 ) {
      return -DBL_MAX * DBL_MAX;
      }*/

    if ( is_exact() ) {
      return CGAL_POLYNOMIAL_TO_DOUBLE(lower_bound());
    }

    Exact_nt xacc(acc);


    while ( ivl.approximate_width() > acc ) {
      subdivide();
      /*if (!refined_) {
	++sturm_refined__;
	refined_=true;
	std::cout << "Refining " << *this << " to compute double " << std::endl;
	}*/
      if ( is_exact() ) { break; }
    }

    return CGAL_POLYNOMIAL_TO_DOUBLE( lower_bound() );
  }

  const std::pair<NT, NT>& isolating_interval() const {
    CGAL_precondition(idx >=0);
    return ivl.to_pair();
  }

  //========================
  // CONVERTOR TO interval
  //========================
  std::pair<double,double> compute_interval() const
  {
    if (*this == infinity()){ 
      return std::make_pair(std::numeric_limits<double>::infinity(),
			    std::numeric_limits<double>::infinity());
    }
    else if (*this == -infinity()) {
      return std::make_pair(-std::numeric_limits<double>::infinity(),
			    -std::numeric_limits<double>::infinity());
    } else {
      compute_double();
	    
      std::pair<double,double> i_low =
	CGAL_POLYNOMIAL_TO_INTERVAL(ivl.lower_bound());
      std::pair<double,double> i_high =
	CGAL_POLYNOMIAL_TO_INTERVAL(ivl.upper_bound());
	    
      return std::pair<double,double>(i_low.first, i_high.second);
    }
  }

  //================
  // STREAM WRITER
  //================
  template<class Stream>
  Stream& write(Stream& os) const
  {
    if ( idx == -2 ) {
      os << "-inf";
    }
    else if ( idx == -1 ) {
      os << "inf";
    }
    else {
      os << "{" << ivl << ", " << idx << "}";
    }
    if (idx != -2 && idx != -1) {
      Self copy = *this;
      os << " = " << copy.compute_double();
    }
    return os;
  }

protected:
  bool is_up_;
  int  
  int                    idx;
  mutable Interval       ivl;
  Polynomial             p_;
  Standard_sequence      sseq;
  mutable Sign           s_lower, s_upper;
  mutable unsigned int   multiplicity_;
  Traits                 tr_;
  mutable bool refined_;
};

template<class Stream, class S, class I>
Stream&
operator<<(Stream& os, const Sturm_root_rep<S,I>& r)
{
  return r.write(os);
}

/*template<class S, class I>
  std::pair<double,double>
  to_interval(const Sturm_root_rep<S,I>& r)
  {
  return r.compute_interval();
  }*/

} } } //namespace CGAL::POLYNOMIAL::internal

namespace CGAL {

/*template<class S, class I>
  double
  to_double(const CGAL_POLYNOMIAL_NS::internal::Sturm_root_rep<S,I>& r)
  {
  return r.compute_double();
  }


  template<class S, class I>
  std::pair<double,double>
  to_interval(const CGAL_POLYNOMIAL_NS::internal::Sturm_root_rep<S,I>& r)
  {
  return r.compute_interval();
  }*/


template <class T, class I>
class Real_embeddable_traits< CGAL::POLYNOMIAL::internal::Sturm_root_rep<T,I> > 
  : public INTERN_RET::Real_embeddable_traits_base< CGAL::POLYNOMIAL::internal::Sturm_root_rep<T,I> , Tag_true > {
public:
  typedef CGAL::POLYNOMIAL::internal::Sturm_root_rep<T,I>  Type;
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
      return x.compute_double();
    }
  };
    
  class To_interval 
    : public std::unary_function< Type, std::pair< double, double > > {
  public:
    std::pair<double, double> operator()( const Type& x ) const {

      return x.compute_interval();
    }          
  };
};


} //namespace CGAL

namespace std
{
  template <class S, class I>
  class numeric_limits<CGAL_POLYNOMIAL_NS::internal::Sturm_root_rep<S,I> >:
    public numeric_limits<typename CGAL_POLYNOMIAL_NS::internal::Sturm_root_rep<S,I>::NT >
  {
  public:
    typedef numeric_limits<typename CGAL_POLYNOMIAL_NS::internal::Sturm_root_rep<S,I>::NT > P;
    typedef CGAL_POLYNOMIAL_NS::internal::Sturm_root_rep<S,I> T;
    static const bool is_specialized = true;
    static T min BOOST_PREVENT_MACRO_SUBSTITUTION () throw() {return T((P::min)());}
    static T max BOOST_PREVENT_MACRO_SUBSTITUTION () throw() {return T((P::max)());}
    /*static const int digits =0;
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
      static const int max_exponent10=0;*/
    static const bool has_infinity=true;
    /*static const bool has_quiet_NaN = false;
      static const bool has_signaling_NaN= false;
      static const float_denorm_style has_denorm= denorm_absent;
      static const bool has_denorm_loss = false;
    */
    static T infinity() throw() {return T::infinity();}
    /*static T quiet_NaN() throw(){return T(0);}
      static T denorm_min() throw() {return T(0);}
      static const bool is_iec559=false;
      static const bool is_bounded =false;
      static const bool is_modulo= false;
      static const bool traps = false;
      static const bool tinyness_before =false;
      static const float_round_style round_stype = round_toward_zero;*/
  };
};
#endif                                            // CGAL_POLYNOMIAL_STURM_ROOT_REP_H
