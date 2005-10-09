#ifndef CGAL_POLYNOMIAL_STURM_ROOT_REP_H
#define CGAL_POLYNOMIAL_STURM_ROOT_REP_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Rational/Sign_Sturm_sequence.h>

#include <cfloat>

//#include <CGAL/Polynomial/internal/Bisection.h>


CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

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
      return Sign(s1 * s2);
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

  bool is_rational() const {
    return is_exact();
  }

  NT to_rational() const {
    Polynomial_precondition( is_rational() );
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
    } else {
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
	} else {
	  multiplicity_ += 2;
	  q = differentiate( differentiate(q) );
	}
      } else {
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

#if 0
  void set_lower(const Self& r, typename Interval::Endpoint b,
		 const Sign& s_l ) const
  {
    typename Interval::NT l;
    if ( b == Interval::LOWER ) {
      l = r.ivl.lower_bound();
    } else {
      l = r.ivl.upper_bound();
    }

    ivl.set_lower(l);
    s_lower = s_l;
  }

  void set_upper(const Self& r, typename Interval::Endpoint b,
		 const Sign& s_u) const
  {
    typename Interval::NT u;
    if ( b == Interval::LOWER ) {
      u = r.ivl.lower_bound();
    } else {
      u = r.ivl.upper_bound();
    }

    ivl.set_upper(u);
    s_upper = s_u;
  }

  void set_interval(const Self& r, typename Interval::Endpoint b) const
  {
    typename Interval::NT x;
    if ( b == Interval::LOWER ) {
      x = r.ivl.lower_bound();
    } else {
      x = r.ivl.upper_bound();
    }
    ivl = Interval(x);
    s_lower = s_upper = CGAL::ZERO;
  }
#endif

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
      } else {
	if ( upper_bound() <= r.lower_bound() ) {
	  return CGAL::SMALLER;
	} else if ( lower_bound() >= r.upper_bound() ) {
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
    } else if ( r.is_exact() ) {
	
      // we have already checked the case that this root is also
      // known exactly
      //	std::cout << "second is exact" << std::endl;
      if ( upper_bound() <= r.lower_bound() ) {
	//	  std::cout << "this.upper bound <= other.lower bound" << std::endl;
	return CGAL::SMALLER;
      } else if ( lower_bound() >= r.upper_bound() ) {
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
    } else {
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
      } else if ( lower_bound() >= r.upper_bound() ) {
	//	  std::cout << "this.lower bound >= other.upper bound" << std::endl;
	return CGAL::LARGER;
      } else if ( upper_bound() > r.lower_bound() &&
		  lower_bound() < r.upper_bound() ) {
	//	  std::cout << "intersecting intervals" << std::endl;
	return compare_intersecting(r);
      }
    }


    bool this_line_should_not_have_been_reached = false;
    CGAL_assertion( this_line_should_not_have_been_reached );
    if (0) this_line_should_not_have_been_reached= false;
    return CGAL::EQUAL;
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
	} else {
	  if ( s == sign_lower() ) {
	    return CGAL::LARGER;
	  } else {
	    //	      set_upper( r, Interval::UPPER, s );
	    set_upper( r.upper_bound(), s );
	    return compare_same_interval(r);
	  }
	}
      } else if ( r.upper_bound() > upper_bound() ) {
	const This* new_this = static_cast<const This*>(this);
	return opposite(r.compare_intersecting( *new_this ));
      } else if ( upper_bound() == r.upper_bound() ) {
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
	} else {
	  if ( s == sign_upper() ) {
	    return CGAL::SMALLER;
	  } else {
	    //	      set_lower( r, Interval::LOWER, s );
	    set_lower( r.lower_bound(), s );
	    return compare_same_interval(r);
	  }
	}
      } else if ( r.lower_bound() < lower_bound() ) {
	const This* new_this = static_cast<const This*>(this);
	return opposite(r.compare_intersecting( *new_this ));
      } else if ( lower_bound() == r.lower_bound() ) {
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
	   s_ru_at_this != CGAL::ZERO	) {
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
  } // end-compare_intersecting


    //--------------------------------------------------------------
  void subdivide() const
  {
#if 1
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
    } else {
      ivl = second_half;
    }
#else
    typename Interval::NT mid = ivl.middle();

    Sign_at_functor sign_at_f(this);

    Sign s_mid = apply(sign_at_f, mid);

    if ( s_mid == CGAL::ZERO ) {
      ivl.set_lower(mid);
      ivl.set_upper(mid);
      return;
    }

    if ( s_mid == sign_upper() ) {
      ivl.set_upper(mid);
    } else {
      ivl.set_lower(mid);
    }
#endif
  }
#if 1
  template<class Child>
  Polynomial compute_simple(Ring_tag, const Child&) const
  {
    Polynomial gcd = sseq.exact( sseq.exact_size() - 1 );
    return p_.pseudo_quotient( gcd );
  }

  template<class Child>
  Polynomial compute_simple(Field_tag, const Child&) const
  {
    Polynomial gcd = sseq.exact( sseq.exact_size() - 1 );
    return p_ / gcd;
  }
#endif
  Polynomial compute_simple(Ring_tag, const Self&) const
  {
    Polynomial gcd = sseq[sseq.size() - 1];
    return tr_.pseudo_quotient_object()(p_, gcd);
  }

    
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
    } else {
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
  static Self infinity(const Traits& tr = Traits())
  {
    return Self(-1, 0, tr);
  }


  //================
  // CONSTRUCTORS
  //================
protected:
  Sturm_root_rep(int i, int, const Traits& tr)
    : idx(i), ivl(), p_(), sseq(), s_lower(CGAL::ZERO),
      s_upper(CGAL::ZERO), multiplicity_(0), tr_(tr) {}

public:
  Sturm_root_rep(const Traits& tr = Traits())
    : idx(-1), ivl(), p_(), sseq(), s_lower(CGAL::ZERO),
      s_upper(CGAL::ZERO), multiplicity_(0), tr_(tr) {}
  //-------
  template<class T>
  Sturm_root_rep(const T& a, const Traits& tr = Traits())
    : idx(1), ivl(a), p_(), sseq(), s_lower(CGAL::ZERO),
      s_upper(CGAL::ZERO), multiplicity_(0), tr_(tr) {}
  //-------
  Sturm_root_rep(const Interval& ivl,
		 const Polynomial& p,
		 const Standard_sequence& sseq, int idx,
		 const Traits& tr)
    : idx(idx), ivl(ivl), p_(p), sseq(sseq),
      s_lower(CGAL::ZERO), s_upper(CGAL::ZERO), multiplicity_(0),
      tr_(tr)
  {
    Sign_at_functor sign_at_p(this);
    s_lower = apply(sign_at_p, ivl.lower_bound());
    s_upper = apply(sign_at_p, ivl.upper_bound());

    //      s_lower = sign_at( *this, Interval::LOWER );
    //      s_upper = sign_at( *this, Interval::UPPER );
  }
  //-------
  Sturm_root_rep(const Self& other)
    : idx(other.idx), ivl(), p_(other.p_), sseq(other.sseq),
      s_lower(other.s_lower),	s_upper(other.s_upper),
      multiplicity_(other.multiplicity_), tr_(other.tr_)
  {
    if ( other.idx >= 1 ) {
      ivl = other.ivl;
    }
  }

  //=======================
  // ASSIGNMENT OPERATOR
  //=======================
  const Self& operator=(const Self& other)
  {
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
    } else if ( idx == -2 ) {
      return Self(-1,0, tr_);
    } else if (idx == 1) {
      // HACK by Daniel
      return Self(-ivl, tr_);
    } else {

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
    }

    return multiplicity_;
  }

  //====================
  // THE PARITY METHOD
  //====================
  bool is_even_multiplicity() const
  {
    if ( idx < 0 ) { return false; }

    if ( is_exact() ) {
      return (multiplicity() % 2 == 0);
    } else {
      Sign_at_functor sign_at_f(this);

      Sign sl = apply(sign_at_f, ivl.lower_bound());
      Sign su = apply(sign_at_f, ivl.upper_bound());
      return (sl == su);
    }
  }


  //==========================
  // ACCESS TO THE INTERVAL
  //==========================
  typename Interval::NT lower_bound() const { return ivl.lower_bound(); }
  typename Interval::NT upper_bound() const { return ivl.upper_bound(); }
  Interval    interval()    const { return ivl; }
  int         index()       const { return idx; }
  const Polynomial&     polynomial() const { return p_; }
  //    Polynomial  simple_polynomial() const { return p; }

public:
  //======================
  // CONVERTOR TO double
  //======================
  double to_double(double acc = 1e-10) const
  {
    if ( idx == -1 ) {
      return DBL_MAX * DBL_MAX;
    }

    if ( idx == -2 ) {
      return -DBL_MAX * DBL_MAX;
    }

    if ( is_exact() ) {
      return CGAL::to_double(lower_bound());
    }

    Exact_nt xacc(acc);

#if 0
    if ( ivl.is_double() ) {
      while ( ivl.middle().is_double() ) {
	subdivide();
	if ( is_exact() ) { break; }
      }
    } else {
      while ( ivl.width() > xacc ) {
	subdivide();
	if ( is_exact() ) { break; }
      }
    }
#else
    while ( ivl.approximate_width() > acc ) {
      subdivide();
      if ( is_exact() ) { break; }
    }
#endif
      
    return CGAL::to_double( lower_bound() );
  }

  //========================
  // CONVERTOR TO interval
  //========================
  std::pair<double,double> to_interval() const
  {
    to_double();

    std::pair<double,double> i_low =
      CGAL::to_interval(ivl.lower_bound());
    std::pair<double,double> i_high =
      CGAL::to_interval(ivl.upper_bound());

    return std::pair<double,double>(i_low.first, i_high.second);
  }

  //================
  // STREAM WRITER
  //================
  template<class Stream>
  Stream& write(Stream& os) const
  {
    if ( idx == -2 ) {
      os << "-inf";
    } else if ( idx == -1 ) {
      os << "inf";
    } else {
      os << "{" << ivl << ", " << idx << "}";
    }
    Self copy = *this;
    os << " = " << copy.to_double();
    return os;
  }

protected:
  int                    idx;
  mutable Interval       ivl;
  Polynomial             p_;
  Standard_sequence      sseq;
  mutable Sign           s_lower, s_upper;
  mutable unsigned int   multiplicity_;
  Traits                 tr_;
};


template<class Stream, class S, class I>
Stream&
operator<<(Stream& os, const Sturm_root_rep<S,I>& r)
{
  return r.write(os);
}



CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE



CGAL_BEGIN_NAMESPACE

template<class S, class I>
double
to_double(const POLYNOMIAL_NS::internal::Sturm_root_rep<S,I>& r)
{
  return r.to_double();
}

template<class S, class I>
std::pair<double,double>
to_interval(const POLYNOMIAL_NS::internal::Sturm_root_rep<S,I>& r)
{
  return r.to_interval();
}




CGAL_END_NAMESPACE

namespace std {
  template <class S, class I>
  struct numeric_limits<POLYNOMIAL_NS::internal::Sturm_root_rep<S,I> >: 
    public numeric_limits<typename POLYNOMIAL_NS::internal::Sturm_root_rep<S,I>::NT > {
    typedef numeric_limits<POLYNOMIAL_NS::internal::Sturm_root_rep<S,I> > P;
    typedef POLYNOMIAL_NS::internal::Sturm_root_rep<S,I> T;
    static const bool is_specialized = true;
    static T min() throw() {return T(P::min());}
    static T max() throw() {return T(P::max());}
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


#endif // CGAL_POLYNOMIAL_STURM_ROOT_REP_H
