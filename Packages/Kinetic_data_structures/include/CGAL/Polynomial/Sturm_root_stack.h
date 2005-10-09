#ifndef CGAL_STURM_LAZY_SOLVER_H
#define CGAL_STURM_LAZY_SOLVER_H

#include <CGAL/Polynomial/basic.h>
//#include <CGAL/Polynomial/utilities.h>
//#include <CGAL/Polynomial/make_simple.h>

// the Standard sequence
#include <CGAL/Polynomial/internal/Rational/Standard_sequence.h>

// representations for roots and intervals
#include <CGAL/Polynomial/internal/Sturm_root_rep.h>
#include <CGAL/Polynomial/internal/Sturm_isolating_interval.h>

// the root counters
//#include <CGAL/Polynomial/internal/Rational/Sturm_root_counter.h>
//#include <CGAL/Polynomial/internal/Rational/Descartes_root_counter.h>
//#include <CGAL/Polynomial/internal/Rational/Bezier_root_counter.h>

#include <list>

#include <CGAL/Polynomial/Kernel.h>

CGAL_POLYNOMIAL_BEGIN_NAMESPACE


//================================================================
//================================================================
// ************** the lazy Sturm solver **************************
//================================================================
//================================================================

#if 0
namespace Polynomiali {

  template<class Tag>
  struct Which_root_counter_chooser
  {
    template<class P, class M>
    struct Which {
      typedef
      CGAL::Polynomial::internal::Sturm_root_counter<P,M>
      Root_count;
    };
  };


  template<>
  struct Which_root_counter_chooser<Descartes_tag>
  {
    template<class P, class M>
    struct Which {
      typedef
      CGAL::Polynomial::internal::Descartes_root_counter<P>
      Root_count;
    };
  };



  template<>
  struct Which_root_counter_chooser<Bezier_tag>
  {
    template<class P, class M>
    struct Which {
      typedef
      CGAL::Polynomial::internal::Bezier_root_counter<P>
      Root_count;
    };
  };


  template<class STag, class P, class M>
  struct Which_root_counter
  {
    typedef Which_root_counter_chooser<STag> Which;
    typedef typename Which::template Which<P,M>::Root_count
    Root_count;
  };


} // namespace Polynomiali
#endif


struct Sturm_tag {};

template<class T, class S = Sturm_tag, class M = Field_tag,
	 bool Drop_even_roots = false>
class Sturm_root_stack
{
public:
  typedef T                             Traits;
  typedef typename Traits::Function     Polynomial;
  typedef typename Traits::Isolating_interval     Interval;
  typedef typename Traits::NT           NT;
  typedef typename Traits::Sign_at      Sign_at;

  typedef M                             Method_tag;
  typedef S                             Subdivision_tag;

  enum { DROP_EVEN_ROOTS = Drop_even_roots };

  typedef typename Traits::Standard_sequence Standard_sequence;

protected:
  //  typedef CGAL::Sign                             Sign;
  //  typedef internal::Sign_at<Polynomial>          Sign_at;
  typedef std::pair<unsigned int, unsigned int>  uint_pair;

  //  typedef internal::Sturm_isolating_interval<NT>           Interval;
  typedef std::list<Interval>                    Interval_container;
  typedef std::pair<Interval,Interval>           Interval_pair;
  typedef typename Traits::Root_count            Root_count;
  //  typedef typename
  //  Polynomiali::Which_root_counter<S,P,M>::Root_count
  //  Root_count;

  typedef Sturm_root_stack<T,S,M,DROP_EVEN_ROOTS>  Self;

public:
  typedef internal::Sturm_root_rep<Self,Interval>  Root;

protected:
  //====================
  // PROTECTED METHODS
  //====================


#if 0
  // ACCESS METHODS TO THE INTERVAL LIST AND SIGN VARIATIONS LIST
  void pop_front( Descartes_tag ) const
  {
    i_list.pop_front();
    nroots.pop_front();
  }
#endif
  void pop_front( Sturm_tag ) const
  {
    i_list.pop_front();
    s_variations.pop_front();
  }

  unsigned int num_roots_in_first_interval( Sturm_tag ) const
  {
    uint_pair ip = s_variations.front();
    return ip.first - ip.second;
  }
#if 0
  unsigned int num_roots_in_first_interval( Descartes_tag ) const
  {
    return nroots.front();
  }
#endif

  void push_front(const Interval& ivl, unsigned int sva,
		  unsigned int svb, Sturm_tag) const
  {
    i_list.push_front(ivl);
    s_variations.push_front( uint_pair(sva, svb) );
  }
#if 0
  void push_front(const Interval& ivl,
		  unsigned int nr, Descartes_tag) const
  {
    i_list.push_front(ivl);
    nroots.push_front(nr);
  }
#endif
  void push_front(const Interval& ivl, Sturm_tag) const
  {
    i_list.push_front(ivl);
    s_variations.push_front( uint_pair(1, 0) );
  }
#if 0
  void push_front(const Interval& ivl, Descartes_tag) const
  {
    i_list.push_front(ivl);
    nroots.push_front(1);
  }
#endif

  //------------------------------------------------------------------

  //------------------------------------------------------------------
  // the method for subdivision; it guarantees that the first interval
  // contains only one root
  template<class SA, class RC>
  void subdivide(Sturm_tag tag, const SA& sign_at_p,
		 const RC& rc) const
  {
    if ( i_list.size() == 0 ) { return; }

    if ( num_roots_in_first_interval(tag) == 1 ) { return; }

    while ( num_roots_in_first_interval(tag) > 1 ) {
      Interval ivl = i_list.front();

#if 1
      Interval_pair ivl_pair = ivl.split();
      Interval first_half = ivl_pair.first;
      Interval second_half = ivl_pair.second;

      unsigned int sv_a = s_variations.front().first;
      unsigned int sv_b = s_variations.front().second;

      Sign s_mid = 
      	first_half.apply_to_endpoint(sign_at_p, Interval::UPPER);

      if ( s_mid == CGAL::ZERO ) {
	Interval mid_ivl = ivl;

	Interval mid = first_half.endpoint_interval(Interval::UPPER);

	Sign s_midp, s_midm;
	unsigned int sv_midm(0), sv_midp(0);

	do {
	  mid_ivl = mid_ivl.middle_half();
	  s_midm = mid_ivl.apply_to_endpoint(sign_at_p, Interval::LOWER);
	  s_midp = mid_ivl.apply_to_endpoint(sign_at_p, Interval::UPPER);

	  if ( s_midm != CGAL::ZERO && s_midp != CGAL::ZERO ) {
	    sv_midm = mid_ivl.apply_to_endpoint(rc, Interval::LOWER);
	    sv_midp = mid_ivl.apply_to_endpoint(rc, Interval::UPPER);
	  }

	} while ( s_midm == CGAL::ZERO || s_midp == CGAL::ZERO ||
		  sv_midm != sv_midp + 1 );

	pop_front( tag );
	if ( sv_midp > sv_b ) {
	  Interval ivl_ = ivl.split_at( mid_ivl.upper_bound() ).second;

	  push_front(ivl_, sv_midp, sv_b, tag);
	}

	push_front(mid, tag);

	if ( sv_a > sv_midm ) {
	  Interval ivl_ = ivl.split_at( mid_ivl.lower_bound() ).first;

	  push_front(ivl_, sv_a, sv_midm, tag);
	}
      } else {
	unsigned int sv_mid =
	  first_half.apply_to_endpoint(rc, Interval::UPPER);

	unsigned int n1 = sv_a - sv_mid; // num roots in (a, mid)
	unsigned int n2 = sv_mid - sv_b; // num roots in (mid, b)

	CGAL_assertion( n1 + n2 == num_roots_in_first_interval(tag) );

	pop_front( tag );
	if ( n2 > 0 ) {
	  push_front(second_half, sv_mid, sv_b, tag);
	}
	if ( n1 > 0 ) {
	  push_front(first_half, sv_a, sv_mid, tag);
	}
      } // end-if
#else
      //      Interval_pair ivl_pair = ivl.split();
      //      Interval first_half = ivl_pair.first;
      //      Interval second_half = ivl_pair.second;

      typename Interval::NT a = ivl.lower_bound();
      typename Interval::NT b = ivl.upper_bound();

      typename Interval::NT mid = ivl.middle();

      typename Interval::NT a1 = a, b1 = b;

      unsigned int sv_a = s_variations.front().first;
      unsigned int sv_b = s_variations.front().second;

      Sign s_mid = apply(sign_at_p, mid);

      if ( s_mid == CGAL::ZERO ) {
	Sign s_midp, s_midm;
	unsigned int sv_midm(0), sv_midp(0);

	do {
	  a1 = midpoint(a1, mid);
	  b1 = midpoint(b1, mid);
	  s_midm = apply(sign_at_p, a1);
	  s_midp = apply(sign_at_p, b1);

	  if ( s_midm != CGAL::ZERO && s_midp != CGAL::ZERO ) {
	    sv_midm = apply(rc, a1);
	    sv_midp = apply(rc, b1);
	  }

	} while ( s_midm == CGAL::ZERO || s_midp == CGAL::ZERO ||
		  sv_midm != sv_midp + 1 );

	pop_front( tag );
	if ( sv_midp > sv_b ) {
	  Interval ivl_(b1, b);
	  push_front(ivl_, sv_midp, sv_b, tag);
	}

	Interval imid(mid);
	push_front(imid, tag);

	if ( sv_a > sv_midm ) {
	  Interval ivl_(a, a1);
	  push_front(ivl_, sv_a, sv_midm, tag);
	}
      } else {
	unsigned int sv_mid = apply(rc, mid);

	unsigned int n1 = sv_a - sv_mid; // num roots in (a, mid)
	unsigned int n2 = sv_mid - sv_b; // num roots in (mid, b)

	CGAL_assertion( n1 + n2 == num_roots_in_first_interval(tag) );

	pop_front( tag );
	if ( n2 > 0 ) {
	  Interval ivl_(mid, b);
	  push_front(ivl_, sv_mid, sv_b, tag);
	}
	if ( n1 > 0 ) {
	  Interval ivl_(a, mid);
	  push_front(ivl_, sv_a, sv_mid, tag);
	}
      } // end-if

#endif

    } // end-while
  } // end subdivide method
#if 0
  //------------------------------------------------------------------
  // the method for subdivision; it guarantees that the first interval
  // contains only one root
  template<class SA, class RC>
  void subdivide(Descartes_tag tag, const SA& sign_at_p,
		 const RC& rc) const
  {
    if ( i_list.size() == 0 ) { return; }

    if ( num_roots_in_first_interval(tag) == 1 ) { return; }

    // MK: WRITE THIS METHOD AS IN THE CASE OF Sturm_tag

    while ( num_roots_in_first_interval(tag) > 1 ) {
      //      std::cout << "." << std::flush;
      Interval ivl = i_list.front();
      Interval_pair ivl_pair = ivl.split();
      Interval first_half = ivl_pair.first;
      Interval second_half = ivl_pair.second;

      //      Sign s_mid = sign_at_p(mid);
      Sign s_mid =
	ivl_pair.first.apply_to_bound(sign_at_p, Interval::UPPER);

      if ( s_mid == CGAL::ZERO ) {
	Interval mid_ivl = ivl;

	Sign s_midp, s_midm;
	unsigned int nr(0);

	do {
	  mid_ivl = mid_ivl.middle_half();
	  s_midm = mid_ivl.apply_to_bound(sign_at_p, Interval::LOWER);
	  s_midp = mid_ivl.apply_to_bound(sign_at_p, Interval::UPPER);

	  //	  s_midm = sign_at_p(midm);
	  //	  s_midp = sign_at_p(midp);

	  if ( s_midm != CGAL::ZERO && s_midp != CGAL::ZERO ) {
	    nr = mid_ivl.apply_to_interval(rc);
	  }
	} while ( s_midm == CGAL::ZERO || s_midp == CGAL::ZERO ||
		  nr != 1 );


	//	Interval_NT  a = ivl.lower_bound();
	//	Interval_NT  b = ivl.upper_bound();

	Interval ivl_upper =
	  ivl.split_at(mid_ivl, Interval::UPPER).second;
	  //	  Interval::create_from_endpoints(mid_ivl, Interval::UPPER,
	  //					  ivl, Interval::UPPER);

	Interval ivl_lower =
	  ivl.split_at(mid_ivl, Interval::LOWER).first;
	  //	  Interval::create_from_endpoints(ivl, Interval::LOWER,
	  //					  mid_ivl, Interval::LOWER);

	unsigned int np = ivl_upper.apply_to_interval(rc);
	unsigned int nm = ivl_lower.apply_to_interval(rc);

	pop_front( tag );

	if ( np > 0 ) {
	  push_front(ivl_upper, np, tag);
	}

	Interval mid = first_half.collapse_to_bound(Interval::UPPER);
	push_front(mid, tag);

	if ( nm > 0 ) {
	  push_front(ivl_lower, nm, tag);
	}

      } else {
	//	Interval_NT  a = ivl.lower_bound();
	//	Interval_NT  b = ivl.upper_bound();

	int n1 = first_half.apply_to_interval(rc);
	int n2 = second_half.apply_to_interval(rc);

	pop_front( tag );
	if ( n2 > 0 ) {
	  push_front(second_half, n2, tag);
	}
	if ( n1 > 0 ) {
	  push_front(first_half, n1, tag);
	}

	if ( n1 + n2 == 0 ) { break; }
      } // end-if
    } // end-while
  } // end subdivide method
#endif

  //------------------------------------------------------------------

  // the method for initializing the interval list
  template<class RC, class SA>
  void initialize_interval_list(const Interval initial_ivl,
				Sturm_tag tag, const SA& sign_at_p,
				const RC& rc)
  {
#if 1
    Sign s_a = initial_ivl.apply_to_endpoint(sign_at_p, Interval::LOWER);
    Sign s_b = initial_ivl.apply_to_endpoint(sign_at_p, Interval::UPPER);

    if ( s_a != CGAL::ZERO && s_b != CGAL::ZERO ) {

      unsigned int sv_a = initial_ivl.apply_to_endpoint(rc, Interval::LOWER);
      unsigned int sv_b = initial_ivl.apply_to_endpoint(rc, Interval::UPPER);

      if ( sv_a > sv_b ) {
	push_front(initial_ivl, sv_a, sv_b, tag);
      }
      return;
    }

    if ( s_a == CGAL::ZERO && s_b != CGAL::ZERO ) {
      Interval ia = initial_ivl.interval_around_endpoint(Interval::LOWER);
      Interval iam(ia.lower_bound(), initial_ivl.lower_bound());
      Interval iap = initial_ivl.split_at(ia.upper_bound()).first;

      Sign s_ap, s_am;
      unsigned int sv_ap(0), sv_am(0);

      do {
	iap = iap.split().first;
	iam = iam.split().second;

	s_ap = iap.apply_to_endpoint(sign_at_p, Interval::UPPER);
	s_am = iam.apply_to_endpoint(sign_at_p, Interval::LOWER);

	if ( s_ap != CGAL::ZERO ) {
	  sv_ap = iap.apply_to_endpoint(rc, Interval::UPPER);
	}
	if ( s_am != CGAL::ZERO ) {
	  sv_am = iam.apply_to_endpoint(rc, Interval::LOWER);
	}

      } while ( s_ap == CGAL::ZERO || s_am == CGAL::ZERO ||
		sv_ap + 1 != sv_am );


      unsigned int sv_b =
	initial_ivl.apply_to_endpoint(rc, Interval::UPPER);

      Interval ivl_upper =
	initial_ivl.split_at(iap.upper_bound()).second;

      push_front(ivl_upper, sv_ap, sv_b, tag);

      Interval ivl_a = initial_ivl.endpoint_interval(Interval::LOWER);
      push_front(ivl_a, tag);
      return;
    }

    if ( s_a != CGAL::ZERO && s_b == CGAL::ZERO ) {
      Interval ib = initial_ivl.interval_around_endpoint(Interval::UPPER);

      Interval ibm =
	initial_ivl.split_at(ib.lower_bound()).second;

      Interval ibp(initial_ivl.upper_bound(), ib.upper_bound());

      Sign s_bp, s_bm;
      unsigned int sv_bp(0), sv_bm(0);
      do {
	ibp = ibp.split().first;
	ibm = ibm.split().second;
	s_bp = ibp.apply_to_endpoint(sign_at_p, Interval::UPPER);
	s_bm = ibm.apply_to_endpoint(sign_at_p, Interval::LOWER);

	if ( s_bp != CGAL::ZERO ) {
	  sv_bp = ibp.apply_to_endpoint(rc, Interval::UPPER);
	}
	if ( s_bm != CGAL::ZERO ) {
	  sv_bm = ibm.apply_to_endpoint(rc, Interval::LOWER);
	}
      } while ( s_bp == CGAL::ZERO || s_bm == CGAL::ZERO ||
		sv_bp != sv_bm - 1 );

      unsigned int sv_a =
	initial_ivl.apply_to_endpoint(rc, Interval::LOWER);

      Interval ivl_b = initial_ivl.endpoint_interval(Interval::UPPER);
      push_front(ivl_b, tag);

      Interval ivl_lower =
	initial_ivl.split_at( ibm.lower_bound() ).first;
					
      push_front(ivl_lower, sv_a, sv_bm, tag);
      return;
    }

    if ( s_a == CGAL::ZERO && s_b == CGAL::ZERO ) {
      // both a and b are roots
      Interval ia =
	initial_ivl.interval_around_endpoint(Interval::LOWER);
      Interval ib =
	initial_ivl.interval_around_endpoint(Interval::UPPER);

      Interval iam(ia.lower_bound(), initial_ivl.lower_bound());
      Interval iap = initial_ivl.split_at(ia.upper_bound()).first;
      Interval ibm = initial_ivl.split_at(ib.lower_bound()).second;
      Interval ibp(initial_ivl.upper_bound(), ib.upper_bound());

      Sign s_ap, s_am, s_bp, s_bm;
      unsigned int sv_am(0), sv_ap(0), sv_bm(0), sv_bp(0);
      unsigned int nbig(0), nsmall(0);

      do {
	iap = iap.split().first;
	iam = iam.split().second;
	s_ap = iap.apply_to_endpoint(sign_at_p, Interval::UPPER);
	s_am = iam.apply_to_endpoint(sign_at_p, Interval::LOWER);

	ibp = ibp.split().first;
	ibm = ibm.split().second;
	s_bp = ibp.apply_to_endpoint(sign_at_p, Interval::UPPER);
	s_bm = ibm.apply_to_endpoint(sign_at_p, Interval::LOWER);

	if ( s_ap != CGAL::ZERO && s_bm != CGAL::ZERO ) {
	  sv_ap = iap.apply_to_endpoint(rc, Interval::UPPER);
	  sv_bm = ibm.apply_to_endpoint(rc, Interval::LOWER);

	  nsmall = sv_ap - sv_bm; // num roots in (ap, bm)
	}
	if ( s_am != CGAL::ZERO && s_bp != CGAL::ZERO) {
	  sv_am = iam.apply_to_endpoint(rc, Interval::LOWER);
	  sv_bp = ibp.apply_to_endpoint(rc, Interval::UPPER);

	  nbig = sv_am - sv_bp;   // num roots in (am, bp)
	}
      } while ( s_ap == CGAL::ZERO || s_am == CGAL::ZERO ||
		s_bp == CGAL::ZERO || s_bm == CGAL::ZERO ||
		nbig != nsmall + 2 );

      Interval ivl_b = initial_ivl.endpoint_interval(Interval::UPPER);
      push_front(ivl_b, tag);

      Interval ivl_interior(iap.upper_bound(), ibm.lower_bound());

      push_front(ivl_interior, sv_ap, sv_bm, tag);

      Interval ivl_a = initial_ivl.endpoint_interval(Interval::LOWER);
      push_front(ivl_a, tag);
      return;
    }
#else
    typename Interval::NT a = initial_ivl.lower_bound();
    typename Interval::NT b = initial_ivl.upper_bound();

    Sign s_a = apply(sign_at_p, a);
    Sign s_b = apply(sign_at_p, b);

    if ( s_a != CGAL::ZERO && s_b != CGAL::ZERO ) {

      unsigned int sv_a = apply(rc, a);
      unsigned int sv_b = apply(rc, b);

      if ( sv_a > sv_b ) {
	push_front(initial_ivl, sv_a, sv_b, tag);
      }
      return;
    }

    if ( s_a == CGAL::ZERO && s_b != CGAL::ZERO ) {
      Interval ia =
	initial_ivl.interval_around_endpoint(Interval::LOWER);
      typename Interval::NT am = ia.lower_bound();
      typename Interval::NT ap = ia.upper_bound();

      Sign s_ap, s_am;
      unsigned int sv_ap(0), sv_am(0);

      do {
	ap = midpoint(ap, a);
	am = midpoint(am, a);

	s_ap = apply(sign_at_p, ap);
	s_am = apply(sign_at_p, am);

	if ( s_ap != CGAL::ZERO ) {
	  sv_ap = apply(rc, ap);
	}
	if ( s_am != CGAL::ZERO ) {
	  sv_am = apply(rc, am);
	}

      } while ( s_ap == CGAL::ZERO || s_am == CGAL::ZERO ||
		sv_ap + 1 != sv_am );

      unsigned int sv_b = apply(rc, b);

      Interval ivl_upper(ap, b);
      push_front(ivl_upper, sv_ap, sv_b, tag);

      Interval ivl_a(a);
      push_front(ivl_a, tag);
      return;
    }

    if ( s_a != CGAL::ZERO && s_b == CGAL::ZERO ) {
      Interval ib =
	initial_ivl.interval_around_endpoint(Interval::UPPER);
      typename Interval::NT bm = ib.lower_bound();
      typename Interval::NT bp = ib.upper_bound();

      Sign s_bp, s_bm;
      unsigned int sv_bp(0), sv_bm(0);
      do {
	bp = midpoint(bp, b);
	bm = midpoint(bm, b);
	s_bp = apply(sign_at_p, bp);
	s_bm = apply(sign_at_p, bm);

	if ( s_bp != CGAL::ZERO ) {
	  sv_bp = apply(rc, bp);
	}
	if ( s_bm != CGAL::ZERO ) {
	  sv_bm = apply(rc, bm);
	}
      } while ( s_bp == CGAL::ZERO || s_bm == CGAL::ZERO ||
		sv_bp != sv_bm - 1 );

      unsigned int sv_a = apply(rc, a);

      Interval ivl_b(b);
      push_front(ivl_b, tag);

      Interval ivl_lower(a, bm);
      push_front(ivl_lower, sv_a, sv_bm, tag);
      return;
    }

    if ( s_a == CGAL::ZERO && s_b == CGAL::ZERO ) {
      // both a and b are roots
      Interval ia =
	initial_ivl.interval_around_endpoint(Interval::LOWER);
      Interval ib =
	initial_ivl.interval_around_endpoint(Interval::UPPER);

      typename Interval::NT am = ia.lower_bound();
      typename Interval::NT ap = ia.upper_bound();
      typename Interval::NT bm = ib.lower_bound();
      typename Interval::NT bp = ib.upper_bound();
      
      Sign s_ap, s_am, s_bp, s_bm;
      unsigned int sv_am(0), sv_ap(0), sv_bm(0), sv_bp(0);
      unsigned int nbig(0), nsmall(0);

      do {
	ap = midpoint(ap, a);
	am = midpoint(am, a);
	bp = midpoint(bp, b);
	bm = midpoint(bm, b);

	s_ap = apply(sign_at_p, ap);
	s_am = apply(sign_at_p, am);
	s_bp = apply(sign_at_p, bp);
	s_bm = apply(sign_at_p, bm);

	if ( s_ap != CGAL::ZERO && s_bm != CGAL::ZERO ) {
	  sv_ap = apply(rc, ap);
	  sv_bm = apply(rc, bm);

	  nsmall = sv_ap - sv_bm; // num roots in (ap, bm)
	}
	if ( s_am != CGAL::ZERO && s_bp != CGAL::ZERO) {
	  sv_am = apply(rc, am);
	  sv_bp = apply(rc, bp);

	  nbig = sv_am - sv_bp;   // num roots in (am, bp)
	}
      } while ( s_ap == CGAL::ZERO || s_am == CGAL::ZERO ||
		s_bp == CGAL::ZERO || s_bm == CGAL::ZERO ||
		nbig != nsmall + 2 );

      Interval ivl_b(b);
      push_front(ivl_b, tag);

      Interval ivl_interior(ap, bm);
      push_front(ivl_interior, sv_ap, sv_bm, tag);

      Interval ivl_a(a);
      push_front(ivl_a, tag);
      return;
    }
    
#endif
    //bool this_line_should_not_have_been_reached = false;
    CGAL_postcondition_msg(0,"this_line_should_not_have_been_reached" );
    //if (0) this_line_should_not_have_been_reached=false;
  }

#if 0
  //------------------------------------------------------------------
  template<class RC, class SA>
  void initialize_interval_list(const Interval& initial_ivl,
				Bezier_tag, const SA& sign_at_p,
				const RC& rc)
  {
    initialize_interval_list(initial_ivl, Descartes_tag(),
			     sign_at_p, rc);
  }

  //------------------------------------------------------------------

  // the method for initializing the interval list
  template<class RC, class SA>
  void initialize_interval_list(const Interval& initial_ivl,
				Descartes_tag tag,
				const SA& sign_at_p,
				const RC& rc)
  {
    // MK: WRITE THIS METHOD AS IN THE CASE OF Sturm_tag

    Sign s_a = initial_ivl.apply_to_bound(sign_at_p, Interval::LOWER);
    Sign s_b = initial_ivl.apply_to_bound(sign_at_p, Interval::UPPER);

    if ( s_a != CGAL::ZERO && s_b != CGAL::ZERO ) {
      int nr = initial_ivl.apply_to_interval(rc);
      
      if ( nr > 0 ) {
	push_front(initial_ivl, nr, tag);
      }
      return;
    }

    if ( s_a == CGAL::ZERO && s_b != CGAL::ZERO ) {
      Interval ia = initial_ivl.interval_around_bound(Interval::LOWER);
      Interval iam(ia.lower_bound(), initial_ivl.lower_bound());
      Interval iap = initial_ivl.split_at(ia, Interval::UPPER).first;
      //      Interval_NT ap = a_plus(a, b);
      //      Interval_NT am = a - NT(1);

      Sign s_ap, s_am;
      int np(0), nm(0);

      do {
	iap = iap.split().first;
	iam = iam.split().second;

	s_ap = iap.apply_to_bound(sign_at_p, Interval::UPPER);
	s_am = iam.apply_to_bound(sign_at_p, Interval::LOWER);

	//	ap = Interval::midpoint(a, ap);
	//	am = Interval::midpoint(a, am);
	//	s_ap = sign_at_p(ap);
	//	s_am = sign_at_p(am);

	if ( s_ap != CGAL::ZERO ) {
	  np = iap.apply_to_interval(rc);
	  //	  np = root_counter_(ap, b);
	}
	if ( s_am != CGAL::ZERO ) {
	  nm = iam.apply_to_interval(rc);
	  //	  nm = root_counter_(am, b);
	}

      } while ( s_ap == ZERO || s_am == CGAL::ZERO || np + 1 != nm );

      Interval ivl_upper =
	initial_ivl.split_at(iap, Interval::UPPER).second;
      push_front(ivl_upper, np, tag);

      Interval ivl_a = initial_ivl.collapse_to_bound(Interval::LOWER);
      push_front(ivl_a, tag);

      //      push_front(ap, b, np, tag);
      //      push_front(a, tag);
      return;
    }

    if ( s_a != CGAL::ZERO && s_b == CGAL::ZERO ) {
      Interval ib = initial_ivl.interval_around_bound(Interval::UPPER);

      Interval ibm =
	initial_ivl.split_at(ib, Interval::LOWER).second;
      Interval ibp(initial_ivl.upper_bound(), ib.upper_bound());
      //      Interval_NT bp = b + NT(1);
      //      Interval_NT bm = b_minus(a, b);

      Sign s_bp, s_bm;
      int np(0), nm(0);
      do {
	ibp = ibp.split().first;
	ibm = ibm.split().second;
	//	bp = Interval::midpoint(b, bp);
	//	bm = Interval::midpoint(b, bm);
	//	s_bp = sign_at_p(bp);
	//	s_bm = sign_at_p(bm);
	s_bp = ibp.apply_to_bound(sign_at_p, Interval::UPPER);
	s_bm = ibm.apply_to_bound(sign_at_p, Interval::LOWER);

	if ( s_bp != CGAL::ZERO ) {
	  np = ibp.apply_to_interval(rc);
	  //	  np = root_counter_(a, bp, check);
	}
	if ( s_bm != CGAL::ZERO ) {
	  nm = ibm.apply_to_interval(rc);
	  //	  nm = root_counter_(a, bm, check);
	}

	
      } while ( s_bp == CGAL::ZERO || s_bm == CGAL::ZERO ||
		np != nm + 1 );

      Interval ivl_b = initial_ivl.collapse_to_bound(Interval::UPPER);
      push_front(ivl_b, tag);

      Interval ivl_lower =
	initial_ivl.split_at(ibm, Interval::LOWER).first;
	//	Interval::create_from_endpoints(initial_ivl, Interval::LOWER,
	//					ibm, Interval::LOWER);
					
      push_front(ivl_lower, nm, tag);
      //      push_front(b, tag);
      //      push_front(a, bm, nm, tag);
      return;
    }

    if ( s_a == CGAL::ZERO && s_b == CGAL::ZERO ) {
      // both a and b are roots
      //      Interval_NT ap = a_plus(a, b);
      //      Interval_NT am = a - NT(1);
      //      Interval_NT bp  = b + NT(1);
      //      Interval_NT bm = b_minus(a, b);
      Interval ia =
	initial_ivl.interval_around_bound(Interval::LOWER);
      Interval ib =
	initial_ivl.interval_around_bound(Interval::UPPER);

      //      Interval_NT ap = Interval::midpoint(a, b);
      //      Interval_NT am = a - NT(1);
      //      Interval_NT bp  = b + NT(1);
      //      Interval_NT bm = Interval::midpoint(a, b);
      Interval iam(ia.lower_bound(), initial_ivl.lower_bound());
      Interval iap = initial_ivl.split_at(ia, Interval::UPPER).first;
      Interval ibm = initial_ivl.split_at(ib, Interval::LOWER).second;
      Interval ibp(initial_ivl.upper_bound(), ib.lower_bound());
	Interval::create_from_endpoints(initial_ivl, Interval::UPPER,
					ib, Interval::UPPER);

      Sign s_ap, s_am, s_bp, s_bm;
      int nbig(0), nsmall(0);

      do {
	iap = iap.split().first;
	iam = iam.split().second;
	s_ap = iap.apply_to_bound(sign_at_p, Interval::UPPER);
	s_am = iam.apply_to_bound(sign_at_p, Interval::LOWER);

	ibp = ibp.split().first;
	ibm = ibm.split().second;
	s_bp = ibp.apply_to_bound(sign_at_p, Interval::UPPER);
	s_bm = ibm.apply_to_bound(sign_at_p, Interval::LOWER);

	//	ap = Interval::midpoint(a, ap);
	//	am = Interval::midpoint(a, am);
	//	s_ap = sign_at_p(ap);
	//	s_am = sign_at_p(am);


	//	bp = Interval::midpoint(b, bp);
	//	bm = Interval::midpoint(b, bm);
	//	s_bp = sign_at_p(bp);
	//	s_bm = sign_at_p(am);


	if ( s_ap != CGAL::ZERO && s_bm != CGAL::ZERO ) {
	  Interval ismall(iap.upper_bound(), ibm.lower_bound());
	  nsmall = ismall.apply_to_interval(rc);
	  //	  nsmall = root_counter_(ap, bm, check);
	}
	if ( s_am != CGAL::ZERO && s_bp != CGAL::ZERO) {
	  Interval ibig(iam.lower_bound(), ibp.upper_bound());
	  nbig = ibig.apply_to_interval(rc);
	  //	  nbig = root_counter_(am, bp, check);
	}
      } while ( s_ap == CGAL::ZERO || s_am == CGAL::ZERO ||
		s_bp == CGAL::ZERO || s_bm == CGAL::ZERO ||
		nbig != nsmall + 2 );

      Interval ivl_b = initial_ivl.collapse_to_bound(Interval::UPPER);
      push_front(ivl_b, tag);

      Interval ivl_interior(iam.upper_bound(), ibm.lower_bound());
      //      push_front(ap, bm, sv_ap, sv_bm, tag);
      push_front(ivl_interior, nsmall, tag);

      Interval ivl_a = initial_ivl.collapse_to_bound(Interval::LOWER);
      push_front(ivl_a, tag);

      //      push_front(b, tag);
      //      push_front(ap, bm, nsmall, tag);
      //      push_front(a, tag);
      return;
    }

    bool this_line_should_not_have_been_reached = false;
    CGAL_assertion( this_line_should_not_have_been_reached );
  }
#endif
  //------------------------------------------------------------------

  void initialize_root_counter(Sturm_tag)
  {
    root_counter_ = Root_count(sseq, traits_);
  }
#if 0
  void initialize_root_counter(Bezier_tag)
  {
    initialize_root_counter(Descartes_tag());
  }

  void initialize_root_counter(Descartes_tag)
  {
    root_counter_ = Root_count(p);
  }
#endif

  Interval compute_initial_interval(Method_tag tag)
  {
    typedef typename Traits::Root_bound RBE;

    // MK: I MAY NEED TO CHANGE THIS SO THAT THE ROOT BOUND EVALUATOR
    //     TAKES A FUNCTION HANDLE; I MAY ALSO WANT TO USE THE
    //     INTERVAL VERSION OF THE FUNCTION TO DO THE COMPUTATIONS
    RBE rbe = traits_.root_bound_object(false, tag);

    typename RBE::result_type root_bound = rbe(p);

    Root B_lower = -root_bound;
    Root B_upper = root_bound;

    Root lower = start_;
    Root upper = end_;

    if ( start_ < B_lower ) {
      lower = B_lower;
    }

    if ( B_upper < end_ ) {
      upper = B_upper;
    }

    typedef typename Interval::NT INT;

    std::pair<double,double> ilower =
      CGAL::to_interval(lower.lower_bound());
    std::pair<double,double> iupper =
      CGAL::to_interval(upper.upper_bound());


    double low = std::min(ilower.first, iupper.first);
    double high = std::max(ilower.second, iupper.second);


    if ( CGAL::is_finite(low) && CGAL::is_finite(high) ) {
      Polynomial_assertion( low <= high );
      return Interval(INT(low), INT(high));
    }

    INT low_nt = std::min(lower.lower_bound(), upper.lower_bound());
    INT high_nt = std::max(lower.upper_bound(), upper.upper_bound());

    Polynomial_assertion( low_nt <= high_nt );

    return Interval(low_nt, high_nt);
  }

  // initialization method: unused check_if_simple
  void initialize(bool, Method_tag tag)
  {
    if ( p.is_zero() )  { return; }

    if ( end_ == -infinity() || start_ == infinity() ) {
      return;
    }

    Subdivision_tag  subdivision_tag;

    sseq = Standard_sequence(p);

#if 0
    normalize(p_simple, tag, success);
#endif

    Sign_at sign_at_p(p);
    
    initialize_root_counter( subdivision_tag );

    Interval initial_ivl = compute_initial_interval(tag);

    initialize_interval_list( initial_ivl, subdivision_tag, sign_at_p,
			      root_counter_ );
  }

public:
  //===================================
  // POSITIVE INFINITY
  //===================================
  static Root infinity() { return Root::infinity(); }


  Sturm_root_stack() {}

public:
  //==============
  // CONSTRUCTORS
  //==============
  Sturm_root_stack(const Polynomial& p,
		    const Root& start = -infinity(),
		    const Root& end = infinity(),
		    const Traits& tr = Traits(),
		    bool check_if_simple = true)
    : start_(start), end_(end), traits_(tr), root_counter_(),
      p(p), root_idx(0)
  {
    initialize( check_if_simple, Method_tag() );
    pop();
  }

  //===============
  // DESTRUCTOR
  //===============
  virtual ~Sturm_root_stack() {}

  //=============
  // METHODS
  //=============
protected:
#if 0
  const Root& top(Sturm_tag tag) const
  {
    return current_;
#if 0
    if ( i_list.size() == 0 ) {
      return infinity();
    } else {
      Interval ivl = i_list.front();
      return Root(ivl, p, sseq, root_idx+1);
    }
#endif
  }
#endif

  void pop(Sturm_tag tag) const
  {
    Sign_at sign_at_p(p);

    do {

      bool is_even_parity(false);
      do {
	if ( i_list.size() == 0 ) {
	  current_ = infinity();
	  return;
	} else {
	  subdivide( tag, sign_at_p, root_counter_ );

	  Interval ivl = i_list.front();
	  pop_front( tag );
	  // now I pass the simple polynomial and the gcd; I do not
	  // need to compute the simple polynomial in the solver?
	  // maybe for the descartes solver...
	  //	  cur = Root(ivl, p_simple, p_orig, ++root_idx);
	  current_ = Root(ivl, p, sseq, ++root_idx, traits_);

	} // end-if

	is_even_parity = current_.is_even_multiplicity();

	if ( DROP_EVEN_ROOTS && is_even_parity ) {
	  root_idx--;
	}

      } while ( DROP_EVEN_ROOTS && current_ != infinity() &&
		is_even_parity );

      if ( current_ <= start_ ) {
	root_idx--;
      }

    } while ( current_ != infinity() && current_ <= start_ );

#if 0
    if ( current_ == infinity() || current_ >= end_ ) {
      return infinity();
    }

    return cur;
#endif
  }


#if 0
  Root next_root(Bezier_tag) const
  {
    return next_root(Descartes_tag());
  }

  Root next_root(Descartes_tag tag) const
  {
    Root cur;
    Sign_at sign_at_p(p);

    do {
      bool is_even_parity(false);
      do {
	if ( i_list.size() == 0 ) {
	  cur = infinity();
	} else {
	  subdivide( tag, sign_at_p, root_counter_ );

	  if ( i_list.size() == 0 ) {
	    cur = infinity();
	  } else {
	    Interval ivl = i_list.front();
	    pop_front( tag );
	    //	    cur = Root(ivl, p_simple, p_orig, ++root_idx);
	    cur = Root(ivl, p, sseq[sseq.size() - 1], ++root_idx);
	  }
	} // end-if

	is_even_parity = cur.is_even_parity();

	if ( DROP_EVEN_ROOTS && is_even_parity ) {
	  root_idx--;
	}

      } while ( DROP_EVEN_ROOTS && cur != infinity() &&
		is_even_parity );

      if ( cur <= start_ ) {
	root_idx--;
      }

    } while ( cur != infinity() && cur <= start_ );

    if ( cur == infinity() || cur >= end_ ) {
      return infinity();
    }

    return cur;
  }
#endif
public:
  const Root& top() const
  {
    return current_;
  }

  bool empty() const {
    //    return ( i_list.size() == 0 || current() >= end_ );
    return ( top() >= end_ );
  }


  void pop() const
  {
    return pop( Subdivision_tag() );
  }

  const Root& start() const { return start_; }
  const Root& end()   const { return end_; }

  const Standard_sequence& standard_sequence() const { return sseq; }
  Standard_sequence&       standard_sequence() { return sseq; }

  const Root_count& root_count() const { return root_counter_; }
  Root_count&       root_count() { return root_counter_; }

  std::list<uint_pair>&     sv_list() { return s_variations; }
  std::list<unsigned int>&  nr_list() { return nroots; }

  Interval_container&        ivl_list() { return i_list; }
  const Interval_container&  ivl_list() const { return i_list; }

  int& root_index() { return root_idx; }

  template<class Stream>
  Stream& write(Stream& os) const
  {
    for (unsigned int i = 0; i < sseq.size(); i++) {
      os << sseq[i] << std::endl;
    }

    typename Interval_container::iterator ivl_it;
    typename std::list<uint_pair>::iterator nr_it;
    for (ivl_it = i_list.begin(), nr_it = s_variations.begin();
	 ivl_it != i_list.end(); ++ivl_it, ++nr_it) {
      os << "{";
      ivl_it->write(os);
      int nr = nr_it->first - nr_it->second;
      os << ", " << nr << "} ";
    }
    os << std::endl;
    return os;
  }


protected:
  Root                             start_, end_;
  Traits                           traits_;
  Root_count                       root_counter_;
  Polynomial                       p;
  Standard_sequence                sseq;
  mutable Root                     current_;
  mutable std::list<uint_pair>     s_variations;
  mutable std::list<unsigned int>  nroots;
  mutable Interval_container       i_list;
  mutable int                      root_idx;
};




template<class Stream, class T, class S, class M, bool D>
Stream& operator<<(Stream& os,
		   const Sturm_root_stack<T,S,M,D>& solver)
{
  return solver.write(os);
}



CGAL_POLYNOMIAL_END_NAMESPACE


#endif // CGAL_STURM_LAZY_SOLVER_H
