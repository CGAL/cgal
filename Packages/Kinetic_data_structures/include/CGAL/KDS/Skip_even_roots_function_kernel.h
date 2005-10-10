#ifndef CGAL_KDS_ROOT_ENUM_FK_H
#define CGAL_KDS_ROOT_ENUM_FK_H
#include <CGAL/KDS/basic.h>
#include <CGAL/Polynomial/internal/Explicit_root.h>
//#include <CGAL/Polynomial/Polynomial.h>
//#include <CGAL/Polynomial/polynomial_converters.h>



//template <class S> class Check_kds_solver;

CGAL_KDS_BEGIN_NAMESPACE;

//! A wrapper around solvers to make them better handle kinetic data structures
/*!  The wrapper modifies the functionality of the solvers from the
  kinetic kernel. It does the following 

  - verify that the certificate function is non-negative
  when created

  - handle the case when the certificate function fails
  immediately. The PolynomialKernel::Solver solvers
  handle open intervals. Sometimes, namely when the
  certificate function is 0 at the time of creation
  and has a negative derivative, then the beginning of
  the interval is actually the first root of interest.


  Using it is slightly akward due to the passing of unbound templates
  into the PolynomialKernel.

  See \example CGAL/KDS/define_exact_kds.h for an example of how to
  use this.
*/
template <class Traits_t>
struct Skip_even_roots_function_kernel: public Traits_t {

  class Root_stack: private Traits_t::Root_stack {
  private:
    typedef typename Traits_t::Root_stack Wrapped_solver;
    typedef typename Traits_t::Function Function;
    
  public:
    typedef typename Wrapped_solver::Root Root;
    typedef Traits_t Traits;
    //! Construct and check preconditions
    /*!
    
    */
    Root_stack(const Function &uf, const Root& lb, 
	       const Root& ub, const Traits_t& k): solver_(k.root_stack_object(uf, lb, ub)), 
						   one_even_(false),
						   rm_(k.is_even_multiplicity_object(uf)){
      CGAL_KDS_LOG(LOG_LOTS, "Solving " << uf << std::endl); //<< " at time " << lb << std::endl);
#if 0
      {
	Wrapped_solver sc(uf, lb, ub, k);
	CGAL_KDS_LOG(LOG_LOTS, "Found roots: ");
	while (!sc.finished()){
	  CGAL_KDS_LOG(LOG_LOTS, sc.top() << " ");
	  sc.pop();
	}
	CGAL_KDS_LOG(LOG_LOTS, std::endl);
      }
#endif
      CGAL_expensive_assertion_code(typename Traits_t::Sign_at sar=k.sign_at_object(uf));
      CGAL_expensive_assertion_code(if (sar(lb)== CGAL::NEGATIVE){std::cerr << "Invalid certificate with function " << uf << std::endl << "In interval from " << lb << std::endl <<"to " << ub << std::endl;});
      CGAL_expensive_assertion(sar(lb)!= CGAL::NEGATIVE);
      
      if (solver_.empty()){
	finish();
      } else {
	CGAL::POLYNOMIAL::Sign sn= k.sign_between_roots_object(lb, solver_.top())(uf);
	if (sn == CGAL::NEGATIVE){
	  root_=lb;
	  Wrapped_solver sp= k.root_stack_object(uf, lb,ub);
	  std::cerr << "First root is " <<  sp.top()<<std::endl;
	  std::cerr << "Degeneracy " << uf << " at " << lb << std::endl;
	} else {
	  // \todo fix hack
	  one_even_=true;
	  find_next_root();
	}
	CGAL_exactness_postcondition(root_ >=lb && root_<=ub);
      }
      CGAL_KDS_LOG(LOG_LOTS, "Solution is " << root_ << std::endl);
    }
  
    Root_stack(){}

    //! Drop even roots
    const Root& top() const {
      return root_;
    }
  
    //! Drop the front root.
    /*!  Even roots are skipped if debugging code is turned off. If
      debugging code is turned on, then the simulator needs to know
      about even roots so as not to ask for an audit during one of them.
    */
    void pop() {
      CGAL_precondition(!empty());
      if (!solver_.finished()){
	find_next_root();
      } else {
	finish();
      }
    }

    bool empty() const {
      return root_== Root::infinity();
    }
  protected:

    void finish() {
      root_= Root::infinity();
    }
    void find_next_root() {
#ifdef NDEBUG
      const bool debug=false;
#else
      const bool debug=true;
#endif

      if ( debug && !one_even_ && rm_(root_) ){
	CGAL_assertion(!one_even_);
	CGAL_KDS_LOG(LOG_LOTS, "Keeping front root" << std::endl);
	one_even_=true;
      } else {
	root_=solver_.top();
	solver_.pop();
	while (!debug && rm_(root_)){
	  if (solver_.empty()){
	    finish();
	    break;
	  } else {
	    CGAL_KDS_LOG(LOG_LOTS, "skipping even root " << root_ << std::endl);
	    assert(!solver_.empty());
	    root_=solver_.top();
	    solver_.pop();
	  }
	}
      }
    }
    /*
      static CGAL::Sign sign_between_roots(const Function &uf,
      const Root &lb, const Root &ub,  const Traits& k) {
      CGAL_expensive_precondition_code(Wrapped_solver ws(uf,lb, ub, k));
      CGAL_expensive_precondition(ws.finished());
      typename Traits::NT rat= rational_between_roots(lb, ub, k);
      CGAL_assertion(Root(rat) > lb && Root(rat) < ub);
      return sign(k.sign_at_object(uf)(rat));
      return 
      }

      static typename Traits::NT rational_between_roots(const Root &r0, const Root &r1, const Traits& ) {
      CGAL_precondition(r0<r1);
      typedef typename Traits::NT result_type;
      typedef std::pair<double, double> Ival;
      Ival i0= to_interval(r0);
      Ival i1= to_interval(r1);
      if (r0== -Root::infinity()) {
      if (r1== Root::infinity()){
      return result_type(0);
      } else {
      return result_type(i1.first-1);
      }
      } else if (r1== Root::infinity()){
      return i0.second+1;
      }
      
      if (i0.second < i1.first) {
      return result_type(.5*(i0.second+ i1.first));
      } else {
      //double d0= POLYNOMIAL_NS::to_double(r0);
      //double d1= POLYNOMIAL_NS::to_double(r1);
      double d0,d1;
      if (0){
      ++d0; ++d1;
      }
      Ival i0= to_interval(r0);
      Ival i1= to_interval(r1);
      if (i0.second < i1.first) {
      return result_type(.5*(i0.second+ i1.first));
      } else {
	  
      result_type lb = i0.first;
      result_type ub = i1.second;
      result_type mid;
      while (true){
      mid = .5*(lb+ub);
      if (Root(mid) <= r0) {
      lb= mid;
      } else if (Root(mid) >= r1){
      ub= mid;
      } else {
      CGAL_KDS_LOG(LOG_LOTS, "Between hard case " << r0 << "----" 
      << r1 << std::endl);
      CGAL_KDS_LOG(LOG_LOTS, "Return " << mid << std::endl);
      return mid;
      }
      }
      }
      }
      }

      template <class R>
      static CGAL::Sign sign_at_root(const Function &f, const R &r, const Traits &k) {
      std::pair<double, double> i= to_interval(r); // was CGAL
      if (i.first==i.second){
      double d= i.second;
      typename Traits::Sign_at sa= k.sign_at_object(f);
      CGAL::Sign sn= CGAL::sign(sa(typename Traits::NT(d)));
      return sn;
      } else {
      Wrapped_solver s= Wrapped_solver(f, typename Wrapped_solver::Root(i.first), 
      Wrapped_solver::Root::infinity(), k);
      while (!s.finished() && s.current() < r){
      s.advance();
      }
      if (s.current()==r){
      return CGAL::ZERO;
      }
      // now we know it is not a root
      //typename K::Root rr= k_.solver_object(p_, r).next_root();
      //POLYNOMIAL_NS::Sign sb;
      return sign_between_roots(f, r, s.current(), k);
      }
      }

      template <class RT>
      static CGAL::Sign sign_at_root(const Function &f, const POLYNOMIAL_NS::Explicit_root<RT> &r, const Traits &k){
      typedef  POLYNOMIAL_NS::Explicit_root<RT> R;
      typename R::Representation rep= r.representation();
      typedef  typename POLYNOMIAL_NS::Polynomial<typename R::Representation> Rep_poly;
      typename POLYNOMIAL_NS::Polynomial_converter<Function, Rep_poly> pc;
      return CGAL::sign(pc(f)(rep));
      }
    */
    Wrapped_solver solver_;
    Root root_;
    bool one_even_;
    typename Traits_t::Is_even_multiplicity rm_;
  };


  Root_stack root_stack_object(const typename Traits_t::Function &f,
			       const typename Traits_t::Root &lb 
			       = -std::numeric_limits<typename Traits_t::Root>::infinity(),
			       const typename Traits_t::Root &ub 
			       = std::numeric_limits<typename Traits_t::Root>::infinity()) {
    return Root_stack(f, lb, ub, *this);
  }
};



CGAL_KDS_END_NAMESPACE;

#endif
