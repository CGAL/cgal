#include <iostream>
#include <cassert>
#include <cstdlib>
#include <CGAL/Polynomial/Kernel.h>
#include <CGAL/Polynomial/Numeric_root_stack.h>
#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/Root_stack_default_traits.h>
//#include <CGAL/Polynomial/Upper_bound_root_stack_Descartes_traits.h>
//#include <CGAL/Polynomial/Upper_bound_root_stack.h>
#include <CGAL/Polynomial/internal/numeric_solvers.h>

#ifdef CGAL_POLYNOMIAL_USE_GSL
#include <CGAL/Polynomial/internal/GSL_numeric_solver.h>
#endif

#include "Check_solver.h"
bool verbose=true;
typedef CGAL_POLYNOMIAL_NS::Polynomial<double> Pd;
typedef CGAL_POLYNOMIAL_NS::Polynomial<CGAL::POLYNOMIAL::Default_field_nt> Pe;
typedef CGAL_POLYNOMIAL_NS::Root_stack_default_traits<Pd> Dt;

struct Interval_root_stack {
  typedef CGAL::POLYNOMIAL::Interval_polynomial Pi;
  typedef CGAL::POLYNOMIAL::Default_field_nt ENT;
  Interval_root_stack(){}
  typedef CGAL::POLYNOMIAL::Interval_nt Root;
  typedef Root Interval;
  typedef Dt Traits;
  Interval_root_stack(const Pd &p, Root lb, Root ub, Traits): lb_(lb),
							      ub_(ub) {
    CGAL_POLYNOMIAL_NS::Polynomial_converter<Pd, Pi, CGAL_POLYNOMIAL_NS::To_interval<double> > pc;
    CGAL_POLYNOMIAL_NS::Polynomial_converter<Pd, Pe, CGAL::NT_converter<double, ENT> > pce;
    p_=pc(p);
    pe_= pce(p);
    ok_=false;
    stack_.push_back(Root(std::numeric_limits<double>::infinity()));
    stack_.push_back(Root(lb.inf(), ub.sup()));
    refine();
  }

  const Root &top() const {
    assert(ok_);
    return stack_.back();
  }

  void pop() {
    if (stack_.size()>1) {
      assert(!stack_.empty());
      stack_.pop_back();
      ok_=false;
      if (stack_.size() >1) refine();
    }
  }

  bool empty() const {
    return stack_.size()==1;
  }

  static bool is_unknown_sign(Interval rv) {
    return rv.inf() <=0 && rv.sup() >=0;
  }

  bool is_odd(double lb, double ub) const {
    // put in filtered evaluation
    Interval lbi= p_(Interval(lb));
    Interval ubi= p_(Interval(ub));
    bool ret=false;
    if (is_unknown_sign(lbi) 
	|| is_unknown_sign(ubi)){
      ret= CGAL::sign(pe_(ENT(lb))) != CGAL::sign(pe_(ENT(ub)))
	|| CGAL::sign(pe_(ENT(lb))) == CGAL::ZERO;
    } else {
      ret= CGAL::sign(lbi) != CGAL::sign(ubi)
	|| CGAL::sign(lbi) == CGAL::ZERO;
    }
    return ret;
  }

  bool no_root(Interval r) const {
    Interval rv= p_(r);
    return !is_unknown_sign(rv);
  }

  static bool is_small(Interval r) {
    return r.sup()- r.inf() < .001;
  }

  double mp(double lb, double ub) const {
    if (lb == -std::numeric_limits<double>::infinity() 
	&& ub == std::numeric_limits<double>::infinity()){
      return 0;
    } else if (lb == -std::numeric_limits<double>::infinity()) {
      return ub-1000000;
    } else if (ub == std::numeric_limits<double>::infinity()){
      return lb + 1000000;
    } else {
      return (lb+ub)/2.0;
    }
  }

 

  void refine() {
    CGAL_POLYNOMIAL_NS::Interval_arithmetic_guard gd;
    assert(stack_.size()>1);
    while (stack_.size() >1 ){
      if (! is_small(stack_.back())) {
	//std::cout << "Splitting " << stack_.back() << std::endl;
	Interval b= stack_.back();
	stack_.pop_back();
	Interval fh= Interval(b.inf(), mp(b.inf(),b.sup()));
	Interval sh(fh.sup(), b.sup());
	stack_.push_back(sh);
	stack_.push_back(fh);
      }
      bool changed=false;
      do {
	changed=false;
	while(stack_.size() >1&&no_root(stack_.back()) ){
	  assert(!stack_.empty());
	  //std::cout << "Discarding " << stack_.back() << std::endl;
	  stack_.pop_back();
	  changed =true;
	}
	
	while (stack_.size() >1 && is_small(stack_.back()) 
	       && !is_odd(stack_.back().inf(), stack_.back().sup())) {
	  assert(!stack_.empty());
	  //std::cout << "Discarding " << stack_.back() << std::endl;
	  stack_.pop_back();
	  changed =true;
	}
      } while (changed);
      if (is_small(stack_.back())) break;
    }
   
    if (!empty()) {
      Interval b= stack_.back();
      assert(is_small(b));
      assert(is_odd(b.inf(),b.sup()));
    }
    ok_=true;
  }
  
  bool ok_;
  std::vector<Root> stack_;
  Root lb_, ub_;
  CGAL::POLYNOMIAL::Interval_polynomial p_;
  Pe pe_;
};


namespace std {
  template <>
  struct numeric_limits<CGAL_POLYNOMIAL_NS::Interval_nt>: public numeric_limits<double> {
    typedef numeric_limits<double> P;
    typedef CGAL_POLYNOMIAL_NS::Interval_nt  T;
    static const bool is_specialized = true;
    static T min BOOST_PREVENT_MACRO_SUBSTITUTION () throw() {return T((P::min)());}
    static T max BOOST_PREVENT_MACRO_SUBSTITUTION () throw() {return T((P::max)());}
    static T infinity() throw() {return P::infinity();}
  };
}


int main(int argc, char* argv[])
{
  //assert(std::numeric_limits<double>::has_infinity());
    if ( argc > 1 ) {
        int is_verbose = std::atoi(argv[1]);
        if ( is_verbose == 0 ) {
            verbose = false;
        } else verbose = true;
    }


#ifdef CGAL_USE_TNT
    {
        if (verbose) std::cout <<"JAMA______________________________________\n";
        else std::cout << "JAMA &\t";
        typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Dt,
            CGAL_POLYNOMIAL_NS::internal::JAMA_numeric_solver> NRE;
        typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
        K k;
        Check_solver<K > cg(k,verbose);
        cg.all();
        std::cout << std::endl;
    }
#endif
    if (0) {
        if (verbose) std::cout <<"Inter______________________________________\n";
        else std::cout << "Interval &\t";
	typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, Interval_root_stack> K;
        K k;
        Check_solver<K > cg(k,verbose);
        cg.all();
        std::cout << std::endl;
    }
#ifdef CGAL_POLYNOMIAL_USE_GSL
    {
        if (verbose) std::cout <<"GSL________________________________________\n";
        else std::cout << "GSL &\t";
        typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Dt, CGAL_POLYNOMIAL_NS::internal::GSL_numeric_solver> NRE;
        typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
        K k;
        Check_solver<K > cg(k,verbose);
        cg.all();
        std::cout << std::endl;
    }
#endif
    {
        if (verbose) std::cout <<"Turk______________________________________\n";
        else std::cout << "Turk &\t";
        typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Dt,
            CGAL_POLYNOMIAL_NS::internal::Turkowski_numeric_solver> NRE;
        typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
        K k;
        Check_solver<K > cg(k,verbose);
        cg.all();
        std::cout << std::endl;
    }

#ifdef CGAL_POLYNOMIAL_USE_GSL
    {
        if (verbose) std::cout <<"CleanGSL__________________________________\n";
        else std::cout << "CleanGSL &\t";
        typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Dt,
            CGAL_POLYNOMIAL_NS::internal::GSL_cleaned_numeric_solver> NRE;
        typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
        K k;
        Check_solver<K > cg(k,verbose);
        cg.cleaned();
        std::cout << std::endl;
    }
#endif

#ifdef CGAL_USE_TNT
    {
        if (verbose) std::cout <<"CleanJAMA________________________________\n";
        else std::cout << "CleanJAMA &\t";
        typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Dt,
            CGAL_POLYNOMIAL_NS::internal::JAMA_cleaned_numeric_solver> NRE;
        typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
        K k;
        Check_solver<K > cg(k,verbose);
        cg.cleaned();
        std::cout << std::endl;
    }
#endif
    {
        if (verbose) std::cout <<"CleanTurk________________________________\n";
        else std::cout << "CleanTurk &\t";
        typedef CGAL_POLYNOMIAL_NS::Numeric_root_stack<Dt,
            CGAL_POLYNOMIAL_NS::internal::Turkowski_cleaned_numeric_solver> NRE;
        typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
        K k;
        Check_solver<K > cg(k,verbose);
        cg.cleaned();
        std::cout << std::endl;
    }

/*{
  if (verbose) std::cout <<"Descartes__________________________________\n";
  else std::cout << "Descartes &\t";
  typedef CGAL_POLYNOMIAL_NS::Upper_bound_enumerator_Descartes_traits<Pd> Dt;
  typedef CGAL_POLYNOMIAL_NS::Upper_bound_root_enumerator<Dt> NRE;
  typedef CGAL_POLYNOMIAL_NS::Kernel<Pd, NRE> K;
  K k;
  Check_solver<K > cg(k,verbose);
  cg.all();
  std::cout << std::endl;
  }*/

    return 0;
}
