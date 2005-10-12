// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

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
*/
#ifdef NDEBUG
template <class Traits_t, bool SKIP_EVEN_ROOTS=true>
#else
template <class Traits_t, bool SKIP_EVEN_ROOTS=false>
#endif
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
      CGAL_exactness_assertion_code(typename Traits_t::Sign_at sar=k.sign_at_object(uf));
      CGAL_exactness_assertion_code(if (sar(lb)== CGAL::NEGATIVE){
				      std::cerr << "Invalid certificate with function " 
						<< uf << std::endl << "In interval from "
						<< lb << std::endl <<"to " << ub 
						<< std::endl;});
      CGAL_exactness_assertion(sar(lb)!= CGAL::NEGATIVE);
      
      CGAL::POLYNOMIAL::Sign sn= k.sign_between_roots_object(lb, solver_.top())(uf);
      if (sn == CGAL::NEGATIVE){
	root_=lb;
      } else {
	if (!solver_.empty()){
	  root_= solver_.top();
	  solver_.pop();
	  find_next_root();
	} else {
	  root_=std::numeric_limits<Root>::infinity();
	}
      }
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
      
      if (SKIP_EVEN_ROOTS || !rm_(root_)) {
	root_=solver_.top();
	solver_.pop();
	find_next_root();
      } else {
	if (one_even_){
	  root_=solver_.top();
	  solver_.pop();
	  one_even_=false;
	  find_next_root();
	} else {
	  one_even_=true;
	}
      }
    }

    bool empty() const {
      return root_==std::numeric_limits<Root>::infinity();
    }
  protected:

    void finish() {
      root_= Root::infinity();
    }
    void find_next_root() {
      if (!SKIP_EVEN_ROOTS) return;
      else {
	while (!solver_.empty() && rm_(root_)){
	  CGAL_KDS_LOG(LOG_LOTS, "skipping even root " << solver_.top() << std::endl);
	  root_=solver_.top();
	  solver_.pop();
	}
	if (rm_(root_)){
	  root_= std::numeric_limits<Root>::infinity();
	}
      }
    }

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
