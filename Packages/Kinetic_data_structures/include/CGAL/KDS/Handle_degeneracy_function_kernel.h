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

#ifndef CGAL_KDS_ROOT_DEGEN_FK_H
#define CGAL_KDS_ROOT_DEGEN_FK_H
#include <CGAL/KDS/basic.h>




CGAL_KDS_BEGIN_NAMESPACE;


template <class Traits_t>
struct Handle_degeneracy_function_kernel: public Traits_t {

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
	       const Root& ub, const Traits_t& k): solver_(k.root_stack_object(uf, lb, ub)){
      
      CGAL::POLYNOMIAL::Sign sn= k.sign_between_roots_object(lb, solver_.top())(uf);
      if (sn == CGAL::NEGATIVE){
	extra_root_=lb;
      } else {
	extra_root_=std::numeric_limits<Root>::infinity();
      }
    }
  
    Root_stack(){}

    //! Drop even roots
    const Root& top() const {
      if (extra_root_== std::numeric_limits<Root>::infinity())
	return solver_.top();
      else return extra_root_;
    }
  
    void pop() {
      if (extra_root_!= std::numeric_limits<Root>::infinity()){
	extra_root_=std::numeric_limits<Root>::infinity();
      } else {
	solver_.pop();
      }
    }

    bool empty() const {
      return extra_root_==std::numeric_limits<Root>::infinity() && solver_.empty();
    }
  protected:
    Wrapped_solver solver_;
    Root extra_root_;
  };


  Root_stack root_stack_object(const typename Traits_t::Function &f,
			       const typename Traits_t::Root &lb,
			       const typename Traits_t::Root &ub) {
    return Root_stack(f, lb, ub, *this);
  }
};



CGAL_KDS_END_NAMESPACE;

#endif
