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

#ifndef CGAL_POLYNOMIAL_CORE_SOLVER_H
#define CGAL_POLYNOMIAL_CORE_SOLVER_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Polynomial/internal/CORE_polynomial.h>
#include <CGAL/Polynomial/internal/Root_stack_traits_base.h>
#include <CGAL/CORE_BigInt.h>

#include <iostream>

/*namespace CGAL {
double to_double(const CORE::BigInt &bi)
{
  return bi.doubleValue();
  }


  } //namespace CGAL*/

namespace CGAL { namespace POLYNOMIAL {


class CORE_Expr_root_stack
{
protected:
  typedef CORE_Expr_root_stack This;
  //typedef CORE::Poly<CORE::BigInt> BIP;
public:

  typedef internal::CORE_polynomial Function;
  typedef Function::NT Coef;
  typedef CORE::Sturm<Coef> CORE_Sturm;
  

  struct Traits: public internal::Root_stack_traits_base<Function> {
    
  };

  typedef CORE::Expr Root;

  CORE_Expr_root_stack(const Function &f,
		       const Root &lb,
		       const Root &ub,
		       const Traits &tr): f_(f), ub_(ub),  cur_valid_(false), tr_(tr), one_even_left_(false){
    initialize(lb);
  }

  CORE_Expr_root_stack():  num_roots_(0){}

  const Root& top() const
  {
    CGAL_precondition(!empty());
    if(!cur_valid_) make_cur_valid();
    return cur_;
  }
  void pop() {
    if (!one_even_left_) {
      --num_roots_;
      cur_valid_=false;
      CGAL_precondition(num_roots_>=0);
    } else {
      one_even_left_=false;
    }
  }

  bool empty() const
  {
    if (num_roots_ != 0 && !cur_valid_) {
      make_cur_valid();
    }
    return num_roots_==0;
  }


  std::ostream &write(std::ostream &out) const {
    return out << f_ << ": " << cur_ << std::endl;
  }
protected:
  void make_cur_valid() const {
    CGAL_precondition(!cur_valid_);
    if (num_roots_==0) {
      no_roots();
    } else {
      make_root();
      enforce_upper_bound();
    }
  }

  Function f_;
  CORE_Sturm sturm_;
  Root  ub_;
  mutable Root cur_;
  mutable bool cur_valid_;
  mutable int num_roots_;
  mutable CORE::BigFloat bflb_, bfub_;
  Traits tr_;
  mutable bool one_even_left_;
  mutable int offset_in_interval_;
  mutable CGAL::Sign last_sign_;

  void initialize(const Root& lb) {
    if (f_.degree()<=0) {
      no_roots();
      return;
    } else {
      //std::cout <<"solving " << f_ << std::endl;
      f_.contract();
      //std::cout << f_.core_polynomial() << std::endl;
      sturm_= CORE_Sturm(f_.core_polynomial()/*, false*/); //BigInt to BigRat
      
      offset_in_interval_=0;
      //CORE::BigFloat bflb, bfub;
      
      /*if (lb == -std::numeric_limits<Root>::infinity()){
	bflb_= -f_.core_polynomial().CauchyUpperBound();
	} else {*/
      CORE::BigFloat offset(.5);
      CGAL_postcondition(offset.isExact());
      bflb_= bf_lower_bound(lb);
      CGAL_postcondition(bflb_.isExact());
      do {
	bflb_ -= offset; // hack to get around assuming core is consistent with 0 endpoint
	last_sign_=CGAL::sign(f_.core_polynomial().eval(bflb_));
      } while (last_sign_==CGAL::ZERO);
      CGAL_postcondition(bflb_.isExact());

      bfub_= bf_upper_bound(ub_);
      CGAL_precondition(bflb_ < bfub_);
      
      
      //std::cout << " in interval: " << bflb_ << " " << bfub_ << std::endl;
      num_roots_= sturm_.numberOfRoots(bflb_, bfub_);
      //std::cout << "nr= " << num_roots_ << std::endl;
      //CORE::Expr testr;
      ++num_roots_;
      do {
	--num_roots_;
	if ( num_roots_ == 0) {
	  no_roots();
	  return;
	}
	make_root();
      } while (cur_ < lb);
      //make_cur_root(testr);
    }
    //std::cout << "There are " << num_roots_ << " roots.\n";
    enforce_upper_bound();
  }

  void enforce_upper_bound() const {
    if (cur_ < ub_) return;
    else {
      //std::cout << "Rejected " << cur_ << std::endl;
      no_roots();
    }
  }

  void make_root() const {
    CGAL_precondition(num_roots_!=0);
    //std::cout << "making root: " << CGAL::to_double(bflb_) << " " << CGAL::to_double(bfub_) << std::endl;
  
    CORE::BFInterval bfi;
    bfi= sturm_.isolateRoot(1+offset_in_interval_, bflb_, bfub_);
   
    //std::cout << "got: " << CGAL::to_double(bfi.first) << " " << CGAL::to_double(bfi.second) << std::endl;
    //int nr= sturm_.numberOfRoots(bfi.first, bfi.second);
    //int nr=1;
    CGAL::Sign cur_sign=CGAL::sign(f_.core_polynomial().eval(bfi.second));
    if (cur_sign==0) {
      Traits::Sign_after sa= tr_.sign_after_object();
      cur_sign= sa(f_, bfi.second);
      ++offset_in_interval_;
    } else {
      offset_in_interval_=0;
      bflb_= bfi.second;
    }
    CGAL_precondition(cur_sign!= CGAL::ZERO);
    CGAL_precondition(last_sign_ != CGAL::ZERO);
    if (last_sign_== cur_sign) {
      one_even_left_=true;
      //std::cout << "it is even" << std::endl;
    } else {
      one_even_left_=false;
    }
    last_sign_= cur_sign;
    //std::cout << nr << " " << bfi.first << " " << bfi.second <<  std::endl;
    cur_= CORE::Expr(f_.core_polynomial(), bfi);
    cur_valid_=true;
    //cur_ =Root(e/*/f_.scale()*/, nr);
    //std::cout << "root= " << cur_ <<  " " << e << std::endl;
  }

  void no_roots() const {
    //ub_= CORE::Expr(0);
    cur_= infinity<Root>();
    num_roots_=0;
    one_even_left_=false;
    cur_valid_=false;
  }

  /*void initialize_counters(const Root &lb) {
    std::cout << "Computing strum of " << poly_ << "..." << std::flush;
    CORE_Sturm sturm(poly_);
    std::cout << "done." << std::endl;
    num_roots_=0;
    CGAL_assertion(-ub_ != infinity<Root>());
    num_roots_= sturm.numberOfRoots();
    if (lb== -infinity<Root>() && ub_== infinity<Root>()) {
      counter_=0;
    }
    else if (ub_ == infinity<Root>()) {
      std::cout << bf_lower_bound(lb.representation()) << std::endl;
      //num_roots_= sturm.numberOfRootsAbove(bf_lower_bound(lb.representation()));
      counter_ = sturm.numberOfRootsBelow(bf_lower_bound(lb.representation()));
    }
    else if (lb == infinity<Root>()) {
      //num_roots_= sturm.numberOfRootsBelow(bf_upper_bound(ub_.representation()));
      counter_ = 0;
    }
    else {
      counter_= sturm.numberOfRootsBelow(bf_lower_bound(lb.representation()));
      
    }
    }*/

  //! There are probably better ways of doing this
  Coef bf_lower_bound(const CORE::Expr &rt) const
  {
    machine_double lb, ub;
    rt.doubleInterval(lb, ub);
    return Coef(lb);
  }

  //! There are probably better ways of doing this
  Coef bf_upper_bound(const CORE::Expr &rt) const
  {
    machine_double lb, ub;
    rt.doubleInterval(lb, ub);
    return Coef(ub);
  }
};

inline std::ostream &operator<<(std::ostream &out, const CORE_Expr_root_stack &o) {
  return o.write(out);
}

} } //namespace CGAL::POLYNOMIAL;
#endif
