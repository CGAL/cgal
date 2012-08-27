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

#ifndef CGAL_POLYNOMIAL_ROOT_BOUND_SOLVER_CORE_H
#define CGAL_POLYNOMIAL_ROOT_BOUND_SOLVER_CORE_H

#include <CGAL/Polynomial/basic.h>
//#include <CGAL/Polynomial/Tools/Isolating_interval.h>
//#include <CGAL/Polynomial/Tools/Simple_interval_root.h>
#include <CGAL/Polynomial/internal/Descartes_root_count.h>
#include <CGAL/Polynomial/polynomial_converters.h>
#include <CGAL/Polynomial/Interval_polynomial.h>
#include <vector>

#define CGAL_DSPRINT(x)

namespace CGAL { namespace POLYNOMIAL {
//! A Solver which uses
template <class Traits_t>
class Upper_bound_root_stack
{
protected:
  typedef Upper_bound_root_stack<Traits_t> This;

  typedef internal::Descartes_root_count Root_count;

  struct Interval_info;
  typedef typename Traits_t::Isolating_interval Interval;

  typedef std::vector<Interval_info> Intervals;
  typedef typename Traits_t::Sturm_root_count Sturm_root_counter;
  typedef typename Traits_t::Function Polynomial;
  typedef typename Traits_t::FT FT;
  enum Interval_classification {GOOD, EMPTY, SOME};
public:
  typedef Traits_t Traits;
  typedef typename Traits_t::Root Root;
  /*struct Root: public Traits_t::Root {
    Root(){}
    template <class TNT>
    Root(TNT nt): Traits_t::Root(nt){}
    Root(const NT &nt, bool is_odd, Traits_t k): Traits_t::Root(nt, is_odd, k){}
    Root(const Interval &nt, bool is_odd=true, Traits_t k= Traits_t()):
    Traits_t::Root(nt, is_odd, k){}
    Root(const Interval &ii, const typename Traits_t::Function &sa,
    Sign slb, Sign sub,
    Traits_t k): Traits_t::Root(ii, sa, slb, sub, k){}
    };
    friend class Root;*/

  Upper_bound_root_stack(): cur_(Root::infinity()) {
  };
  Upper_bound_root_stack(const Polynomial &f,
			 const Root &lb,
			 const Root &ub,
			 const Traits &tr): kernel_(tr), f_(f),
					    lb_(lb), ub_(ub),
					    rc_(kernel_.root_count_object(f_)),
					    has_ss_(false) {
    //std::cout << f << std::endl;
    initialize();

    //do {
    isolate_top();
    //} while (!(cur_ >lb_));

    //compute_estimate(lb);
   
  };
  Upper_bound_root_stack(const Polynomial &f,
			 const Root &lb,
			 const Root &ub,
			 const Traits &tr,
			 bool dont_isolate): kernel_(tr), f_(f),
					     lb_(lb), ub_(ub),
					     rc_(kernel_.root_count_object(f_)),
					     has_ss_(false) {
    initialize();
  };

  const Root& top() const
  {
    CGAL_Polynomial_precondition(cur_ != std::numeric_limits<Root>::quiet_NaN());
    return cur_;
  }

  void pop() {
    isolate_top();
  }


  /*double estimate() const {
    return estimate_;
    }*/

  bool empty() const
  {
    return cur_==std::numeric_limits<Root>::infinity();
  }

  /*Isolating_interval spanned_interval() const {
    CGAL_precondition(!empty());
    if (intervals_.empty()){
    return cur_.isolating_interval();
    } else {
    Isolating_interval lb;
    if (cur_ != Root()){
    return lb= cur_.isolating_interval();
    } else {
    lb= intervals_.back().interval;
    }

    Isolating_iterval ub= intervals_.front().interval;
    return lb| ub;
    }
    }*/

  // the following methods should be protected, but...
  bool top_is_isolated() const
  {
    return cur_ != std::numeric_limits<Root>::quiet_NaN();

    /*intervals_.back().interval.is_singular()
      || intervals_.back().num_roots.is_single()
      && intervals_.back().left_sign != ZERO
      && intervals_.back().right_sign != ZERO
      && intervals_.back().right_sign != intervals_.back().left_sign;*/
  }

  bool refine_top(FT ub) {
    while (!intervals_.empty() && intervals_.back().interval.contains(ub) < ub) {
      Interval_classification c=classify_top();
      if (c==GOOD) {
	cur_=make_next_root();
	if (cur_ <= lb_) {
	  cur_= std::numeric_limits<Root>::quiet_NaN();
	  continue;
	}
	else {
	  return true;
	}
      }
      else if (c== EMPTY) {
	intervals_.pop_back();
      }
      else {
	subdivide_top();
      }
    }
    return false;
  }

  bool refine_bottom(double lb) {
    CGAL_Polynomial_assertion(0);
    return false;
  }

  void pop_no_isolate() {
    if (intervals_.empty()) {
      cur_= std::numeric_limits<Root>::infinity();
    }
    else {
      cur_= std::numeric_limits<Root>::quiet_NaN();
    }
  }

  /*Isolating_interval interval() const {

  }*/

  std::pair<double, double> double_interval() const
  {
    //CGAL_Polynomial_precondition(!empty());
    if (empty()) return std::pair<double, double>(std::numeric_limits<double>::infinity(),
						  std::numeric_limits<double>::infinity());
    double lb;
    if (top_is_isolated()) {
      lb= cur_.double_interval(std::numeric_limits<double>::infinity()).first;
    }
    else if (!intervals_.empty()) {
      lb= intervals_.back().interval.double_interval().first;
    }
    else {
      return std::pair<double, double>(std::numeric_limits<double>::infinity(),
				       std::numeric_limits<double>::infinity());
    }
    double ub;
    if (intervals_.empty()) {
      ub= cur_.double_interval(std::numeric_limits<double>::infinity()).second;
    }
    else {
      ub= intervals_.front().interval.double_interval().second;
    }
    return std::pair<double, double>(lb, ub);
  }

  //! establish the invariant
  /*!
    The front interval contains only one root.
  */
  void isolate_top() {
    //write_intervals(std::cout);
    do {
      while (!intervals_.empty()) {
	Interval_classification cl= classify_top();
	if (cl== EMPTY) intervals_.pop_back();
	else if (cl== SOME) {
	  subdivide_top();
	}
	else {
	  break;
	}
      }                                 // while (the first has more than one interval);
      //std::cout << "Initialized.\n";
      //if (!intervals_empty()){
      cur_= make_next_root();
    } while (cur_<= lb_);
    // }
  }

  /*void compute_estimate(const Root &lb) {
    Interval_arithmetic_guard gd;
    Polynomial_converter<Polynomial, Interval_polynomial,  To_interval<FT> > pc;
    Interval_polynomial ip= pc(f_);
    Interval_nt iv(CGAL::to_interval(lb_));
    if (iv.inf()==iv.sup()) {
      iv= std::make_pair(iv.inf(), iv.inf()+.1);
    }
    while (ip(iv).inf() >0 && iv.sup() != std::numeric_limits<double>::infinity()) {
      iv= std::pair<double, double>(iv.inf(), iv.inf()+ 2*(iv.sup()-iv.inf()));
    }
    
    if (iv.sup() !=  std::numeric_limits<double>::infinity()) {
      double top= iv.sup();
      iv= std::make_pair(iv.inf(), iv.inf()+ (iv.sup()-iv.inf())/2.0);
      for (int i=0; i< 4; ++i){
	double mp= (top+iv.sup())/2.0;
	Interval_nt ivt= std::make_pair(iv.inf(), mp);
	if ( ip(ivt).inf() >0) {
	  iv= ivt;
	} else {
	  top= ivt.sup();
	}
      }
    } else {
      ++filtered__;
    }

    estimate_= iv.inf()+ (iv.sup()-iv.inf())/2.0;
  }*/

protected:

 

  Root make_next_root() {
    if (intervals_.empty()) {
      //std::cout << "Out of intervals.\n";
      return Root::infinity();
    }

    CGAL_precondition(!intervals_.back().num_roots.is_zero());
    Root ret;

    if (intervals_.back().interval.is_singular()) ret= Root(intervals_.back().interval,
							    intervals_.back().num_roots.is_odd());
    else ret= Root(intervals_.back().interval,f_, intervals_.back().left_sign,
		   intervals_.back().right_sign,  kernel_);
    CGAL_DSPRINT(std::cout << "Trying root " << ret << std::endl);
    CGAL_DSPRINT(std::cout << "Root count was " << intervals_.back().num_roots << std::endl);
    if ( !(ub_> ret)) {
      CGAL_DSPRINT(std::cout << "Rejection due to >=" << ub_ << std::endl);
      intervals_.clear();
      return Root::infinity();
    }
    else {
      intervals_.pop_back();
      return ret;
    }
  }

  void initialize() {
    CGAL_Polynomial_expensive_precondition(lb_<ub_);

    if (f_.is_constant()) return;

    typename Traits::Root_bound rbe= kernel_.root_bound_object();
    FT rb; if (lb_==-Root::infinity() || ub_== Root::infinity()) rb= rbe(f_);
    Interval lbi, ubi;
    if (lb_ == -Root::infinity()) {
      lbi= Interval(-rb);
    }
    else {
      lbi= lb_.isolating_interval_object();    //power_of_two(lb_.interval().lower_bound());
    }
    if (ub_== Root::infinity()) {
      ubi= Interval(rb);
    }
    else {
      ubi= ub_.isolating_interval_object();
    }
    Interval ii= lbi || ubi;

    intervals_.push_back(Interval_info(ii,// skip rc
				       ii.apply_to_endpoint(kernel_.sign_at_object(f_),
							    Interval::LOWER),
				       ii.apply_to_endpoint(kernel_.sign_at_object(f_), Interval::UPPER)));
    CGAL_DSPRINT(std::cout << "Initial interval is " << ii <<std::endl);
  }

  void subdivide_top() {

    Interval_info ii= intervals_.back();
    CGAL_DSPRINT(std::cout << "Investigating " << ii.interval << std::endl);
    intervals_.pop_back();

    const Interval uh= ii.interval.second_half();
    const Interval mi= uh.lower_endpoint_interval();
    const Interval lh= ii.interval.first_half();
    CGAL::Sign sm= uh.apply_to_endpoint(kernel_.sign_at_object(f_), Interval::LOWER);

    bool mid_is_root=(sm == CGAL::ZERO);

    CGAL_DSPRINT(std::cout << "Produced u interval of " << uh << std::endl);
    intervals_.push_back(Interval_info(uh, sm, ii.right_sign));
    //}
    if (mid_is_root) {
      intervals_.push_back(Interval_info(mi, CGAL::ZERO, CGAL::ZERO));
    }

    intervals_.push_back(Interval_info(lh, ii.left_sign, sm));

    CGAL_DSPRINT(write_intervals(std::cout));
  }

  /*void subdivide_top(double v){
    Interval_info ii= intervals_.back();
    CGAL_DSPRINT(std::cout << "Investigating " << ii.interval << std::endl);
    intervals_.pop_back();

    }*/

  //return true if this pops
  unsigned int sturm_top() {
    if (!has_ss_) {
      ss_= kernel_.Sturm_root_count_object(f_);
      has_ss_=true;
    }
    CGAL_DSPRINT(std::cout << "Computing sturm for " << intervals_.back().interval << std::endl);
    //std::cout << "Computing sturm for " << intervals_.back().interval << std::endl;
    return intervals_.back().interval.apply_to_interval(ss_);
  }

  void upperbound_top() {
    intervals_.back().num_roots=compute_root_count(intervals_.back().interval,
						   intervals_.back().left_sign,
						   intervals_.back().right_sign);    //ii.apply_to_interval(rc_);
  }

  void count_singular_top() {
    CGAL_Polynomial_assertion(intervals_.back().num_roots.is_unknown());
    int deg=intervals_.back().interval.apply_to_endpoint(kernel_.multiplicity_object(f_), Interval::LOWER);
    intervals_.back().num_roots=deg;
  }

  Interval_classification classify_top() {
    CGAL_Polynomial_assertion(intervals_.back().left_sign
			      == intervals_.back().interval.apply_to_endpoint(kernel_.sign_at_object(f_), Interval::LOWER));
    CGAL_Polynomial_assertion(intervals_.back().right_sign
			      == intervals_.back().interval.apply_to_endpoint(kernel_.sign_at_object(f_), Interval::UPPER));

    // compute a root count
    if (intervals_.back().interval.is_singular()) {
      count_singular_top();
      //std::cout << "Breaking in singular." << std::endl;
      return GOOD;
    }
    else if (intervals_.back().interval.approximate_width() <= min_interval_width()) {
      unsigned int ct= sturm_top();
      if (ct == 0) {
	//intervals_.pop_back();
	return EMPTY;
      }
      else if (ct==1) {
	intervals_.back().num_roots =1;
	if (intervals_.back().num_roots.is_single()
	    && intervals_.back().left_sign != ZERO
	    && intervals_.back().right_sign != ZERO) return GOOD;
	else return SOME;
      }
      else {
	intervals_.back().num_roots=Root_count(ct);
	return SOME;
      }
    }
    else  if (intervals_.back().num_roots.is_unknown()) {
      upperbound_top();

      if (intervals_.back().num_roots.is_zero()) {
	//intervals_.pop_back();
	return EMPTY;
      }

      if (intervals_.back().num_roots.is_single()
	  && intervals_.back().left_sign != ZERO
	  && intervals_.back().right_sign != ZERO) {
	CGAL_assertion( intervals_.back().left_sign != intervals_.back().right_sign);
	//std::cout << "Breaking in normal." << std::endl;
	return GOOD;
      }
      else {
	return SOME;
      }
    }
    CGAL_Polynomial_postcondition(0);
    return SOME;
  }

  Root_count compute_root_count(const Interval &ii,
				CGAL::Sign sl,
				CGAL::Sign sr) {
    return Root_count(ii.apply_to_interval(rc_, sl, sr));
  }

  void write_intervals(std::ostream &out) {
    for (unsigned int i=0; i< intervals_.size(); ++i) {
      out << "(" << intervals_[i].interval << ", " << intervals_[i].left_sign << intervals_[i].right_sign << ")";
    }
    std::cout << std::endl;
  }
  static double min_interval_width() {
    return .000001;
  }

  //! Just return small intervals if doubles are used
  template <class NTT>
  static bool is_small(const NTT &, double) {
    return false;
  }

  static bool is_small(double, double wid) {
    return wid < .00001;
  }

  struct Interval_info
  {
    Root_count num_roots;
    CGAL::Sign left_sign;
    CGAL::Sign right_sign;
    Interval_info(Interval in, /*Root_count num,*/ CGAL::Sign ls,
		  CGAL::Sign rs):
      num_roots(Root_count::UNKNOWN), left_sign(ls), right_sign(rs), interval(in){}
    Interval interval;
  };

  Traits kernel_;
  Polynomial f_;
  Intervals intervals_;
  Root lb_, ub_;
  Root cur_;
  typename Traits::Root_count rc_;
  Sturm_root_counter ss_;
  bool has_ss_;
};
} } //namespace CGAL::POLYNOMIAL

#endif
