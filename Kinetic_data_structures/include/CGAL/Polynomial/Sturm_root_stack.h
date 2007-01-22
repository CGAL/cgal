// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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

#ifndef CGAL_STURM_LAZY_SOLVER_H
#define CGAL_STURM_LAZY_SOLVER_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Simple_interval_root.h>
#include <CGAL/Polynomial/internal/Sturm_isolating_interval.h>


#include <list>

#include <CGAL/Polynomial/Kernel.h>

CGAL_POLYNOMIAL_BEGIN_NAMESPACE

//================================================================
//================================================================
// ************** the lazy Sturm solver **************************
//================================================================
//================================================================

#define CGAL_KINETIC_STURM_DEBUG(x) 
//std::cout << x << std::endl;
#define CGAL_KINETIC_STURM_DEBUG_WRITE(x) 
//x << std::endl;

template<class T>
class Sturm_root_stack
{
public:
  typedef T                             Traits;
  typedef internal::Simple_interval_root<Traits>  Root;
  typedef typename Traits::Function     Polynomial;
  typedef typename Traits::Isolating_interval     Interval;
  typedef typename Traits::FT           NT;
  typedef typename Traits::Sign_at      Sign_at;
  typedef typename Traits::Standard_sequence Standard_sequence;

protected:

  typedef typename Traits::Root_count            Root_count;
  typedef Sturm_root_stack<T>  Self;

  struct ID {
    std::pair<NT,NT> interval_;
    unsigned int lbc_, ubc_;
    ID(const std::pair<NT,NT> &iv,
       unsigned int lc, unsigned int uc):  interval_(iv), lbc_(lc), ubc_(uc){}
    std::ostream &write(std::ostream &out) const {
      if (interval_.first == interval_.second) out << interval_.first;
      else out << "(" << CGAL::to_double(interval_.first) << ", " << CGAL::to_double(interval_.second) << ")";
      out << "(" << lbc_ << "," << ubc_ << ")";
      return out;
    }
  };


protected:

  //------------------------------------------------------------------
  // the method for subdivision; it guarantees that the first interval
  // contains only one root
  void subdivide() {
    if (intervals_.empty() ) { return; }

    while ( intervals_.back().lbc_ - intervals_.back().ubc_!=1  ) {
      CGAL_KINETIC_STURM_DEBUG("Next interval has " << intervals_.back().lbc_ - intervals_.back().ubc_ << " roots");
      CGAL_precondition( intervals_.back().lbc_ - intervals_.back().ubc_ >1);
      ID ivl = intervals_.back();
      intervals_.pop_back();
      NT mp= (ivl.interval_.first+ ivl.interval_.second)/NT(2.0);

      CGAL::Sign mps= sign_at(p_, mp);
      if (mps==0) {
	CGAL_KINETIC_STURM_DEBUG("Hit root at " << CGAL::to_double(mp));
	NT offset= (mp-ivl.interval_.first)/NT(10.0);
	
	unsigned int mlc, muc;
	do {
	  // OK, a bit silly
	  offset= offset/NT(2.0);
	  mlc= root_counter_(mp-offset);
	  muc= root_counter_(mp+offset);
	  CGAL_precondition(mlc > muc);
	} while ( mlc-muc != 1);
	CGAL_KINETIC_STURM_DEBUG("Settled on offset of " << CGAL::to_double(offset));
	CGAL_assertion(ivl.ubc_<= muc);
	if (muc- ivl.ubc_ >0) {
	  intervals_.push_back(ID(std::make_pair(mp+offset, ivl.interval_.second), muc, ivl.ubc_));
	  CGAL_KINETIC_STURM_DEBUG("Adding interval of ");
	  CGAL_KINETIC_STURM_DEBUG_WRITE(intervals_.back().write(std::cout)); 
	}
	CGAL_assertion(mlc-muc == 1);
	intervals_.push_back(ID(std::make_pair(mp, mp), mlc, muc));
	CGAL_KINETIC_STURM_DEBUG("Adding interval of ");
	CGAL_KINETIC_STURM_DEBUG_WRITE(intervals_.back().write(std::cout)); 
	CGAL_assertion(ivl.lbc_ >= mlc);
	if (ivl.lbc_ - mlc >0) {
	  intervals_.push_back(ID(std::make_pair(ivl.interval_.first, mp-offset),ivl.lbc_, mlc));
	  CGAL_KINETIC_STURM_DEBUG("Adding interval of ");
	  CGAL_KINETIC_STURM_DEBUG_WRITE(intervals_.back().write(std::cout)); 
	}
      } else {
	unsigned int mpc= root_counter_(mp);
	CGAL_precondition( ivl.ubc_ <=mpc);
	if (mpc-ivl.ubc_!= 0) {
	  intervals_.push_back(ID(std::make_pair(mp, ivl.interval_.second), mpc, ivl.ubc_));
	  CGAL_KINETIC_STURM_DEBUG(std::cout << "Adding interval of ");
	  CGAL_KINETIC_STURM_DEBUG_WRITE(intervals_.back().write(std::cout)); 
	}
	CGAL_precondition(mpc <= ivl.lbc_);
	if ( ivl.lbc_-mpc != 0) {
	  intervals_.push_back(ID(std::make_pair( ivl.interval_.first, mp), ivl.lbc_, mpc));
	  CGAL_KINETIC_STURM_DEBUG(std::cout << "Adding interval of ");
	  CGAL_KINETIC_STURM_DEBUG_WRITE(intervals_.back().write(std::cout)); 
	}
      }
    }                                     // end-while
  }                                         // end subdivide method

  void initialize() {
    if ( p_.is_zero() )  { 
      CGAL_KINETIC_STURM_DEBUG("Zero polynomial ");
      return; 
    }
    if (start_==finish_) {
      CGAL_KINETIC_STURM_DEBUG("Empty interval ");
      return;
    }

    Root ninf= -Root::infinity();
    
    if ( finish_ == ninf || start_ == Root::infinity() ) {
      CGAL_KINETIC_STURM_DEBUG("Empty interval " << start_ << " to " << finish_);
      return;
    }

    sseq_ = Standard_sequence(p_);
    CGAL_postcondition(sseq_[sseq_.size()-1].degree() > -1);
    if (sseq_[sseq_.size()-1].degree() > 0) {
      CGAL_KINETIC_STURM_DEBUG("Non-square free: " << sseq_[sseq_.size()-1]);
      non_square_free_part_= sseq_[sseq_.size()-1];
      typename Traits::Quotient quo= traits_.quotient_object();
      p_= quo(p_,non_square_free_part_);
      sseq_ = Standard_sequence(p_);
    } else {
      non_square_free_part_=NT(1);
    }


    typedef typename Traits::Root_bound RBE;

    RBE rbe = traits_.root_bound_object(false);
    
    NT lb, ub;
    if (start_== -Root::infinity()) {
      lb= -rbe(p_);
    } else {
      lb= start_.isolating_interval().first;
      while (sign_at(p_, lb) == CGAL::ZERO) {
	//std::cout << "Having to step off start from " << ub << " to ";
	lb -= .000000001;
	//std::cout << lb << std::endl;
      }
    }
    if (finish_== Root::infinity()) {
      ub= rbe(p_);
    } else {
      ub= finish_.isolating_interval().second;
      while (sign_at(p_, ub) == CGAL::ZERO) {
	ub += .000000001;
      }
    }
    
    CGAL_KINETIC_STURM_DEBUG("initial interval is " << lb << "..." << ub << " which is " 
			     << CGAL::to_double(lb) << " to " << CGAL::to_double(ub));
    
    for (unsigned int i=0; i< sseq_.size(); ++i){
      CGAL_KINETIC_STURM_DEBUG(sseq_[i]);
    }

    root_counter_ = Root_count(sseq_, traits_);
    unsigned int lc= root_counter_(lb);
    unsigned int uc= root_counter_(ub);
    CGAL_KINETIC_STURM_DEBUG("There are " << lc-uc << " roots");

    if (lc != uc) {
      intervals_.push_back(ID(std::make_pair(lb, ub), lc, uc));
    }
  }

public:
  Sturm_root_stack() {}

  //==============
  // CONSTRUCTORS
  //==============
  Sturm_root_stack(const Polynomial& p,
		   const Root& start = -Root::infinity(),
		   const Root& end = Root::infinity(),
		   const Traits& tr = Traits())
    : start_(start), finish_(end), traits_(tr), root_counter_(),
      p_(p), another_even_(false) {
    initialize();
    if (!intervals_.empty()) {
      do_pop(true);
    } else {
      current_= Root();
    }
    start_=Root();
  }
 
  //=============
  // METHODS
  //=============
protected:


  void do_pop(bool first=false) {
    if (intervals_.empty()) {
      CGAL_precondition(current_ != Root());
      current_=Root();
    } else {
      CGAL_precondition(!intervals_.empty());
      do {
	subdivide( );
	std::pair<NT, NT> iv= intervals_.back().interval_;
	intervals_.pop_back();
	CGAL_postcondition(iv.first == iv.second || sign_at(p_, iv.first) != sign_at(p_, iv.second));
	if (iv.first == iv.second) current_= Root(iv.first);
	else current_ = Root(iv, p_, sign_at(p_, iv.first), /*sign_at(p_, iv.second)*/ CGAL::ZERO, traits_);
      } while ( first && current_ < start_);
      CGAL_KINETIC_STURM_DEBUG("Returning " << current_);
      if (current_>= finish_) {
	CGAL_KINETIC_STURM_DEBUG("Out of roots");
	current_= Root();
      } else if (current_.isolating_interval().first == current_.isolating_interval().second ) {
	int mult = traits_.multiplicity_object()(non_square_free_part_, current_.isolating_interval().first);
	CGAL_KINETIC_STURM_DEBUG("Rational multiplicity of " << mult);
	if (mult%2 ==1) another_even_=true; // +1 for the original function
      } else {
	CGAL_postcondition(sign_at(non_square_free_part_, current_.isolating_interval().first)  != CGAL::ZERO);
	CGAL_postcondition(sign_at(non_square_free_part_, current_.isolating_interval().second)  != CGAL::ZERO);
	if (sign_at(non_square_free_part_, current_.isolating_interval().first) 
	    != sign_at(non_square_free_part_, current_.isolating_interval().second)) {
	  CGAL_KINETIC_STURM_DEBUG("Created even root" << current_);
	  another_even_=true;
	}
      }

    }
  }


public:
  const Root& top() const
  {
    return current_;
  }

  bool empty() const
  {
    return current_ == Root();
  }

  void pop()
  {
    if (another_even_) {
      another_even_=false;
    } else {
      do_pop();
    }
  }

  std::ostream& write(std::ostream& os) const
  {
    /*for (unsigned int i = 0; i < sseq.size(); i++) {
      os << sseq[i] << std::endl;
      }*/

    /*typename Interval_container::iterator ivl_it;
      typename std::list<uint_pair>::iterator nr_it;
      for (ivl_it = i_list.begin(), nr_it = s_variations.begin();
      ivl_it != i_list.end(); ++ivl_it, ++nr_it) {
      os << "{";
      ivl_it->write(os);
      int nr = nr_it->first - nr_it->second;
      os << ", " << nr << "} ";
      }
      os << std::endl;*/
    os << sseq_[0] << "->" ;
    os << top();
	  
    return os;
  }

  /*  bool operator==(const This &o) const {
    
  }*/

protected:
  CGAL::Sign sign_at(const Polynomial &p, const NT &nt) const {
    return traits_.sign_at_object()(p, nt);
  }

  Root                             start_, finish_;
  Traits                           traits_;
  Root_count                       root_counter_;
  Polynomial                       p_;
  Polynomial                       non_square_free_part_;
  Standard_sequence                sseq_;
  Root                             current_;
  bool                             another_even_;
  std::vector<ID>                  intervals_;

};

template<class T>
std::ostream& operator<<(std::ostream& os,
			 const Sturm_root_stack<T>& solver)
{
  return solver.write(os);
}


CGAL_POLYNOMIAL_END_NAMESPACE
#endif                                            // CGAL_STURM_LAZY_SOLVER_H
