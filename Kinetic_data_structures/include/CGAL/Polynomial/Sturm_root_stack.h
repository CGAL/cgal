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

#ifndef CGAL_STURM_LAZY_SOLVER_H
#define CGAL_STURM_LAZY_SOLVER_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Simple_interval_root.h>
#include <CGAL/Polynomial/internal/Sturm_isolating_interval.h>


#include <list>

#include <CGAL/Polynomial/Kernel.h>

namespace CGAL { namespace POLYNOMIAL {

//================================================================
//================================================================
// ************** the lazy Sturm solver **************************
//================================================================
//================================================================

#define CGAL_KINETIC_STURM_DEBUG(x)
// std::cout << x << std::endl;
#define CGAL_KINETIC_STURM_DEBUG_WRITE(x)
// x << std::endl;

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
  void subdivide() const {
    if (intervals_.empty() ) { return; }

    while ( intervals_.back().lbc_ - intervals_.back().ubc_!=1  ) {
      int num_roots= intervals_.back().lbc_ - intervals_.back().ubc_;
      CGAL_KINETIC_STURM_DEBUG("Next interval has " << intervals_.back().lbc_ - intervals_.back().ubc_ << " roots");
      CGAL_precondition( intervals_.back().lbc_ - intervals_.back().ubc_ >1);
      ID ivl = intervals_.back();
      intervals_.pop_back();
      NT mp= (ivl.interval_.first+ ivl.interval_.second)/NT(2.0);
      CGAL_precondition(sign_at(p_, ivl.interval_.first) != CGAL::ZERO);
      CGAL_precondition(sign_at(p_, ivl.interval_.second) != CGAL::ZERO);
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
	  num_roots -= (muc- ivl.ubc_);
	}
	CGAL_assertion(mlc-muc == 1);
	intervals_.push_back(ID(std::make_pair(mp, mp), mlc, muc));
	--num_roots;
	CGAL_KINETIC_STURM_DEBUG("Adding interval of ");
	CGAL_KINETIC_STURM_DEBUG_WRITE(intervals_.back().write(std::cout)); 
	CGAL_assertion(ivl.lbc_ >= mlc);
	if (ivl.lbc_ - mlc >0) {
	  intervals_.push_back(ID(std::make_pair(ivl.interval_.first, mp-offset),ivl.lbc_, mlc));
	  CGAL_KINETIC_STURM_DEBUG("Adding interval of ");
	  CGAL_KINETIC_STURM_DEBUG_WRITE(intervals_.back().write(std::cout));
	  num_roots -= (ivl.lbc_ - mlc ); 
	}
	
	if (num_roots != 0) {
	  std::cerr << "Initial interval is " << ivl.interval_.first  << ", " << ivl.interval_.second << std::endl;
	  std::cerr << "Midpoint interval is is " << mp-offset << ", " << mp+offset << std::endl;
	  std::cerr << "Counts are " << ivl.lbc_ << " " << mlc << " " << muc << " " << ivl.ubc_ << std::endl;
	}
	CGAL_postcondition(!intervals_.empty() && num_roots==0);
      } else {
	unsigned int mpc= root_counter_(mp);
	CGAL_precondition( ivl.ubc_ <=mpc);
	if (mpc-ivl.ubc_!= 0) {
	  intervals_.push_back(ID(std::make_pair(mp, ivl.interval_.second), mpc, ivl.ubc_));
	  CGAL_KINETIC_STURM_DEBUG(std::cout << "Adding interval of ");
	  CGAL_KINETIC_STURM_DEBUG_WRITE(intervals_.back().write(std::cout)); 
	  num_roots -= (mpc-ivl.ubc_);
	}
	CGAL_precondition(mpc <= ivl.lbc_);
	if ( ivl.lbc_-mpc != 0) {
	  intervals_.push_back(ID(std::make_pair( ivl.interval_.first, mp), ivl.lbc_, mpc));
	  CGAL_KINETIC_STURM_DEBUG(std::cout << "Adding interval of ");
	  CGAL_KINETIC_STURM_DEBUG_WRITE(intervals_.back().write(std::cout)); 
	  num_roots -= ivl.lbc_-mpc;
	}
	if (num_roots != 0) {
	  std::cerr << "Initial interval is " << ivl.interval_.first  << ", " << ivl.interval_.second << std::endl;
	  std::cerr << "Midpoint is " << mp << std::endl;
	  std::cerr << "Counts are " << ivl.lbc_ << " " << mpc << ivl.ubc_ << std::endl;
	}
	CGAL_postcondition(!intervals_.empty() && num_roots==0);
      }
    }                                     // end-while
  }                                         // end subdivide method

  void initialize(const Root &start) {
    if ( p_.is_zero() )  { 
      CGAL_KINETIC_STURM_DEBUG("Zero polynomial ");
      return; 
    }
    if (p_.degree() == 1) {
      NT v= -p_[0]/p_[1];
      CGAL_assertion(p_(v) == 0);
      CGAL_KINETIC_STURM_DEBUG("Linear solution of " << v);
      intervals_.push_back(ID(std::make_pair(v,v), 1,0));
      non_square_free_part_= Polynomial(1);
      return;
    }
    if (start==finish_) {
      CGAL_KINETIC_STURM_DEBUG("Empty interval ");
      return;
    }

    Root ninf= -Root::infinity();

    Standard_sequence sseq = Standard_sequence(p_);
    CGAL_postcondition(sseq[sseq.size()-1].degree() > -1);
    if (sseq[sseq.size()-1].degree() > 0) {
      CGAL_KINETIC_STURM_DEBUG("Non-square free: " << sseq[sseq.size()-1]);
      non_square_free_part_= sseq[sseq.size()-1];
      typename Traits::Quotient quo= traits_.quotient_object();
      p_= quo(p_,non_square_free_part_);
      sseq = Standard_sequence(p_);
    } else {
      non_square_free_part_=NT(1);
    }


    typedef typename Traits::Root_bound RBE;

    RBE rbe = traits_.root_bound_object(false);
    
    NT lb, ub;
    if (start== -Root::infinity()) {
      lb= -rbe(p_);
    } else {
      lb= start.isolating_interval().first;
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
    
    for (unsigned int i=0; i< sseq.size(); ++i){
      CGAL_KINETIC_STURM_DEBUG(sseq[i]);
    }

    root_counter_ = Root_count(sseq, traits_);
    unsigned int lc= root_counter_(lb);
    unsigned int uc= root_counter_(ub);
    CGAL_KINETIC_STURM_DEBUG("There are " << lc-uc << " roots");

    if (lc != uc) {
      intervals_.push_back(ID(std::make_pair(lb, ub), lc, uc));
    }
  }

public:
  Sturm_root_stack(): done_(true) {}

  //==============
  // CONSTRUCTORS
  //==============
  Sturm_root_stack(const Polynomial& p,
		   const Root& start = -Root::infinity(),
		   const Root& end = Root::infinity(),
		   const Traits& tr = Traits())
    : finish_(end), traits_(tr), root_counter_(),
      p_(p), current_ok_(false), another_even_(false) {
    CGAL_precondition(start <= end);
    initialize(start);
    done_=false;
    do_pop();
    while (!empty() && current_ < start) {
      CGAL_KINETIC_STURM_DEBUG("Dropping root " << current_ << " because of " << start);
      current_ok_=false;
      do_pop();
    }
  }
 
  //=============
  // METHODS
  //=============
protected:

  void do_pop() const {
    CGAL_precondition(!done_);
    CGAL_precondition(!current_ok_);
    if (intervals_.empty()) {
      CGAL_KINETIC_STURM_DEBUG("No more roots" << std::endl);
      clean();
    } else {
      CGAL_precondition(!intervals_.empty());
      
      subdivide( );
      std::pair<NT, NT> iv= intervals_.back().interval_;
      intervals_.pop_back();
      CGAL_postcondition(iv.first == iv.second || sign_at(p_, iv.first) != sign_at(p_, iv.second));
	current_ok_=true;
      if (iv.first == iv.second) {
      current_= Root(iv.first);
      } else {
	current_ = Root(iv, p_, sign_at(p_, iv.first), /*sign_at(p_, iv.second)*/ CGAL::ZERO, traits_);
      }
      
      CGAL_KINETIC_STURM_DEBUG("Created root " << current_);
      if (current_>= finish_) {
	CGAL_KINETIC_STURM_DEBUG("Root past end: " << current_);
	clean();
      } else if (current_.isolating_interval().first == current_.isolating_interval().second ) {
	int mult = 1+traits_.multiplicity_object()(non_square_free_part_, current_.isolating_interval().first);
	CGAL_KINETIC_STURM_DEBUG("Rational multiplicity of " << mult);
	if (mult%2 ==0) {
	  CGAL_KINETIC_STURM_DEBUG("Created even rational root" << current_);
	  another_even_=true; 
	} 
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

  void clean() const {
    current_=Root();
    //finish_=Root();
    //p_=Polynomial();
    //root_counter_=Root_count();
    done_=true;
  }


public:
  const Root& top() const
  {
    CGAL_precondition(!empty());
    return current_;
  }

  bool empty() const
  {
    if (!done_ && !current_ok_) {
      do_pop();
      CGAL_postcondition(done_ || current_ok_);
    }
    return done_;
  }

  void pop()
  {
    CGAL_precondition(!empty());
    if (another_even_) {
      another_even_=false;
    } else {
      current_ok_=false;
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
    os << p_ << "(" << non_square_free_part_ << ")" << "->" ;
    os << top();
	  
    return os;
  }

  /*  bool operator==(const This &o) const {
    
  }*/

protected:
  CGAL::Sign sign_at(const Polynomial &p, const NT &nt) const {
    return traits_.sign_at_object()(p, nt);
  }

  Root                             finish_;
  Traits                           traits_;
  Root_count                       root_counter_;
  Polynomial                       p_;
  Polynomial                       non_square_free_part_;
  mutable Root                             current_;
  mutable bool                             current_ok_;
  mutable bool                             another_even_;
  mutable bool                             done_;
  mutable std::vector<ID>                  intervals_;

};

template<class T>
std::ostream& operator<<(std::ostream& os,
			 const Sturm_root_stack<T>& solver)
{
  return solver.write(os);
}


} } //namespace CGAL::POLYNOMIAL
#endif                                            // CGAL_STURM_LAZY_SOLVER_H
