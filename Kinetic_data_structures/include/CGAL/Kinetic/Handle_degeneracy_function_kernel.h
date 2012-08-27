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

#ifndef CGAL_KINETIC_ROOT_DEGEN_FK_H
#define CGAL_KINETIC_ROOT_DEGEN_FK_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/internal/debug_counters.h>

namespace CGAL { namespace Kinetic {


template <class Traits_t, bool SLOPPY>
class HDRS{
  private:
    typedef typename Traits_t::Root_stack Wrapped_solver;
    typedef typename Traits_t::Function Function;

  public:
    typedef typename Wrapped_solver::Root Root;
    typedef Traits_t Traits;
    //! Construct and check preconditions
    /*!

    */
    HDRS(const Function &uf, const Root& lb,
	 const Root& ub, const Traits_t& k): solver_(k.root_stack_object(uf, lb, ub)) {
      CGAL_LOG(Log::LOTS, "Function= " << uf << std::endl);
      CGAL_expensive_precondition(solver_.empty() || solver_.top() >= lb);
      if (uf.degree() == -1) {
	CGAL_LOG(Log::SOME, "Zero function found at time " << lb << std::endl);	
	++ internal::zero_certificates__;
      }
#ifndef NDEBUG
      if (!SLOPPY && k.sign_at_object()(uf, lb) == CGAL::NEGATIVE) {
	CGAL_ERROR( "Invalid certificate constructed for function " << uf << " between " << lb 
			    << " and " << ub << " will fail immediately." << std::endl);
	CGAL_exactness_precondition(k.sign_at_object()(uf, lb) != CGAL::NEGATIVE);
      }
#endif
      if (solver_.empty()) {
	CGAL_LOG(Log::LOTS, "No failure" << std::endl);
	 //sn = k.sign_between_roots_object()(uf, lb, ub);
      } else if (solver_.top() == lb) {
	CGAL_LOG(Log::LOTS, "Degeneracy at " << solver_.top() << std::endl);
	CGAL::Sign sn = k.sign_after_object()(uf, lb);
	if (sn == CGAL::NEGATIVE) {
	  ++internal::function_degeneracies__;
	  CGAL_LOG(Log::LOTS, "Extra root at lower bound of " << lb << std::endl);
	} else {
	  CGAL_LOG(Log::LOTS, "Popping extra root at lower bound of " << lb << std::endl);
	  solver_.pop();
	}
      } 
    }

    HDRS(){}

    //! Drop even roots
    const Root& top() const
    {
      return solver_.top();
    }

    void pop() {
      solver_.pop();
    }

    bool empty() const
    {
      return solver_.empty();
    }
  /*double estimate() const {
      return solver_.estimate();
      }*/

    std::ostream &write(std::ostream& out) const {
      out << solver_;
      return out;
    }

  protected:
    Wrapped_solver solver_;
  };

template <class Traits_t, bool SLOPPY>
struct Handle_degeneracy_function_kernel: public Traits_t
{

  typedef HDRS<Traits_t, SLOPPY> Root_stack;

  Root_stack root_stack_object(const typename Traits_t::Function &f,
			       const typename Traits_t::Root &lb,
			       const typename Traits_t::Root &ub) const {
    return Root_stack(f, lb, ub, *this);
  }
};

template <class T, bool SLOPPY>
std::ostream &operator<<(std::ostream &out, const HDRS<T, SLOPPY> &k) {
  return k.write(out);
}

} } //namespace CGAL::Kinetic
#endif
