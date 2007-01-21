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

#ifndef CGAL_KINETIC_ROOT_DEGEN_FK_H
#define CGAL_KINETIC_ROOT_DEGEN_FK_H
#include <CGAL/Kinetic/basic.h>

CGAL_KINETIC_BEGIN_NAMESPACE;


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
      CGAL_expensive_precondition(solver_.empty() || solver_.top() > lb);
#ifndef NDEBUG
      std::cout << "Function= " << uf << std::endl;
#endif
      CGAL::Sign sn;
      if (solver_.empty()) {
	sn = k.sign_between_roots_object(lb, ub)(uf);
      } else {
	sn = k.sign_between_roots_object(lb, solver_.top())(uf);
      }
      if (sn == CGAL::NEGATIVE) {
	if (!SLOPPY) {
#ifndef NDEBUG
	  if (k.sign_at_object(uf)(lb) == CGAL::NEGATIVE) {
	    CGAL_KINETIC_ERROR( "Invalid certificate constructed for function " << uf << " between " << lb 
				<< " and " << ub << " will fail immediately." << std::endl);
	    CGAL_exactness_assertion(sn != CGAL::NEGATIVE);
	  }
#endif
	  extra_root_=lb;
	  has_extra_=true;
	} else {
	  if (!solver_.empty()) solver_.pop();
	  has_extra_=false;
	}
      }
      else {
	has_extra_=false;
      }
    }

    HDRS(){}

    //! Drop even roots
    const Root& top() const
    {
      if (has_extra_) return extra_root_;
      else return solver_.top();
    }

    void pop() {
      if (has_extra_) {
	extra_root_=Root();
	has_extra_=false;
      } else {
	solver_.pop();
      }
    }

    bool empty() const
    {
      return !has_extra_ && solver_.empty();
    }
    double estimate() const {
      return solver_.estimate();
    }

    std::ostream &write(std::ostream& out) const {
      out << solver_;
      return out;
    }

  protected:
    Wrapped_solver solver_;
    Root extra_root_;
    bool has_extra_;
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

CGAL_KINETIC_END_NAMESPACE;
#endif
