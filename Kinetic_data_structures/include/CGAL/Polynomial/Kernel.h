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

#ifndef CGAL_POLYNOMIAL_POLYNOMIAL_KERNEL_H
#define CGAL_POLYNOMIAL_POLYNOMIAL_KERNEL_H
#include <CGAL/Polynomial/basic.h>

#include <CGAL/Polynomial/internal/Kernel/Multiplicity.h>
#include <CGAL/Polynomial/internal/Kernel/Rational_between_roots.h>
#include <CGAL/Polynomial/internal/Kernel/Root_container.h>
#include <CGAL/Polynomial/internal/Kernel/Isolating_interval.h>
#include <CGAL/Polynomial/internal/Kernel/Sign_above.h>
#include <CGAL/Polynomial/internal/Kernel/Sign_at.h>
#include <CGAL/Polynomial/internal/Kernel/Sign_below.h>
#include <CGAL/Polynomial/internal/Kernel/Sign_between_roots.h>
#include <CGAL/Polynomial/internal/Kernel/Is_even_multiplicity.h>
#include <CGAL/Polynomial/internal/Kernel/Is_rational.h>
#include <CGAL/Polynomial/internal/Kernel/To_rational.h>
#include <CGAL/Polynomial/internal/Rational/Rational_traits_base.h>
#include <CGAL/Polynomial/internal/Kernel/Lower_bound_root.h>

namespace CGAL { namespace POLYNOMIAL {

//! The polynomial kernel.
/*!  Operations on rationals are handled by
  internal::Polynomial_rational_kernel. This kernel must be kept
  seperate from the rational kernel as the solver must be able to use
  (and store) the rational kernel and it will not compile if there is
  only one kernel. I am not sure why.

  Basically the following does not compile:
  template <class K>
  struct R{
  K k;
  };
  template <class K>
  struct S{
  typedef R<K> Rt;
  Rt r;
  };
  struct K {
  typedef S<K> St;
  typedef St::Rt Rt;
  St s_o(){
  St s;
  return s;
  }
  };

  The other reason is that the Filtered_kernel just needs rational
  kernels for the non-filtered types. At the moment I don't use this
  since I don't feel like exposing the rational kernel.
*/
template <class Polynomial_t, class Root_stack_t, class NT_t= typename Polynomial_t::NT>
class Kernel: public internal::Rational_traits_base<Polynomial_t>
{
  typedef Kernel<Polynomial_t, Root_stack_t, NT_t> This;
  typedef typename internal::Rational_traits_base<Polynomial_t> P;
public:
  typedef Root_stack_t Root_stack;
  typedef typename Root_stack_t::Root Root;
  typedef Polynomial_t Function;
  typedef NT_t FT;
  typedef typename Root_stack_t::Traits Root_stack_traits;

  //! \todo do something with tr
  Kernel(const Root_stack_traits &tr=Root_stack_traits()):
    solver_traits_(tr){}

  typedef internal::Sign_at<Root, This> Sign_at;
  Sign_at sign_at_object() const
  {
    return Sign_at(*this);
  }

  //! Compute the multiplicity of a zero.
  /*!
    The value passed must be a rational number. Is there a better name?

    \todo fix the functor to make it work on roots
  */
  typedef internal::Multiplicity<This> Multiplicity;
  Multiplicity multiplicity_object(const Function &p0) const
  {
    return Multiplicity(p0, *this);
  }

  //! Compute the sign of p immediately after a root of another function (or of p)
  typedef internal::Sign_above<Root, This> Sign_after;
  Sign_after sign_after_object() const
  {
    return Sign_after(*this);
  }

  //! Compute the sign of p immediately after a root of another function (or of p)
  /*typedef internal::Sign_below<Root, This> Sign_below;
  Sign_below sign_be_object() const
  {
    return Sign_below(*this);
    }*/

  //! Find a rational number between two non-equal roots
  typedef internal::Rational_between_roots<This> Rational_between_roots;
  Rational_between_roots rational_between_roots_object() const
  {
    return Rational_between_roots(*this);
  }

  //! Compute the sign between two roots
  typedef internal::Sign_between_roots<This> Sign_between_roots;
  Sign_between_roots sign_between_roots_object() const
  {
    return Sign_between_roots(*this);
  }

  //! Return true if the root has even multiplicity
  /*typedef internal::Is_even_multiplicity<This> Is_even_multiplicity;
  Is_even_multiplicity is_even_multiplicity_object(const Function &) const
  {
    return Is_even_multiplicity();
    }*/

  //! Return true if the root is an exact rational
  /*typedef internal::Is_rational<This> Is_rational;
  Is_rational is_rational_object() const
  {
    return Is_rational();
    }*/

  /*typedef internal::Lower_bound_root<This> Lower_bound_root;
  Lower_bound_root lower_bound_root_object() const {
    return Lower_bound_root();
    }*/

  //! Return the rational value of the root, assuming it is rational
  typedef internal::To_rational<This> To_rational;
  To_rational to_rational_object() const
  {
    return To_rational();
  }

  typedef internal::To_isolating_interval<This> To_isolating_interval;
  To_isolating_interval to_isolating_interval_object() const {
    return To_isolating_interval();
  }

  //! Return a container for roots in an interval
  /*!
    \todo make sure that the iterator has all the right types.
  */
  typedef internal::Root_container<This> Root_container;
  friend class internal::Root_container<This>;
  Root_container root_container_object(const Function &f,
				       const Root &lb,
				       const Root &ub) const
  {
    return Root_container(f, lb, ub, root_stack_traits_object());
  }

  //! Return a root stack
  /*!
    \todo make sure that the iterator has all the right types.
  */
  Root_stack root_stack_object(const Function &f,
			       const Root &lb,
			       const Root &ub) const
  {
    return Root_stack(f, lb, ub, root_stack_traits_object());
  }

  Root_stack_traits root_stack_traits_object() const
  {
    return solver_traits_;
  }

protected:
  Root_stack_traits solver_traits_;

};

} } //namespace CGAL::POLYNOMIAL
#endif
