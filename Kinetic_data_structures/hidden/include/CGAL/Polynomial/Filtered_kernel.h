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

#ifndef CGAL_POLYNOMIAL_FILTERED_POLYNOMIAL_KERNEL_H
#define CGAL_POLYNOMIAL_FILTERED_POLYNOMIAL_KERNEL_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Filtered_rational/Filtered_rational_traits.h>
#include <CGAL/Polynomial/internal/Filtered_kernel/Filtered_sign_at.h>
#include <CGAL/Polynomial/internal/Filtered_kernel/Filtered_root_multiplicity.h>
#include <CGAL/Polynomial/internal/Kernel/Root_container.h>
#include <CGAL/Polynomial/internal/Kernel/Is_even_multiplicity.h>
#include <CGAL/Polynomial/internal/Kernel/Is_rational.h>
#include <CGAL/Polynomial/internal/Kernel/To_rational.h>

namespace CGAL { namespace POLYNOMIAL {

//! A filtered polynomial kernel.
/*!
 */
template < class Traits_t, class RE>
class Filtered_kernel: public internal::Filtered_rational_traits<Traits_t>
{
    typedef Filtered_kernel< Traits_t, RE > This;
    typedef internal::Filtered_rational_traits<Traits_t> P;
    public:
        typedef RE Root_stack;
        typedef typename RE::Root Root;
        typedef typename P::Function Function;
        typedef typename Function::NT NT;
        typedef typename Root_stack::Traits Root_stack_traits;

        Filtered_kernel(const Root_stack_traits &tr = Root_stack_traits()): ret_(tr){}

        typedef typename internal::Filtered_sign_at<This> Sign_at;
        Sign_at sign_at_object(const Function &p) const
        {
            return Sign_at(p, *this);
        }

        typedef internal::Filtered_root_multiplicity<This> Multiplicity;
        Multiplicity multiplicity_object(const Function &p0) const
        {
            return Multiplicity(p0, *this);
        }

        typedef internal::Is_even_multiplicity<This> Is_even_multiplicity;
        Is_even_multiplicity is_even_multiplicity_object(const Function &) const
        {
            return Is_even_multiplicity();
        }

        typedef internal::Is_rational<This> Is_rational;
        Is_rational is_rational_object() const
        {
            return Is_rational();
        }

        typedef internal::To_rational<This> To_rational;
        To_rational to_rational_object() const
        {
            return To_rational();
        }

//! Return a container for roots in an interval
/*!
  \todo make sure that the iterator has all the right types.
*/
        typedef internal::Root_container<This> Root_container;
        friend class internal::Root_container<This>;
        Root_container root_container_object(const typename P::Function &f,
            const Root &lb=-Root::infinity(),
            const Root &ub= Root::infinity()) const
        {
            return Root_container(f, lb, ub, *this);
        }

        Root_stack root_stack_object(const typename P::Function &f,
            const Root &lb=-Root::infinity(),
            const Root &ub= Root::infinity()) const
        {
            return Root_stack(f, lb, ub, root_stack_traits_object());
        }

        Root_stack_traits root_stack_traits_object() const
        {
            return ret_;
        }
    protected:
        Root_stack_traits ret_;
};

} } //namespace CGAL::POLYNOMIAL
#endif
