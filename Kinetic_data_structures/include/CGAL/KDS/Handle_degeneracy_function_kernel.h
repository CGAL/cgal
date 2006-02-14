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
// $Id$ $Date$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_ROOT_DEGEN_FK_H
#define CGAL_KDS_ROOT_DEGEN_FK_H
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_NAMESPACE;

template <class Traits_t>
struct Handle_degeneracy_function_kernel: public Traits_t
{

    class Root_stack
    {
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
            iem_(k.is_even_multiplicity_object(uf)) {
                CGAL_expensive_precondition(solver_.top() > lb);
                CGAL::POLYNOMIAL::Sign sn= k.sign_between_roots_object(lb, solver_.top())(uf);
                if (sn == CGAL::NEGATIVE) {
                    std::cout << "Degeneracy for " << uf << " at " << lb << std::endl;
                    extra_root_=lb;
                    has_extra_=true;
                }
                else {
                    has_extra_=false;
                }
                one_even_=false;
            }

            Root_stack(){}

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
                }
                else if (!one_even_ && iem_(solver_.top())) {
                    one_even_=true;
                }
                else {
                    solver_.pop();
                    one_even_=false;
                }
            }

            bool empty() const
            {
                return !has_extra_ && solver_.empty();
            }
        protected:
            Wrapped_solver solver_;
            Root extra_root_;
            bool one_even_;
            bool has_extra_;
            typename Traits::Is_even_multiplicity iem_;
    };

    Root_stack root_stack_object(const typename Traits_t::Function &f,
        const typename Traits_t::Root &lb,
    const typename Traits_t::Root &ub) const {
        return Root_stack(f, lb, ub, *this);
    }
};

CGAL_KDS_END_NAMESPACE;
#endif
