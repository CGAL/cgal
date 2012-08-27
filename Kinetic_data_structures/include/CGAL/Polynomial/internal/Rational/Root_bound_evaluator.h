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

#ifndef CGAL_POLYNOMIAL_ROOT_BOUND_EVALUATOR_H
#define CGAL_POLYNOMIAL_ROOT_BOUND_EVALUATOR_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/interval_arithmetic.h>
//#include <CGAL/Polynomial/internal/nt_converters.h>
#include <cmath>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template<class Polynomial, class M_t = CGAL::Field_tag>
class Root_bound_evaluator
{
    public:
        typedef typename Polynomial::NT   NT;
        typedef M_t              Method_tag;

        typedef NT               result_type;
        typedef bool             argument_type;
        typedef Method_tag       argument_type1;
        typedef bool             argument_type2;

    private:
        static NT
            compute_bound(const NT& max_abs, const NT& abs_lead_cf,
        CGAL::Field_tag) {
            NT bound = max_abs / abs_lead_cf + NT(1);
//    bound = NT( ceil(CGAL::to_double(bound)) );
            return bound;
        }

        static NT
            compute_bound(const NT& max_abs, const NT& abs_lead_cf,
			  CGAL::Integral_domain_without_division_tag) {
#if 1
            double d1 = to_double(max_abs);
            double d2 = to_double(abs_lead_cf);
            double bound = std::ceil(d1 / d2) + 1.0;
#else
// MK: I MAY WANT TO CHANGE TO THE FOLLOWING CODE OR EVEN DO
//     COMPUTATIONS WITH INTERVALS
            double d1 = POLYNOMIAL_NS::to_interval(max_abs).second;
            double d2 = POLYNOMIAL_NS::to_interval(abs_lead_cf).first;
            double bound = std::ceil(d1 / d2) + 1.0;
#endif
            return NT(bound);
        }

    public:
        Root_bound_evaluator( bool pow= true,
            Method_tag t = Method_tag()): power_of_2(pow), tag(t) {}

        result_type operator()(const Polynomial& p) const
        {
            if ( p.is_zero() ) { return NT(0); }

            int deg = p.degree();

            NT M = CGAL::abs( p[0] );

            for (int i = 1; i < deg; i++) {
                NT abs = CGAL::abs( p[i] );
                if ( abs > M ) { M = abs; }
            }

            NT abs_lead_cf = CGAL::abs( p[deg] );

            NT bound = compute_bound(M, abs_lead_cf, tag);

            if ( power_of_2 ) {
                NT pow2 = NT(2);
                while ( bound >= pow2 ) {
                    pow2 *= NT(2);
                }
                return pow2;
            }

            return bound;
        }
    protected:
        bool power_of_2;
        Method_tag tag;
};

#if 0
template<class Kernel, class M_t = CGAL::Field_tag>
class Filtered_root_bound_evaluator
{
    public:
        typedef typename Kernel::Function Polynomial;
        typedef typename Polynomial::NT   NT;
        typedef M_t              Method_tag;

        typedef NT result_type;
        typedef Method_tag       second_argument_type;
        typedef bool             first_argument_type;

    private:
        template <class NTT>
            static NTT
            compute_bound(const NTT& max_abs, const NTT& abs_lead_cf,
        CGAL::Field_tag) {
            NTT bound = max_abs / abs_lead_cf + NTT(1);
//    bound = NT( ceil(CGAL::to_double(bound)) );
            return bound;
        }

        template <class NTT>
            static NTT
            compute_bound(const NTT& max_abs, const NTT& abs_lead_cf,
        CGAL::Integral_domain_without_division_tag) {
#if 1
            double d1 = CGAL::to_double(max_abs);
            double d2 = CGAL::to_double(abs_lead_cf);
            double bound = ceil(d1 / d2) + 1.0;
#else
// MK: I MAY WANT TO CHANGE TO THE FOLLOWING CODE OR EVEN DO
//     COMPUTATIONS WITH INTERVALS
            double d1 = CGAL::to_interval(max_abs).second;
            double d2 = CGAL::to_interval(abs_lead_cf).first;
            double bound = ceil(d1 / d2) + 1.0;
#endif
            return NTT(bound);
        }

        template <class Poly>
        static typename Poly::NT do_stuff(const Poly &pp, bool power_of_2, Method_tag tag) {
            typedef typename Poly::NT NT;
            int deg = pp.degree();

            NT M = CGAL::abs( pp[0] );

            for (int i = 1; i < deg; i++) {
                NT abs = CGAL::abs( pp[i] );
                if ( abs > M ) { M = abs; }
            }

            NT abs_lead_cf = CGAL::abs( pp[deg] );

            NT bound = compute_bound(M, abs_lead_cf, tag);

            if ( power_of_2 ) {
                NT pow2 = NT(2);
                while ( bound >= pow2 ) {
                    pow2 *= NT(2);
                }
                return pow2;
            } else return bound;

        }

    public:
        Filtered_root_bound_evaluator(bool pow = true,
            Method_tag t = Method_tag()): power_of_2(pow), tag(t)  {}

        result_type operator()(const Polynomial& p) const
        {
            if ( p.degree() ==-1) { return result_type(0); }

            FPU_CW_t backup = FPU_get_and_set_cw(CGAL_FE_UPWARD);
            typename Kernel::Interval_kernel::Function::NT db= do_stuff(p.interval_function(), power_of_2, tag);
            FPU_set_cw(backup);

            double rub= db.sup();
            if (rub == infinity<double>()) {
                return result_type(do_stuff(p.exact_function(), power_of_2, tag));
            }
            else {
                return result_type(rub);
            }
        }
    protected:
        bool power_of_2;
        Method_tag tag;
};
#endif

} } } //namespace CGAL::POLYNOMIAL::internal
#endif                                            // CGAL_POLYNOMIAL_ROOT_BOUND_EVALUATOR_H
