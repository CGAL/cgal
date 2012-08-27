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

#ifndef CGAL_POLYNOMIAL_FILTERED_POLYNOMIAL_RATIONAL_KERNEL_H
#define CGAL_POLYNOMIAL_FILTERED_POLYNOMIAL_RATIONAL_KERNEL_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/Tools/interval_arithmetic.h>

#include <CGAL/Polynomial/Tools/Filtered_function.h>
#include <CGAL/Polynomial/internal/Virtual_function_explicit.h>
#include <CGAL/Polynomial/internal/Virtual_function_constant.h>
#include <CGAL/Polynomial/internal/Virtual_function_generator.h>

#include <CGAL/Polynomial/Filtered_kernel/Filtered_sign_at.h>
#include <CGAL/Polynomial/Filtered_kernel/Filtered_are_negations.h>
#include <CGAL/Polynomial/Filtered_kernel/Filtered_Descartes_has_root.h>
#include <CGAL/Polynomial/Filtered_kernel/Filtered_root_bound_evaluator.h>
#include <CGAL/Polynomial/Filtered_kernel/Filtered_sign_at.h>
#include <CGAL/Polynomial/Filtered_kernel/Filtered_root_multiplicity.h>

// not yet ported
#include <CGAL/Polynomial/Kernel/Bezier_root_counter.h>
#include <CGAL/Polynomial/Kernel/Sturm_root_counter.h>
#include <CGAL/Polynomial/Kernel/Compare_isolated_roots_in_interval.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template <class SFK_t, class IFK_t, class EFK_t, class IFC_t, class EFC_t, class EIFC_t>
class Filtered_polynomial_rational_kernel
{
    typedef Filtered_polynomial_rational_kernel<SFK_t, IFK_t, EFK_t,
        IFC_t, EFC_t, EIFC_t> This;
    public:

        Filtered_polynomial_rational_kernel(){}

        typedef typename SFK_t::Function::NT NT;
        typedef Filtered_function<typename SFK_t::Function,
            typename IFK_t::Function,
            typename EFK_t::Function,
            IFC_t, EFC_t, EIFC_t> Function;

        typedef typename SFK_t::Function Explicit_function;

        typedef SFK_t Explicit_kernel;
        typedef IFK_t Interval_kernel;
        typedef EFK_t Exact_kernel;

        const Explicit_kernel &explicit_kernel_object() const { return sk_;}
        const Interval_kernel &interval_kernel_object() const {return ik_;}
        const Exact_kernel &exact_kernel_object() const {return ek_;}

        typedef EFC_t Exact_function_converter;
        typedef IFC_t Interval_function_converter;
        typedef EIFC_t Exact_interval_function_converter;

        const Exact_function_converter &exact_function_converter_object() const
        {
            return efc_;
        }

        const Interval_function_converter &interval_function_converter_object() const
        {
            return ifc_;
        }

        const Exact_interval_function_converter &exact_interval_function_converter_object() const
        {
            return eifc_;
        }

        typedef internal::Filtered_sign_at<This> Sign_at;
        Sign_at sign_at_object(const Function &p) const
        {
            return Sign_at(p, *this);
        }

        typedef internal::Filtered_Descartes_root_counter<This> Descartes_root_counter;
        Descartes_root_counter Descartes_root_counter_object(const Function &p) const
        {
            return Descartes_root_counter(p, *this);
        }

        typedef internal::Filtered_are_negations<This> Are_negations;
        Are_negations are_negations_object() const
        {
            return Are_negations(*this);
        }

        typedef internal::Filtered_Descartes_has_root<This> Descartes_has_root;
        Descartes_has_root Descartes_has_root_object(const Function &p) const
        {
            return Descartes_has_root(p, *this);
        }

        typedef internal::Sturm_root_counter<Exact_kernel> Sturm_root_counter;
        Sturm_root_counter Sturm_root_counter_object(const Function &p) const
        {
            return Sturm_root_counter(p.exact_function());
        }

        typedef internal::Bezier_root_counter<Exact_kernel> Bezier_root_counter;
        Bezier_root_counter Bezier_root_counter_object(const Function &p) const
        {
            return Bezier_root_counter(p.exact_function());
        }

        typedef internal::Compare_isolated_roots_in_interval<Exact_kernel>
            Compare_isolated_roots_in_interval;
        Compare_isolated_roots_in_interval compare_isolated_roots_in_interval_object(const Function &p0,
            const Function &p1)const
        {
            return Compare_isolated_roots_in_interval(p0.exact_function(), p1.exact_function(),
                exact_kernel_object());
        }

        typedef internal::Filtered_root_bound_evaluator<This> Root_bound;
        Root_bound root_bound_object(bool power_of_two=true)const
        {
            return Root_bound(power_of_two, *this);
        }

        typedef internal::Filtered_root_multiplicity<This> Zero_multiplicity;
        Zero_multiplicity zero_multiplicity_object(const Function &p0) const
        {
            return Zero_multiplicity(p0, *this);
        }

    protected:
        typedef internal::Virtual_function_explicit<
            typename Explicit_kernel::Function,
            typename Interval_kernel::Function,
            typename Exact_kernel::Function,
            Interval_function_converter,
            Exact_function_converter,
            Exact_interval_function_converter> EVF;
    public:

        template <class UF>
            Function function_from_generator(const UF &fc) const
        {
            typename Function::VFP vfp= new internal::Virtual_function_generator<UF, This,
                typename Explicit_kernel::Function,
                typename Interval_kernel::Function,
                typename Exact_kernel::Function,
                Interval_function_converter,
                Exact_function_converter,
                Exact_interval_function_converter>(fc, *this);
            return Function(vfp);
        }

//! construct high degree polynomials
        Function function_object(const NT& a0, const NT& a1=0) const
        {
            typename Explicit_kernel::Function f= sk_.function_object(a0, a1);
            typename Function::VFP vfp= new EVF(f, ifc_, efc_,eifc_);
            return Function(vfp);
        }

//! construct high degree polynomials
        Function function_object(const NT& a0, const NT& a1,
            const NT& a2, const NT& a3=0) const
        {
            typename Explicit_kernel::Function f= sk_.function_object(a0, a1,a2,a3);
            typename Function::VFP vfp= new EVF(f, ifc_, efc_, eifc_);
            return Function(vfp);
        }

//! construct high degree polynomials
        Function function_object(const NT& a0, const NT& a1,
            const NT& a2, const NT& a3,
            const NT& a4, const NT& a5=0,
            const NT& a6=0, const NT& a7=0) const
        {
            typename Explicit_kernel::Function f= sk_.function_object(a0, a1,a2,a3,a4,a5,a6,a7);
            typename Function::VFP vfp= new EVF(f, ifc_, efc_, eifc_);
            return Function(vfp);
        }

//! construct high degree polynomials
        Function function_object(const NT& a0, const NT& a1,
            const NT& a2, const NT& a3,
            const NT& a4, const NT& a5,
            const NT& a6, const NT& a7,
            const NT& a8, const NT& a9=0,
            const NT& a10=0, const NT& a11=0,
            const NT& a12=0, const NT& a13=0,
            const NT& a14=0, const NT& a15=0,
            const NT& a16=0, const NT& a17=0,
            const NT& a18=0, const NT& a19=0) const
        {
            typename Explicit_kernel::Function f= sk_.function_object(a0, a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,
                a11,a12,a13,a14,a15,a16,a17,a18,a19);
            typename Function::VFP vfp= new EVF(f, ifc_, efc_, eifc_);
            return Function(vfp);
        }

    protected:
        Explicit_kernel sk_;
        Interval_kernel ik_;
        Exact_kernel ek_;
        Interval_function_converter ifc_;
        Exact_function_converter efc_;
        Exact_interval_function_converter eifc_;
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
