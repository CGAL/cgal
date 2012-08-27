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

#ifndef CGAL_FILTERED_SIGN_STURM_SEQUENCE_H
#define CGAL_FILTERED_SIGN_STURM_SEQUENCE_H

#include <CGAL/Polynomial/basic.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template<class Traits >
class Filtered_sign_Sturm_sequence
: public Traits::Sturm_sequence
{
    protected:
        typedef typename Traits::Function::Exact_function SP;
        typedef SP                EP;
        typedef typename Traits::Function   FH;

    public:
        typedef SP                  Storage_function;
        typedef EP                  Exact_function;
        typedef FH                  Function_handle;
        typedef typename Traits::Method_tag                 Method_tag;

        typedef int                 result_type;

    protected:
        typedef typename FH::Interval_function::NT  Interval_nt;
        typedef typename FH::Exact_function::NT     Exact_nt;
//typedef typename Traits::Exact_to_interval_function_converter Exact_to_interval_function_converter;

        typedef typename FH::Interval_function  Interval_function;

        typedef typename Traits::Interval_traits::Sturm_sequence  ISturm;
        typedef typename Traits::Exact_traits::Sturm_sequence     ESturm;

        typedef typename Traits::Sturm_sequence  Base;

    protected:

        virtual void compute_exact() const
        {
            if ( this->know_exact ) { return; }

            Exact_function eq = Base::tr_.exact_traits_object().differentiate_object()(Base::fhp.exact_function())
                * Base::fhq.exact_function();
            this->eseq =  Base::tr_.exact_traits_object().Sturm_sequence_object(Base::fhp.exact_function(), eq);
            Base::update_interval_Sturm_sequence();
            this->know_exact = true;
        }

        void initialize() {
            FPU_CW_t backup = FPU_get_and_set_cw(CGAL_FE_UPWARD);

            Interval_function iq = Base::tr_.interval_traits_object().differentiate_object()(Base::fhp.interval_function())
                * Base::fhq.interval_function();
            this->iseq = Base::tr_.interval_traits_object().Sturm_sequence_object(Base::fhp.interval_function(), iq);

            FPU_set_cw(backup);

            this->know_exact = false;
        }

    public:
        Filtered_sign_Sturm_sequence() {}

        Filtered_sign_Sturm_sequence(const Function_handle& fhp,
            const Function_handle& fhq,
        const Traits &tr): Base(fhp, fhq, tr,false) {
            initialize();
        }

    public:
        template<class T>
            unsigned int sign_variations(const T& x) const
        {
            return sign_variations_base(x);
        }

        template<class T>
            result_type operator()(const T& a, const T& b) const
        {
            return sum_of_signs(a, b);
        }

        template<class T>
            int sum_of_signs(const T& a, const T& b) const
        {
            CGAL_precondition( b >= a );

            unsigned int Va = sign_variations_base(a);
            if ( Va == 0 ) { return 0; }

            unsigned int Vb = sign_variations_base(b);

            int diff = static_cast<int>(Va) - static_cast<int>(Vb);

            return diff;
        }

    protected:
//Function_handle fhq_;
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif                                            // CGAL_FILTERED_SIGN_STURM_SEQUENCE_H
