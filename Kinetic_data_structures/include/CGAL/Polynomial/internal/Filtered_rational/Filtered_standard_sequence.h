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

#ifndef CGAL_FILTERED_STANDARD_SEQUENCE_H
#define CGAL_FILTERED_STANDARD_SEQUENCE_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Filtered_rational/Filtered_Sturm_sequence.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template<class Traits>
class Filtered_standard_sequence
: public Filtered_Sturm_sequence<Traits>
{
    protected:
        typedef Filtered_Sturm_sequence<Traits> P;
        typedef typename P::SP                SP;
        typedef typename P::EP                EP;
        typedef typename P::FH   FH;

    public:
        typedef SP                  Storage_function;
        typedef EP                  Exact_function;
        typedef FH                  Function_handle;
//typedef typename P::Method_tag                 Method_tag;

    protected:
        typedef typename P::Interval_nt  Interval_nt;
        typedef typename P::Exact_nt     Exact_nt;

        typedef typename P::Interval_function  Interval_function;

        typedef typename P::ISturm  ISturm;
        typedef typename P::ESturm     ESturm;

        typedef typename P::Exact_to_interval_function_converter
            Exact_to_interval_function_converter;

    protected:

        virtual void compute_exact() const
        {
            if ( this->know_exact ) { return; }

            Exact_function eq = P::tr_.exact_traits_object().differentiate_object()(P::fhp.exact_function());
            this->eseq = P::tr_.exact_traits_object().Sturm_sequence_object(P::fhp.exact_function(), eq);
            P::update_interval_Sturm_sequence();
            this->know_exact = true;
        }

        void initialize() {
            FPU_CW_t backup = FPU_get_and_set_cw(CGAL_FE_UPWARD);

            Interval_function iq = P::tr_.interval_traits_object().differentiate_object()(P::fhp.interval_function());
            this->iseq = P::tr_.interval_traits_object().Sturm_sequence_object(P::fhp.interval_function(), iq);

            FPU_set_cw(backup);

            this->know_exact = false;

#if 1
            bool need_exact = false;
            for (int i = P::iseq.size()-1; i >= 0; i--) {
                bool is_zero = true;
                for (int j = 0; j <= P::iseq[i].degree(); j++) {
                    if ( !CGAL::is_finite(P::iseq[i][j]) ) {
                        need_exact = true;
                        break;
                    }
                    if ( P::iseq[i][j].inf() != 0 || P::iseq[i][j].sup() != 0 ) {
                        is_zero = false;
                    }
                }
                if ( is_zero ) {
                    need_exact = true;
                    break;
                }
            }

            if ( !need_exact ) {
                P::know_exact = false;
                return;
            }

            compute_exact();
#endif
        }

    public:
//===============
// CONSTRUCTORS
//===============
        Filtered_standard_sequence() {}

        Filtered_standard_sequence(const Function_handle& fh, const Traits &tr):
        P(fh, Function_handle(0), tr, false) {
//this->fhp = fh;
            initialize();
        }

        virtual ~Filtered_standard_sequence() {}

        const Interval_function& interval(int i) const  { return P::iseq[i]; }
        const Exact_function&    exact(int i)    const
        {
            if ( !this->know_exact ) {
                compute_exact();
            }
            return P::eseq[i];
        }

    public:

        template<class T>
            unsigned int
            number_of_real_roots(const T& a, const T& b) const
        {
            CGAL_precondition( b >= a );

            unsigned int Va = sign_variations(a);
            if ( Va == 0 ) { return 0; }

            unsigned int Vb = sign_variations(b);

            CGAL_assertion( Va > Vb );

            return Va - Vb;
        }

};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif                                            // CGAL_FILTERED_STANDARD_SEQUENCE_H
