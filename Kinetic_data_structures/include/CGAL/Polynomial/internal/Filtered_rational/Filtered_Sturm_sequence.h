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

#ifndef CGAL_FILTERED_STURM_SEQUENCE_H
#define CGAL_FILTERED_STURM_SEQUENCE_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/NT_converter.h>
#include <CGAL/Polynomial/internal/Sign_variations_counter.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

#define CGAL_POLYNOMIAL_NORMALIZE_GCD_IF_CONSTANT

template<class Traits>
class Filtered_Sturm_sequence
{
    protected:
        typedef typename Traits::Function::Exact_function                SP;
        typedef typename Traits::Function::Exact_function                EP;
        typedef typename Traits::Function   FH;

    public:
        typedef SP                  Storage_function;
        typedef EP                  Exact_function;
        typedef FH                  Function_handle;
        typedef typename Traits::Method_tag                 Method_tag;

    protected:
        typedef typename FH::Interval_function::NT  Interval_nt;
        typedef typename FH::Exact_function::NT     Exact_nt;

        typedef typename FH::Interval_function  Interval_function;
        typedef Interval_function               IP;

        typedef typename Traits::Exact_to_interval_converter Exact_to_interval_function_converter;

        typedef typename Traits::Interval_traits::Sturm_sequence  ISturm;
        typedef typename Traits::Exact_traits::Sturm_sequence     ESturm;

    protected:

        void update_interval_Sturm_sequence() const
        {
//Exact_to_interval_function_converter e2i;

            iseq.set_size( eseq.size() );

            for (unsigned int i = 0; i < eseq.size(); i++) {
                iseq[i] = e2i( eseq[i] );
            }

#ifdef CGAL_POLYNOMIAL_NORMALIZE_GCD_IF_CONSTANT
            unsigned int igcd = iseq.size() - 1;
            if ( iseq[igcd].degree() == 0 ) {
                Sign s = CGAL::sign( eseq[igcd][0] );
                if ( s == CGAL::POSITIVE ) {
                    iseq[igcd].set_coef(0, Interval_nt(1));
                }
                else {
                    iseq[igcd].set_coef(0, Interval_nt(-1));
                }
            }
#endif
        }

        virtual void compute_exact() const
        {
            if ( know_exact ) { return; }

            eseq = tr_.exact_traits_object().Sturm_sequence_object(fhp.exact_function(), fhq.exact_function());
            update_interval_Sturm_sequence();
            know_exact = true;
        }

        void initialize() {
            FPU_CW_t backup = FPU_get_and_set_cw(CGAL_FE_UPWARD);

            typename Traits::Interval_traits::Sturm_sequence ss= tr_.interval_traits_object().Sturm_sequence_object(fhp.interval_function(), fhq.interval_function());
            iseq = ss;

            FPU_set_cw(backup);

            know_exact = false;
#if 1
            bool need_exact = false;
            for (int i = iseq.size() - 1; i >= 0; i--) {
                bool is_zero = true;
                for (int j = 0; j <= iseq[i].degree(); j++) {
                    if ( !CGAL::is_finite(iseq[i][j]) ) {
                        need_exact = true;
                        break;
                    }
                    if ( iseq[i][j].inf() != 0 || iseq[i][j].sup() != 0 ) {
                        is_zero = false;
                    }
                }
                if ( is_zero ) {
                    need_exact = true;
                    break;
                }
            }

            if ( !need_exact ) {
                know_exact = false;
                return;
            }

            compute_exact();
#endif
        }

    public:
//===============
// CONSTRUCTORS
//===============
        Filtered_Sturm_sequence() {}

        Filtered_Sturm_sequence(const Function_handle& fhp,
            const Function_handle& fhq,
            const Traits &tr, bool init=true)
        : fhp(fhp), fhq(fhq), e2i(tr.exact_to_interval_converter_object()), tr_(tr) {
            if (init )initialize();
        }

        virtual ~Filtered_Sturm_sequence() {}

    public:
        template<class T>
            unsigned int sign_variations(const T& x) const
        {
            CGAL::To_interval<T>  to_interval;

            Interval_nt ix = to_interval(x);

            std::vector<CGAL::Sign>  signs( iseq.size() );

            bool need_exact = false;

            int k = 1;
            while ( k <= 2 ) {
                FPU_CW_t backup = FPU_get_and_set_cw(CGAL_FE_UPWARD);

                for (unsigned int i = 0; i < iseq.size(); i++) {
                    Interval_nt f_ix = iseq[i](ix);

                    if ( ix.inf() == ix.sup() && ix.inf() == 0 ) {
                        f_ix = iseq[i][0];
                    }

                    if ( !CGAL::is_valid(f_ix.inf()) ||
                    !CGAL::is_valid(f_ix.sup()) ) {
                        need_exact = true;
                        FPU_set_cw(backup);
                        break;
                    }

                    if ( f_ix.inf() > 0 ) {
                        signs[i] = CGAL::POSITIVE;
                    }
                    else if ( f_ix.sup() < 0 ) {
                        signs[i] =  CGAL::NEGATIVE;
                    }
                    else if ( f_ix.inf() == f_ix.sup() ) {
                        signs[i] = CGAL::sign(f_ix.inf());
                    }
                    else {
                        need_exact = true;
                        FPU_set_cw(backup);
                        break;
                    }                             // end-if
                }                                 // end-for

                FPU_set_cw(backup);

                if ( !need_exact ) { break; }

                if ( need_exact && k == 1 ) {
                    compute_exact();
                    need_exact = false;
                    signs.resize( eseq.size() );
                }

                k++;
            }

            if ( need_exact ) {
                CGAL::NT_converter<T,Exact_nt>  to_exact;

                Exact_nt ex = to_exact(x);

                for (unsigned int i = 0; i < eseq.size(); i++) {
                    signs[i] = CGAL::sign( eseq[i](ex) );
                }                                 // end-for
            }

            return
                Sign_variations_counter::sign_variations(signs.begin(),
                signs.end());
        }

        unsigned int size() const
        {
            return iseq.size();
        }

        unsigned int exact_size() const
        {
            if ( !know_exact ) { compute_exact(); }
            return eseq.size();
        }

        template<class T>
            CGAL::Sign sign_at(const T& x, unsigned int i) const
        {
            if ( i < iseq.size() ) {
                CGAL::To_interval<T> to_interval;

                Interval_nt ix = to_interval(x);

                FPU_CW_t backup = FPU_get_and_set_cw(CGAL_FE_UPWARD);

                Interval_nt f_ix = iseq[i](ix);

                FPU_set_cw(backup);

                if ( f_ix.sup() < 0 ) { return CGAL::NEGATIVE; }
                if ( f_ix.inf() > 0 ) { return CGAL::POSITIVE; }
                if ( f_ix.inf() == f_ix.sup() ) {
                    return CGAL::sign(f_ix.inf());
                }

                if ( !know_exact ) {
                    compute_exact();

                    FPU_CW_t backup = FPU_get_and_set_cw(CGAL_FE_UPWARD);

                    Interval_nt f_ix = iseq[i](ix);

                    FPU_set_cw(backup);

                    if ( f_ix.sup() < 0 ) { return CGAL::NEGATIVE; }
                    if ( f_ix.inf() > 0 ) { return CGAL::POSITIVE; }
                    if ( f_ix.inf() == f_ix.sup() ) {
                        return CGAL::sign(f_ix.inf());
                    }
                }
            }

            CGAL::NT_converter<T,Exact_nt>  to_exact;

            Exact_nt ex = to_exact(x);

            if ( !know_exact ) { compute_exact(); }

            return CGAL::sign( eseq[i](ex) );
        }

        template<class T>
            CGAL::Sign sign_at_gcd(const T& x) const
        {
            CGAL::To_interval<T> to_interval;

            Interval_nt ix = to_interval(x);

            FPU_CW_t backup = FPU_get_and_set_cw(CGAL_FE_UPWARD);

            Interval_nt f_ix = iseq[iseq.size()-1](ix);

            FPU_set_cw(backup);

            if ( f_ix.sup() < 0 ) { return CGAL::NEGATIVE; }
            if ( f_ix.inf() > 0 ) { return CGAL::POSITIVE; }
            if ( f_ix.inf() == f_ix.sup() ) {
                return CGAL::sign(f_ix.inf());
            }

            if ( !know_exact ) {
                compute_exact();

                FPU_CW_t backup = FPU_get_and_set_cw(CGAL_FE_UPWARD);

                Interval_nt f_ix = iseq[iseq.size()-1](ix);

                FPU_set_cw(backup);

                if ( f_ix.sup() < 0 ) { return CGAL::NEGATIVE; }
                if ( f_ix.inf() > 0 ) { return CGAL::POSITIVE; }
                if ( f_ix.inf() == f_ix.sup() ) {
                    return CGAL::sign(f_ix.inf());
                }
            }

            CGAL::NT_converter<T,Exact_nt>  to_exact;

            Exact_nt ex = to_exact(x);

            if ( !know_exact ) { compute_exact(); }

            return CGAL::sign( eseq[eseq.size()-1](ex) );
        }

//  IP interval(int i) const { return iseq[i]; }
//  EP exact   (int i) const { return eseq[i]; }

    protected:
        Function_handle  fhp, fhq;
        Exact_to_interval_function_converter e2i;
        Traits tr_;
        mutable ISturm  iseq;
        mutable ESturm  eseq;
        mutable bool    know_exact;
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif                                            // CGAL_FILTERED_STURM_SEQUENCE_H
