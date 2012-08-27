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

#ifndef CGAL_POLYNOMIAL_INTERNAL_PRIMITIVE_PART_STURM_SEQUENCE_H
#define CGAL_POLYNOMIAL_INTERNAL_PRIMITIVE_PART_STURM_SEQUENCE_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/internal/Sign_variations_counter.h>
#include <CGAL/Polynomial/internal/Rational/Sturm_sequence_base.h>
/*!
  \file Monic_Sturm_sequence.h A non-filtered Sturm sequence class
  that consists monic polynomials.
*/

namespace CGAL { namespace POLYNOMIAL { namespace internal {
template<class Kernel_t>
class Primitive_part_Sturm_sequence
: public Sturm_sequence_base<Kernel_t>
{
    protected:
        typedef Sturm_sequence_base<Kernel_t>    Base;

    public:
        typedef typename Base::Kernel            Kernel;

        typedef typename Base::Polynomial        Polynomial;
        typedef typename Base::NT                NT;

    protected:

        template<class RT>
        static RT abs(const RT& x) {
            return CGAL::abs(x);
        }

        template<class RT>
        static RT gcd(const RT& x, const RT& y) {
            return CGAL::gcd(x, y);
        }

        template<class RT>
        static RT lcm(const RT& x, const RT& y) {
            return (x * y) / CGAL::gcd(x, y);
        }

        template<class RT>
        static RT compute_gcd(const CGAL_POLYNOMIAL_NS::Polynomial<RT>& p) {
            int deg = p.degree();
            if ( deg < 0 ) { return RT(1); }
            if ( deg == 0 ) { return abs(p[0]); }

            RT gcd_ = gcd(abs(p[0]), abs(p[1]));
            for (int i = 2; i <= deg; i++) {
                gcd_ = gcd(gcd_, abs(p[i]));
            }
            return gcd_;
        }

        template<class RT>
        static RT compute_lcm(const CGAL_POLYNOMIAL_NS::Polynomial<RT>& p) {
            int deg = p.degree();
            CGAL_assertion( deg >= 1 );

            RT lcm_ = lcm(abs(p[0]), abs(p[1]));
            for (int i = 2; i <= deg; i++) {
                lcm_ = lcm(lcm_, abs(p[i]));
            }
            return lcm_;
        }

        template<class RT>
        static RT compute_lazy_lcm(const CGAL_POLYNOMIAL_NS::Polynomial<RT>& p) {
            int deg = p.degree();
            CGAL_assertion( deg >= 1 );

            RT lcm_ = abs(p[0] * p[1]);
            for (int i = 2; i <= deg; i++) {
                lcm_ *= abs(p[i]);
            }
            return lcm_;
        }

        static Polynomial compute_integer_polynomial(const Polynomial& p) {
            typedef typename Rational_traits<NT>::RT    RT;
            typedef CGAL_POLYNOMIAL_NS::Polynomial<RT>       RT_Polynomial;
            int deg = p.degree();
            if ( deg < 0 ) { return p; }
            if ( deg == 0 ) {
                if ( p[0] > 0 ) {
                    return Polynomial(NT(1));
                }
                else {
                    return Polynomial(NT(-1));
                }
            }

            Rational_traits<NT> rational_traits;
            std::vector<RT> denominators(deg+1);
            for (int i = 0; i <= deg; i++) {
                denominators[i] = rational_traits.denominator(p[i]);
            }
            RT lcm_ = compute_lcm( RT_Polynomial(denominators.begin(),
                denominators.end()) );
            NT lcmq_ = rational_traits.make_rational(lcm_, RT(1));
            return p * lcmq_;
        }

        static Polynomial compute_primitive_polynomial(const Polynomial& p) {
            typedef typename Rational_traits<NT>::RT    RT;
            typedef CGAL_POLYNOMIAL_NS::Polynomial<RT>       RT_Polynomial;
// we assume that p has only integer coefficients
            int deg = p.degree();
            if ( deg < 0 ) { return p; }
            if ( deg == 0 ) {
                if ( p[0] > 0 ) {
                    return Polynomial(NT(1));
                }
                else {
                    return Polynomial(NT(-1));
                }
            }

            Rational_traits<NT> rational_traits;
            std::vector<RT> numerators(deg+1);
            for (int i = 0; i <= deg; i++) {
                numerators[i] = rational_traits.numerator(p[i]);
            }
            RT gcd_ = compute_gcd( RT_Polynomial(numerators.begin(),
                numerators.end()) );
            NT gcdq_ = rational_traits.make_rational(gcd_, RT(1));
            return p / gcdq_;
        }

        void compute_sequence(const Polynomial& p, const Polynomial& q) {
// I HAVE TO FIX THE BUG THAT EXISTS HERE; THE FOLLOWING CODE MAY
// NOT WORK CORRECTLY IF p IS THE ZERO POLYNOMIAL AND q IS NOT
// IN GENERAL I HAVE TO CONSIDER ALL THE LIMITING CASES
            Polynomial ip = compute_integer_polynomial(p);
            Polynomial iq = compute_integer_polynomial(q);
            Polynomial p_prim = compute_primitive_polynomial(ip);
            Polynomial q_prim = compute_primitive_polynomial(iq);

            if ( p.degree() >= 0 ) {
                add( p_prim );
                this->size_++;
            }

            if ( q.degree() == -1 ) { return; }

            this->add( q_prim );
            this->size_++;

            if ( p.degree() < q.degree() ) {
                this->add( -this->seq_[0] );
                this->size_++;
            }

            Polynomial r;

            while ( true ) {
                r = -this->k_.remainder_object()(this->seq_[this->size_ - 2],
                    this->seq_[this->size_ - 1]);
                if ( r.is_zero() ) { break; }

// THE FOLLOWING HACK HAS BEEN DONE SO THAT MP_Float HOPEFULLY
// DOES NOT RUN OUT OF EXPONENT BITS WHEN THE STURM SEQUENCE IS
// COMPUTED
                this->normalize(r, NT());

                Polynomial ir = compute_integer_polynomial(r);
                this->add( compute_primitive_polynomial(ir) );
                this->size_++;
            }
        }

    public:
        Primitive_part_Sturm_sequence() : Base() {}

        Primitive_part_Sturm_sequence(const Polynomial& p, const Polynomial& q,
            const Kernel &k)
        : Base(p, q, k) {
            compute_sequence(p, q);
        }
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif                                            // CGAL_POLYNOMIAL_INTERNAL_PRIMITIVE_PART_STURM_SEQUENCE_H
