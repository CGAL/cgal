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

#ifndef CGAL_POLYNOMIAL_INTERNAL_SUBRESULTANT_STURM_SEQUENCE_H
#define CGAL_POLYNOMIAL_INTERNAL_SUBRESULTANT_STURM_SEQUENCE_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Sign_variations_counter.h>
#include <CGAL/Polynomial/internal/Rational/Sturm_sequence_base.h>
/*!
  \file Monic_Sturm_sequence.h A non-filtered Sturm sequence class
  that consists monic polynomials.
*/

namespace CGAL { namespace POLYNOMIAL { namespace internal {
template<class Kernel_t>
class Subresultant_Sturm_sequence
: public Sturm_sequence_base<Kernel_t>
{
    protected:
        typedef Sturm_sequence_base<Kernel_t>    Base;

    public:
        typedef typename Base::Kernel            Kernel;

        typedef typename Base::Polynomial        Polynomial;
        typedef typename Base::NT                NT;

    protected:

        static NT abs(const NT& x) {
            return CGAL::abs(x);
        }

        static NT power(const NT& x, int n) {
            if ( n < 0 ) { return NT(1) / power(x, -n); }

            NT y = x;
            NT z(1);

//  z = x;
            int m = n;
            while (n > 1) {
                m = m / 2;
                if ( n > 2 * m ) {
                    z = z * y;
                }
                y = y * y;
                n = m;
            }

            z = z * y;
            return z;
        }

        void compute_sequence(const Polynomial& p, const Polynomial& q) {
// I HAVE TO FIX THE BUG THAT EXISTS HERE; THE FOLLOWING CODE MAY
// NOT WORK CORRECTLY IF p IS THE ZERO POLYNOMIAL AND q IS NOT
// IN GENERAL I HAVE TO CONSIDER ALL THE LIMITING CASES
            if ( p.degree() >= 0 ) {
                add( p );
                this->size_++;
            }

            if ( q.degree() == -1 ) { return; }

            this->add( q );
            this->size_++;

            if ( p.degree() < q.degree() ) {
                Polynomial r = -this->seq_[0];
                this->add( r );
                this->size_++;
            }

            NT a, b, psi;
            NT psi_old(-1);
            Polynomial r;

            bool is_first = true;
            while ( true ) {
                int deg2 = this->seq_[this->size_ - 2].degree();
                int deg1 = this->seq_[this->size_ - 1].degree();
                int degdiff21 = deg2 - deg1;
                a = power(this->seq_[this->size_ - 1][deg1], degdiff21 + 1);

                r = this->k_.remainder_object()(a * this->seq_[this->size_ - 2],
                    this->seq_[this->size_ - 1]);
                if ( r.is_zero() ) { break; }

                if ( is_first ) {
                    is_first = false;

                    psi = NT(-1);
                    b = NT(1);
                }
                else {
                    int deg3 = this->seq_[this->size_ - 3].degree();
                    int degdiff32 = deg3 - deg2;
                    NT lcoeff2 = this->seq_[this->size_ - 2][deg2];
                    psi = power(-lcoeff2, degdiff32) * power(psi_old, 1 - degdiff32);
                    b = -lcoeff2 * power(psi, degdiff21);
                }
                if ( a * b > 0 ) { b = -b; }

// THE FOLLOWING HACK HAS BEEN DONE SO THAT MP_Float HOPEFULLY
// DOES NOT RUN OUT OF EXPONENT BITS WHEN THE STURM SEQUENCE IS
// COMPUTED
                this->normalize(r, NT());

                this->add( r / b );
                psi_old = psi;
                this->size_++;
            }
        }

    public:
        Subresultant_Sturm_sequence() : Base() {}

        Subresultant_Sturm_sequence(const Polynomial& p, const Polynomial& q,
            const Kernel &k)
        : Base(p, q, k) {
            compute_sequence(p, q);
        }
};

} } } //namespace CGAL::POLYNOMIAL::internal

#endif // CGAL_POLYNOMIAL_INTERNAL_SUBRESULTANT_STURM_SEQUENCE_H
