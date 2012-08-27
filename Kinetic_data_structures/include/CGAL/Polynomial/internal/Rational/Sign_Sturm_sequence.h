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

#ifndef CGAL_SIGN_STURM_SEQUENCE_H
#define CGAL_SIGN_STURM_SEQUENCE_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Sign_variations_counter.h>
#include <vector>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template<class Sturm_sequence_t>
class Sign_Sturm_sequence : public Sturm_sequence_t
{
    public:
        typedef Sturm_sequence_t                       Sturm_sequence;
        typedef typename Sturm_sequence::Kernel        Kernel;
        typedef typename Sturm_sequence::Polynomial    Polynomial;
        typedef int                                    result_type;

    protected:
        typedef Sturm_sequence                         Base;
        typedef CGAL::Sign                    Sign;
        typedef typename Kernel::Sign_at               Sign_at;

    public:
        Sign_Sturm_sequence() : Base() {}
        Sign_Sturm_sequence(const Polynomial& p, const Polynomial& q,
            const Kernel &k)
            : Base(p, k.differentiate_object()(p) * q, k),
            pder_(k.differentiate_object()(p)), q_(q) {}

    protected:
        template<class NTRep>
            unsigned int sign_variations_base(const NTRep& x) const
        {
	  Sign s0 = Sign_at( )(this->seq_[0] , x);

            CGAL_Polynomial_precondition( s0 != CGAL::ZERO );

            std::vector<Sign> signs(this->size_);
            signs[0] = s0;

            if ( this->size_ > 1 ) {
	      Sign s1 = Sign_at(  )(pder_, x);
	      Sign s2 = Sign_at( )(q_,  x);
                signs[1] = s1 * s2;
            }

            for (unsigned int i = 2; i < this->size_; i++) {
	      signs[i] = Sign_at( )(this->seq_[i] , x);
            }
            return Sign_variations_counter::sign_variations(signs.begin(), signs.end());
        }

        template<class NTRep>
            int sum_of_signs_base(const NTRep& a, const NTRep& b) const
        {
            CGAL_Polynomial_precondition( b >= a );

            unsigned int Va = sign_variations_base(a);
            if ( Va == 0 ) { return 0; }

            unsigned int Vb = sign_variations_base(b);

            int diff = static_cast<int>(Va) - static_cast<int>(Vb);

            return diff;
        }

    public:

        template<class T>
            unsigned int sign_variations(const T& x) const
        {
            return sign_variations_base(x);
        }

// the following operator() should go away; the only reason it is
// there is to be able to apply the sum_of_signs method to an
// interval using the apply functions (see Isolating_interval.h)
        template<class T>
            int operator()(const T& a, const T& b) const
        {
            return sum_of_signs(a, b);
        }

        template<class T>
            int sum_of_signs(const T& a, const T& b) const
        {
            return sum_of_signs_base(a, b);
        }

    protected:
        Polynomial pder_, q_;
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif                                            // CGAL_SIGN_STURM_SEQUENCE_H
