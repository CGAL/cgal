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

#ifndef CGAL_POLYNOMIAL_INTERNAL_STURM_ROOT_COUNTER_H
#define CGAL_POLYNOMIAL_INTERNAL_STURM_ROOT_COUNTER_H

#include <CGAL/Polynomial/basic.h>

#include <CGAL/Polynomial/internal/Sign_variations_counter.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template<class Kernel>
struct Sturm_root_counter
{
    public:

        typedef typename Kernel::Function::NT argument_type;
        typedef typename Kernel::Function::NT first_argument_type;
        typedef typename Kernel::Function::NT second_argument_type;
        typedef unsigned int            result_type;

    protected:
        typedef typename Kernel::Function Polynomial;
        typedef typename Kernel::Function::NT   NT;

        template<class NTRep>
            unsigned int number_of_real_roots_base(const NTRep& a,
            const NTRep& b) const
        {
            CGAL_precondition( b >= a );

            unsigned int Va = sign_variations_near(a, true);
            if ( Va == 0 ) { return 0; }

            unsigned int Vb = sign_variations_near(b, false);

//CGAL_assertion( Va > Vb );

            return Va - Vb;
        }

        template<class Iterator>
            static
            unsigned int sign_variations(const Iterator& first,
        const Iterator& beyond) {
            return Sign_variations_counter::sign_variations(first, beyond);
        }

        template<class NTRep>
            unsigned int sign_variations_near(const NTRep& x, bool above) const
        {
//CGAL::Sign s0 = k_.sign_at_object( sseq[0] )(x);

//CGAL_exactness_precondition( s0 != CGAL::ZERO );

// MK:: We need to optimize that; if f(a) is not zero then I don't
// care about zeros later on in the sequence.

            std::vector<Sign> signs(sseq.size());

	    //typename K::Sign_at sa= k_.sign_at_object();
            for (unsigned int i = 0; i < sseq.size(); i++) {
	      if (above) {
		signs[i] = k_.sign_above_object( sseq[i] )(x);
	      }
	      else {
		signs[i] = k_.sign_below_object( sseq[i] )(x);
	      }
		//signs[i]= sa(sseq[i], x);
            }

            return sign_variations(signs.begin(), signs.end());
        }

    public:
        Sturm_root_counter():sseq(){}
        Sturm_root_counter(const typename Kernel::Function& p,
            const Kernel& k)
            : k_(k), sseq(k_.Sturm_sequence_object(p, k_.differentiate_object()(p))) {}

        Sturm_root_counter(const typename Kernel::Standard_sequence& sseq,
            Kernel k= Kernel())
            : k_(k), sseq(sseq) {}

        template<class T>
            result_type
            operator()(const T& a) const
        {
            return sseq.sign_variations(a);
        }

        template <class T>
            result_type
            operator()(const T& a, const T& b) const
        {
            return number_of_real_roots_base(a, b);
        }

    protected:
        Kernel k_;
        typename Kernel::Sturm_sequence sseq;
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif                                            // CGAL_STURM_ROOT_COUNTER_H
