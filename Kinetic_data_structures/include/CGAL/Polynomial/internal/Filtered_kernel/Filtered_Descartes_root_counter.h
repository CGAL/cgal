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

#ifndef CGAL_POLYNOMIAL_FILTERED_DESCARTES_ROOT_COUNTER_H
#define CGAL_POLYNOMIAL_FILTERED_DESCARTES_ROOT_COUNTER_H

#include <CGAL/Polynomial/basic.h>

#include <CGAL/Polynomial/internal/Alternation_counter.h>
#include <CGAL/Polynomial/internal/interval_arithmetic.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {;

//! A class for computing filtered exact and approximate descartes root bounds.
/*!
  This can be used in two modes, first, just used filtering to perform an exact
  descartes root bound faster. The second will return as soon as it can figure
  out that there is a root in the interval.  This second mode can then return
  the value SOME for the number roots meaning there is some unknown number
  greater than 0 (at most).

  A descartes root bound if performed exact computes whether there are 0, 1 or
  an odd number or an even number of roots. The presence of an odd number of
  roots can be more quickly determined by looking at the signs of the ends of
  the interval (assuming they are non-zero). This bound does not do this
optimization since the values can probably be cached by the calling function
(since a pair of intervals share an endpoint).

The filtering is used in a variety of ways. First, calculation are initially
performed using intervals. If no exact answer is required and there is at
least one definite alternation in the signs of the interval coefficients, then
SOME is returned. If the intervals do not result in sufficient information
begin computed, then the calculations fall back to repeating the process with
an exact number type. Only ring operations are used, so any exact number type
is fine.

I think I should be able to return more exact bounds even if the intervals are
not clear on their sign.  For example:
+?- the value of the ? does not matter
+?+ there are 0 or 2 from the ?, so the parity is still known
+??+ there are at most 1, so this is useless
+??- there are at most 3, so this is useless.
*/

template <class BNT, class Kernel>
unsigned int filtered_Descartes_root_counter(const typename Kernel::Function &fh,
const BNT &begin,
const BNT &end,
bool need_exact,
Kernel k)
{
    typename Kernel::Interval_kernel::NT bi, ei;
    typename Kernel::Interval_kernel::Function fim;
    POLYNOMIAL_NS::Alternation_counter<typename Kernel::Interval_kernel::NT> ac;
    {
        Interval_arithmetic_guard iag;

        if (fh.interval_function().is_constant()) {
//std::cerr << "Why are you solving constant function?\n";
            return 0;
        }

        bi= k.interval_function_converter_object().nt_converter()(begin);
        ei= k.interval_function_converter_object().nt_converter()(end);

//! Could add their optimization here if shift is positive
        fim = k.interval_kernel_object().map_interval_to_positive_object(fh.interval_function())( bi, ei);

        for (int i=0; i<= fim.degree(); ++i) {
            ac.push_back(fim[i]);
            if (!need_exact && ac.number_of_alternations() >0) {
//std::cout << "Quick out in root count.\n";
                return 1;
            }
        }
    }

    if (!ac.is_uncertain()) {
/*std::cout << "No uncertainty in root count: " << fh << " [" << begin << ", " << end << "] is "
  << ac.number_of_alternations() << "\n";*/
        return ac.number_of_alternations();
    }
    else {
/*if (ac.number_of_alternations() >0 && !ac.parity_uncertain()){
  if (ac.number_of_alternations() %2==0) return Descartes_root_count::even();
  else return Descartes_root_count::odd();
  } else*/
//typename CGAL::NT_converter<BNT, typename Kernel::Exact_kernel::Function::NT> ec;

        typename Kernel::Exact_function_converter efc= k.exact_function_converter_object();
        typename Kernel::Exact_kernel::NT be= efc.nt_converter()(begin);
        typename Kernel::Exact_kernel::NT ee= efc.nt_converter()(end);
        typename Kernel::Exact_kernel ek=k.exact_kernel_object();

        typename Kernel::Exact_kernel::Function fem= ek.map_interval_to_positive_object(fh.exact_function())(be, ee);

        POLYNOMIAL_NS::Alternation_counter<typename Kernel::Exact_kernel::NT> ac;
        for (int i=0; i<= fem.degree(); ++i) {
            ac.push_back(fem[i]);
        }
//std::cout << "Exact return in root count.\n";
        return ac.number_of_alternations();
    }
}


template <class Kernel>
class Filtered_Descartes_root_counter
{
    public:
        Filtered_Descartes_root_counter(){}

        Filtered_Descartes_root_counter(const typename Kernel::Function &fh, Kernel k= Kernel()): h_(fh), kernel_(k) {
        }

        ~Filtered_Descartes_root_counter() {
        }

        typedef unsigned int result_type;

        template <class NTT>
            result_type operator()(const NTT &begin, const NTT &end,
            POLYNOMIAL_NS::Sign,
            POLYNOMIAL_NS::Sign) const
        {
            return filtered_Descartes_root_counter(h_, begin, end, true, kernel_);
        }

        template <class NTT>
            result_type operator()(const NTT &begin, const NTT &end) const
        {
            return filtered_Descartes_root_counter(h_, begin, end, true, kernel_);
        }

    protected:
        typename Kernel::Function h_;
        Kernel kernel_;
};

} } } //namespace CGAL::POLYNOMIAL::internal;
#endif
