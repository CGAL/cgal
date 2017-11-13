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

#ifndef CGAL_POLYNOMIAL_INTERNAL_FILTERED_SIGN_ATQ_H
#define CGAL_POLYNOMIAL_INTERNAL_FILTERED_SIGN_ATQ_H

#include <CGAL/Polynomial/basic.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template <class Kernel, class NTT>
Sign filtered_sign_at(const typename Kernel::Function &fh,
NTT &t, const Kernel &k)
{
    typename Kernel::Interval_kernel::Function::NT i,v;
    {
        Interval_arithmetic_guard iag;

        if (fh.interval_function().degree()<1) {
            return ZERO;
        }

        i= k.interval_function_converter_object().nt_converter()(t);
        v= fh.interval_function().value_at(i);
    }

//std::cout << "Interval value is " << v << std::endl;

    if (v.inf() >0) return POSITIVE;
    else if (v.sup() < 0) return NEGATIVE;
    else if (v.inf()==v.sup()) {
        return CGAL::sign(v);
    }
    else {
//std::cout << "Falling back on exact.\n";

//typename CGAL::NT_converter<NTT, typename Kernel::Exact_kernel::Function::NT> ec;
        typename Kernel::Exact_kernel::NT et= k.exact_function_converter_object().nt_converter()(t);

        return CGAL::sign(fh.exact_function().value_at(et));
    };
}


template <class Kernel>
class Filtered_sign_at
{
    typedef Filtered_sign_at<Kernel> This;
    public:

        Filtered_sign_at(){}

        Filtered_sign_at(const typename Kernel::Function &fh, Kernel k= Kernel()): h_(fh), kernel_(k) {
        }

        ~Filtered_sign_at() {
        }

        typedef Sign result_type;
//! not really, but close enough
        typedef typename Kernel::Function::NT argument_type;

        template <class NTT>
            result_type operator()(const NTT &t) const
        {
            return filtered_sign_at(h_, t, kernel_);
        }

    protected:
        typename Kernel::Function h_;
        Kernel kernel_;
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
