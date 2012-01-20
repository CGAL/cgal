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

#ifndef CGAL_POLYNOMIAL_INTERNAL_EVALUATE_H
#define CGAL_POLYNOMIAL_INTERNAL_EVALUATE_H
#include <CGAL/export/CGAL.h>
#include <CGAL/Polynomial/basic.h>
#include <vector>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

CGAL_EXPORT double evaluate_polynomial(const double *b, const double *e, double t);

template <class NT>
inline NT evaluate_polynomial(const std::vector<NT> &coefs, const NT &t)
{
    if (coefs.empty()) return NT(0);
    typename std::vector<NT>::const_reverse_iterator rit= coefs.rbegin();
    NT result = *rit;
    ++rit;
    for (; rit != coefs.rend(); ++rit) {
        result *= t;
        result += (*rit);
    }
    return result;
}


inline double evaluate_polynomial(const std::vector<double>& coefs, double t)
{
    if (coefs.empty()) return 0;
    return evaluate_polynomial(&coefs.front(), &coefs.front()+coefs.size(), t);
}


template<class K, class NT >
class Evaluate_polynomial
{
    typedef typename K::Function            Polynomial;

    public:
        Evaluate_polynomial(){}
        Evaluate_polynomial(const Polynomial& p): coefs_(p.begin(), p.end()) {
/*coefs_.resize(p.degree()+1);
for (unsigned int i=0; i< p.degree(); ++i){
  coefs_[i]= NT(p[i]);
  }*/
        }

        typedef NT argument_type;
        typedef NT result_type;

        result_type operator()(const argument_type& x) const
        {
            return evaluate_polynomial(coefs_, x);
        }

    protected:
        std::vector<result_type> coefs_;
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
