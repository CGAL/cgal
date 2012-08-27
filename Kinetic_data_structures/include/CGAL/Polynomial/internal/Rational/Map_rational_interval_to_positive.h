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

#ifndef CGAL_POLYNOMIAL_INTERNAL_MAP_INTERVAL_TO_POSITIVE_H
#define CGAL_POLYNOMIAL_INTERNAL_MAP_INTERVAL_TO_POSITIVE_H

#include <CGAL/Polynomial/basic.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

//------------------------------------------------------------------

template <class K>
class Map_rational_interval_to_positive
{
    public:
        typedef typename K::Function Polynomial;
        typedef typename Polynomial::NT NT;
        typedef NT first_argument_type;
        typedef NT second_argument_type;
        typedef Polynomial result_type;
        Map_rational_interval_to_positive(){}
        Map_rational_interval_to_positive(const Polynomial &f, const K &k): f_(f), k_(k){}

        Polynomial operator()(const NT &lb, const NT& ub) const
        {
            typedef typename K::Rational_translate_zero  Rational_translate_zero;

//T_a
            Rational_translate_zero tr= k_.rational_translate_zero_object(lb);
            Polynomial t0= tr(f_);

            NT diff= ub-lb;
            NT pdiff= diff;

//H_{b-a}
            int t0_size = t0.degree()+1;
            std::vector<NT> t1_coef(t0_size);

            t1_coef[0] = t0[0];
            for (int i=1; i < t0_size; ++i) {
                t1_coef[i] = pdiff*t0[i];
                pdiff= pdiff*diff;
            }
            Polynomial t1(t1_coef.begin(), t1_coef.end());

// R
            typename K::Invert_variable iv= k_.invert_variable_object();
            Polynomial t2= iv(t1);

// T_1
            Rational_translate_zero tr2 = k_.rational_translate_zero_object(1);

            Polynomial t3 = tr2(t2);
            return t3;
        }
    protected:
        Polynomial f_;
        K k_;
};

template <class K>
class Map_rational_interval_to_positive_2
{
    public:
        typedef typename K::Function Polynomial;
        typedef typename Polynomial::NT NT;
        typedef Polynomial argument_type;
        typedef Polynomial result_type;
        Map_rational_interval_to_positive_2(){}
        Map_rational_interval_to_positive_2(const NT &a, const NT &b, const K &k): lb_(a), ub_(b), k_(k){}

        Polynomial operator()(const Polynomial &f) const
        {
            typedef typename K::Rational_translate_zero  Rational_translate_zero;

//T_a
            Rational_translate_zero tr= k_.rational_translate_zero_object(lb_);
            Polynomial t0= tr(f);

            NT diff= ub_-lb_;
            NT pdiff= diff;

//H_{b-a}
            int t0_size = t0.degree()+1;
            std::vector<NT> t1_coef(t0_size);

            t1_coef[0] = t0[0];
            for (int i=1; i < t0_size; ++i) {
                t1_coef[i] = pdiff*t0[i];
                pdiff= pdiff*diff;
            }
            Polynomial t1(t1_coef.begin(), t1_coef.end());

// R
            typename K::Invert_variable iv= k_.invert_variable_object();
            Polynomial t2= iv(t1);

// T_1
            Rational_translate_zero tr2 = k_.rational_translate_zero_object(1);

            Polynomial t3 = tr2(t2);
            return t3;
        }
    protected:
        NT lb_, ub_;
        K k_;
};
} } } //namespace CGAL::POLYNOMIAL::internal
#endif
