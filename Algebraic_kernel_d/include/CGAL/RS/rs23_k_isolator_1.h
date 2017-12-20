// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.

// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_RS23_K_ISOLATOR_1_H
#define CGAL_RS_RS23_K_ISOLATOR_1_H

// This file includes an isolator. Its particularity is that is isolates the
// roots with RS2, and the refines them until reaching Kantorovich criterion.
// This can take long, but the later refinements will be extremely fast with
// RS3. The functor is not in RS2 neither in RS3 namespace, because it uses
// functions from both.

#include "rs2_calls.h"
#include "rs3_k_refiner_1.h"
#include <CGAL/Gmpfi.h>
#include <vector>

namespace CGAL{

template <class Polynomial_,class Bound_>
class RS23_k_isolator_1{
        public:
        typedef Polynomial_                             Polynomial;
        typedef Bound_                                  Bound;
        private:
        typedef Gmpfi                                   Interval;
        public:
        RS23_k_isolator_1(const Polynomial&);
        Polynomial polynomial()const;
        int number_of_real_roots()const;
        bool is_exact_root(int i)const;
        Bound left_bound(int i)const;
        Bound right_bound(int i)const;
        private:
        Polynomial _polynomial;
        std::vector<Interval> _real_roots;
};

template <class Polynomial_,class Bound_>
RS23_k_isolator_1<Polynomial_,Bound_>::
RS23_k_isolator_1(const Polynomial_ &p){
        CGAL_error_msg("not implemented for these polynomial/bound types");
}

template <>
RS23_k_isolator_1<CGAL::Polynomial<CGAL::Gmpz>,Gmpfr>::
RS23_k_isolator_1(const CGAL::Polynomial<CGAL::Gmpz> &p):_polynomial(p){
        typedef CGAL::Polynomial<CGAL::Gmpz>                    Pol;
        typedef CGAL::Gmpfr                                     Bound;
        typedef CGAL::RS3::RS3_k_refiner_1<Pol,Bound>           KRefiner;
        int numsols;
        std::vector<Gmpfi> intervals;
        RS2::RS2_calls::init_solver();
        RS2::RS2_calls::create_rs_upoly(p,rs_get_default_up());
        set_rs_precisol(0);
        set_rs_verbose(0);
        rs_run_algo((char*)"UISOLE");
        RS2::RS2_calls::insert_roots(std::back_inserter(intervals));
        // RS2 computed the isolating intervals. Now, we use RS3 to refine each
        // root until reaching Kantorovich criterion, before adding it to the
        // root vector.
        numsols=intervals.size();
        for(int j=0;j<numsols;++j){
                Gmpfr left(intervals[j].inf());
                Gmpfr right(intervals[j].sup());
                CGAL_assertion(left<=right);
                KRefiner()(p,left,right,53);
                _real_roots.push_back(intervals[j]);
        }
}

template <>
RS23_k_isolator_1<CGAL::Polynomial<CGAL::Gmpq>,Gmpfr>::
RS23_k_isolator_1(const CGAL::Polynomial<CGAL::Gmpq> &qp):_polynomial(qp){
        typedef CGAL::Polynomial<CGAL::Gmpz>                    ZPol;
        typedef CGAL::Gmpfr                                     Bound;
        typedef CGAL::RS3::RS3_k_refiner_1<ZPol,Bound>          ZKRefiner;
        int numsols;
        std::vector<Gmpfi> intervals;
        CGAL::Polynomial<CGAL::Gmpz> zp=CGAL::RS_AK1::Polynomial_converter_1<
                                                CGAL::Polynomial<Gmpq>,
                                                CGAL::Polynomial<Gmpz> >()(qp);
        RS2::RS2_calls::init_solver();
        RS2::RS2_calls::create_rs_upoly(zp,rs_get_default_up());
        set_rs_precisol(0);
        set_rs_verbose(0);
        rs_run_algo((char*)"UISOLE");
        RS2::RS2_calls::insert_roots(std::back_inserter(intervals));
        // RS2 computed the isolating intervals. Now, we use RS3 to refine each
        // root until reaching Kantorovich criterion, before adding it to the
        // root vector.
        numsols=intervals.size();
        for(int j=0;j<numsols;++j){
                Gmpfr left(intervals[j].inf());
                Gmpfr right(intervals[j].sup());
                ZKRefiner()(zp,left,right,53);
                _real_roots.push_back(intervals[j]);
        }
}

template <class Polynomial_,class Bound_>
Polynomial_
RS23_k_isolator_1<Polynomial_,Bound_>::polynomial()const{
        return _polynomial;
}

template <class Polynomial_,class Bound_>
int
RS23_k_isolator_1<Polynomial_,Bound_>::number_of_real_roots()const{
        return _real_roots.size();
}

template <class Polynomial_,class Bound_>
bool
RS23_k_isolator_1<Polynomial_,Bound_>::is_exact_root(int i)const{
        return _real_roots[i].inf()==_real_roots[i].sup();
}

template <class Polynomial_,class Bound_>
Bound_
RS23_k_isolator_1<Polynomial_,Bound_>::left_bound(int i)const{
        return _real_roots[i].inf();
}

template <class Polynomial_,class Bound_>
Bound_
RS23_k_isolator_1<Polynomial_,Bound_>::right_bound(int i)const{
        return _real_roots[i].sup();
}

} // namespace CGAL

#endif // CGAL_RS_RS23_K_ISOLATOR_1_H
