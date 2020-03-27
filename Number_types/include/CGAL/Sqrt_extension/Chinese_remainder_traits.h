// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer <hemmer@informatik.uni-mainz.de>

#ifndef CGAL_SQRT_EXTENSION_CHINESE_REMAINDER_TRAITS
#define CGAL_SQRT_EXTENSION_CHINESE_REMAINDER_TRAITS

#include <CGAL/basic.h>
#include <CGAL/Chinese_remainder_traits.h>
#include <CGAL/Sqrt_extension/Sqrt_extension_type.h>

namespace CGAL {

template <class NT, class ROOT, class ACDE_TAG,class FP_TAG>
class Chinese_remainder_traits<CGAL::Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG> >{
private:
  typedef CGAL::Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG> EXT;
    typedef Chinese_remainder_traits<NT> CRT_NT;
    typedef Chinese_remainder_traits<ROOT> CRT_ROOT;

public :
    typedef EXT Type;
    typedef typename CRT_NT::Scalar_type Scalar_type;

    struct Chinese_remainder{
        void operator() (
                const Scalar_type& m1, const Scalar_type& m2,
                const Scalar_type& m,
                const Scalar_type& s,  const Scalar_type& t,
                const Type& u1, const Type& u2,
                Type& u) const
        {
            NT   a0,a1;
            ROOT root;
            typename CRT_NT::Chinese_remainder chinese_remainder_nt;
            chinese_remainder_nt(m1,m2,m,s,t,u1.a0(),u2.a0(),a0);

            CGAL_postcondition_msg(!((u1.is_extended() && !u2.is_extended())
                            || (!u1.is_extended() && u2.is_extended())),
                    " chinese remainder not possible for these numbers");

            if(u1.is_extended() || u2.is_extended()){
                chinese_remainder_nt(m1,m2,m,s,t,u1.a1(),u2.a1(),a1);
                typename CRT_ROOT::Chinese_remainder chinese_remainder_root;
                chinese_remainder_root(m1,m2,m,s,t,u1.root(),u2.root(),root);
                u=Type(a0,a1,root);
            }else{
                u=Type(a0);
            }
        }
        void operator() (
                const Scalar_type& m1, const Type& u1,
                const Scalar_type& m2, const Type& u2,
                Scalar_type& m, Type& u) const {
            Scalar_type s,t;
            CGAL::extended_euclidean_algorithm(m1,m2,s,t);
            m = m1 * m2;
            this->operator()(m1,m2,m,s,t,u1,u2,u);
        }
    };
};

} // namespace CGAL

#endif // CGAL_SQRT_EXTENSION_CHINESE_REMAINDER_TRAITS
