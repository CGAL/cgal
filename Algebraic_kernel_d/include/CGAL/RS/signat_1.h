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

#ifndef CGAL_RS_SIGNAT_1_H
#define CGAL_RS_SIGNAT_1_H

#include <CGAL/Gmpfi.h>
#include <CGAL/Polynomial_traits_d.h>
#include "exact_signat_1.h"
//#include <boost/mpl/assert.hpp>
#include <gmp.h>

namespace CGAL{
namespace RS_AK1{

template <class Polynomial_,class Bound_>
struct Signat_1{
        typedef Polynomial_                                     Polynomial;
        typedef Bound_                                          Bound;
        typedef CGAL::Polynomial_traits_d<Polynomial>           PT;
        typedef typename PT::Degree                             Degree;
        Polynomial pol;
        Signat_1(const Polynomial &p):pol(p){};
        CGAL::Sign operator()(const Bound&)const;
}; // struct Signat_1

template <class Polynomial_,class Bound_>
inline CGAL::Sign
Signat_1<Polynomial_,Bound_>::operator()(const Bound_ &x)const{
        typedef Bound_                                          Bound;
        typedef Real_embeddable_traits<Bound>                   REtraits;
        typedef typename REtraits::Sgn                          BSign;
        //typedef Algebraic_structure_traits<Bound>               AStraits;
        // This generic signat works only when Bound_ is an exact type. For
        // non-exact types, an implementation must be provided.
        //BOOST_MPL_ASSERT((boost::is_same<AStraits::Is_exact,Tag_true>));
        int d=Degree()(pol);
        Bound h(pol[d]);
        for(int i=1;i<=d;++i)
                h=h*x+pol[d-i];
        return BSign()(h);
}

template <>
inline CGAL::Sign
Signat_1<Polynomial<Gmpz>,Gmpfr>::operator()(const Gmpfr &x)const{
        // In 32-bit systems, using Gmpfr arithmetic to perform exact
        // evaluations can overflow. For that reason, we only use Gmpfr
        // arithmetic in 64-bit systems.
#if (GMP_LIMB_BITS==64)
        typedef ExactSignat_1<Polynomial,Gmpfr>                 Exact_sign;
#else
        typedef Signat_1<Polynomial,Gmpq>                       Exact_sign;
#endif
        // This seems to work faster for small polynomials:
        // return Exact_sign(pol)(x);
        int d=Degree()(pol);
        if(d==0)
                return pol[0].sign();
        Gmpfi h(pol[d],x.get_precision()+2*d);
        Uncertain<CGAL::Sign> indet=Uncertain<CGAL::Sign>::indeterminate();
        if(h.sign().is_same(indet))
                return Exact_sign(pol)(x);
        for(int i=1;i<=d;++i){
                h*=x;
                h+=pol[d-i];
                if(h.sign().is_same(indet))
                        return Exact_sign(pol)(x);
        }
        CGAL_assertion(!h.sign().is_same(indet));
        return h.sign();
}

// This is the same code as above.
template <>
inline CGAL::Sign
Signat_1<Polynomial<Gmpq>,Gmpfr>::operator()(const Gmpfr &x)const{
        typedef Signat_1<Polynomial,Gmpq>                       Exact_sign;
        int d=Degree()(pol);
        if(d==0)
                return pol[0].sign();
        Gmpfi h(pol[d],x.get_precision()+2*d);
        Uncertain<CGAL::Sign> indet=Uncertain<CGAL::Sign>::indeterminate();
        if(h.sign().is_same(indet))
                return Exact_sign(pol)(x);
        for(int i=1;i<=d;++i){
                h*=x;
                h+=pol[d-i];
                if(h.sign().is_same(indet))
                        return Exact_sign(pol)(x);
        }
        CGAL_assertion(!h.sign().is_same(indet));
        return h.sign();
}

} // namespace RS_AK1
} // namespace CGAL

#endif // CGAL_RS_SIGNAT_1_H
