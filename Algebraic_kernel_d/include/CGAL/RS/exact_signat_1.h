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

#ifndef CGAL_RS_EXACT_SIGNAT_1_H
#define CGAL_RS_EXACT_SIGNAT_1_H

#include <CGAL/Polynomial_traits_d.h>
#include "dyadic.h"

namespace CGAL{
namespace RS_AK1{

template <class Polynomial_,class Bound_>
struct ExactSignat_1{
        typedef Polynomial_                                     Polynomial;
        typedef Bound_                                          Bound;
        typedef CGAL::Polynomial_traits_d<Polynomial>           PT;
        typedef typename PT::Degree                             Degree;
        Polynomial pol;
        ExactSignat_1(const Polynomial &p):pol(p){};
        CGAL::Sign operator()(const Bound&)const;
}; // struct ExactSignat_1

template <>
inline CGAL::Sign
ExactSignat_1<Polynomial<Gmpz>,Gmpfr>::operator()(const Gmpfr &x)const{
        int d=Degree()(pol);
        if(d==0)
                return pol[0].sign();
        // Construct a Gmpfr containing exactly the leading coefficient.
        Gmpfr h(pol[d],pol[d].bit_size());
        CGAL_assertion(h==pol[d]);
        // Horner's evaluation.
        for(int i=1;i<d;++i){
                CGALRS_dyadic_mul(h.fr(),h.fr(),x.fr());
                CGALRS_dyadic_add_z(h.fr(),h.fr(),pol[d-i].mpz());
        }
        // TODO: We can avoid doing the last addition.
        return h.sign();
}

} // namespace RS_AK1
} // namespace CGAL

#endif // CGAL_RS_EXACT_SIGNAT_1_H
