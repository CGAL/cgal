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

#ifndef CGAL_RS_POLYNOMIAL_CONVERTER_1_H
#define CGAL_RS_POLYNOMIAL_CONVERTER_1_H

namespace CGAL{
namespace RS_AK1{

template <class InputPolynomial_,class OutputPolynomial_>
struct Polynomial_converter_1:
public CGAL::unary_function<InputPolynomial_,OutputPolynomial_>{
        typedef InputPolynomial_                        InpPolynomial_1;
        typedef OutputPolynomial_                       OutPolynomial_1;
        OutPolynomial_1 operator()(const InpPolynomial_1&)const;
}; // class Polynomial_converter_1

template <>
Polynomial<Gmpz>
Polynomial_converter_1<Polynomial<Gmpq>,Polynomial<Gmpz> >::operator()(
                const Polynomial<Gmpq> &p)const{
        std::vector<Gmpz> outcoeffs;
        unsigned degree=p.degree();
        mpz_t lcm;
        mpz_init(lcm);
        mpz_lcm(lcm,mpq_denref(p[0].mpq()),mpq_denref(p[degree].mpq()));
        for(unsigned i=1;i<degree;++i)
                mpz_lcm(lcm,lcm,mpq_denref(p[i].mpq()));
        for(unsigned i=0;i<=degree;++i){
                Gmpz c;
                mpz_divexact(c.mpz(),lcm,mpq_denref(p[i].mpq()));
                mpz_mul(c.mpz(),c.mpz(),mpq_numref(p[i].mpq()));
                outcoeffs.push_back(c);
        }
        mpz_clear(lcm);
        return Polynomial<Gmpz>(outcoeffs.begin(),outcoeffs.end());
}

} // namespace RS_AK1
} // namespace CGAL

#endif // CGAL_RS_POLYNOMIAL_CONVERTER_1_H
