// Copyright (c) 2006 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
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
// 
//
// Author(s)     : Luis Peñaranda <penarand@loria.fr>

#include <CGAL/Gbrs_algebraic_kernel.h>
#include <gmp.h>
#include <vector>

typedef CGAL::GBRS_algebraic_kernel<CGAL::Gmpz>	AlgKernel;
typedef AlgKernel::Coefficient			Coefficient;
typedef AlgKernel::Algebraic_real_1		Algebraic;
typedef AlgKernel::Polynomial_1			Polynomial;

int main () {
	AlgKernel ker;
	// construct the polynomial x^3-2x
	std::vector<Coefficient> coefsp;
	coefsp.push_back (Coefficient (0));	// x^0
	coefsp.push_back (Coefficient (-2));	// x^1
	coefsp.push_back (Coefficient (0));	// x^2
	coefsp.push_back (Coefficient (1));	// x^3
	coefsp.push_back (Coefficient (0));	// x^4
	coefsp.push_back (Coefficient (0));	// x^5
	// the container has now all the coefficients in increasing monomial
	// order: <0, -2, 0, 1, 0, 0>
	Polynomial p = ker.construct_polynomial_1_object()
		(coefsp.begin (), coefsp.end ());
	std::cout << "p(x) = " << p << std::endl;
	// we test some overloaded operators:
	std::cout << "p(x) * (int)3 = " << p*3 << std::endl;
	mpz_t x;
	mpz_init (x);
	mpz_set_ui (x, 3);
	std::cout << "p(x) * (mpz_t)3 = " << p*x << std::endl;
	mpz_clear (x);

	// now we create the polynomial q = x-1
	std::vector<Coefficient> coefsq;
	coefsq.push_back (Coefficient(-1414));
	coefsq.push_back (Coefficient(1000));
	Polynomial q = ker.construct_polynomial_1_object()
		(coefsq.begin (), coefsq.end ());
	std::cout << "\nq(x) = " << q << std::endl;

	std::cout<<"p*q = "<<p*q;
	std::cout<<"\np-q = "<<p-q;
	std::cout<<"\np' = "<<p.derive();
	std::cout<<"\nq' = "<<q.derive();
	std::cout<<"\np*=q = "<<(p*=q);
	std::cout<<"\nq*=2 = "<<(q*=2)<<std::endl;

	double aaa=8.5;
	Algebraic alg(aaa);
	std::cout<<"alg="<<alg<<" is root of "<<alg.pol()<<std::endl;

	return 0;
}

