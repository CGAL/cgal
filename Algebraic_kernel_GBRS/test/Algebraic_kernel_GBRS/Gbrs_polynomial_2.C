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

#include <CGAL/Gbrs_polynomial_2.h>

typedef CGAL::Rational_polynomial_2	Polynomial;

int main () {
	Polynomial p (2, 1);	// d_x=2, d_y=1
	Polynomial q (7, 8);
	Polynomial zero (1, 4);	// zero(x) = 0

	p.set_coef (2, 1, 2);	// 2*x^2*y
	p.set_coef (1, 0, 3);	// 3*x
	p.set_coef (0, 0, 7);	// 7
	// p = 2*x^2*y + 3*x + 7
	q = -p;
	p.set_coef (0, 0, 5);
	std::cout << "p = " << p << "\nq = " << q << "\n0 = " << zero << "\n";
	std::cout << "p==q : " << (p==q) << "\np!=q : " << (p!=q) << "\n";

	Polynomial r;
	Polynomial s (1, 0);
	s.set_coef (1, 0, 1);	// s = 1 * x^1
	r = p*s;
	std::cout<<"s = "<<s<<std::endl;
	std::cout<<"r = p*s = "<<r<<std::endl;
	std::cout<<"p*=(s*2) = "<<(p*=(s*2))<<std::endl;
	std::cout << "p+s+q = " << (p+s+q) << std::endl;
	std::cout << "s+p+q = " << (s+p+q) << std::endl;
	std::cout<<"s-=q = "<<(s-=q)<<std::endl;

	return 0;
}

