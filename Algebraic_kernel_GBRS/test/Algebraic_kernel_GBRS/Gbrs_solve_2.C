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
// $URL: $
// $Id: $
// 
//
// Author(s)     : Luis Peñaranda <penarand@loria.fr>

#include <CGAL/Gbrs_polynomial_2.h>
#include <CGAL/Gbrs_solve_2.h>

typedef CGAL::Rational_polynomial_2	Polynomial;

int main () {
	// p=1+2y+3x+4xy
	Polynomial p(1,1);
	p.set_coef (0,0,1);
	p.set_coef (0,1,2);
	p.set_coef (1,0,3);
	p.set_coef (1,1,4);

	// q=1-y+x-xy
	Polynomial q(1,1);
	q.set_coef (0,0,1);
	q.set_coef (0,1,-1);
	q.set_coef (1,0,1);
	q.set_coef (1,1,-1);

	std::cout<<"p="<<p<<"\nq="<<q<<std::endl;

	solve_2(p,q);

	return 0;
}

