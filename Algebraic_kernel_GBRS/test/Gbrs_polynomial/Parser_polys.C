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
// Author(s)     : Elias Tsigaridas <Elias.Tsigaridas@loria.fr>
//                 Luis Peñaranda <penarand@loria.fr>

#include <CGAL/Arr_poly_traits_1.h>
#include <CGAL/Arrangement_2.h>
#ifdef __TEST_ARR
#include <ctime>
#else
#include "../../../../CGAL-3.2.1/examples/Arrangement_2/arr_print.h"
#endif
#include "parsers.h"

typedef CGAL::GBRS_algebraic_kernel<CGAL::Gmpz>	AlgKernel;
typedef AlgKernel::Coefficient			Coefficient;
typedef AlgKernel::Polynomial_1			Polynomial;
typedef CGAL::Arr_poly_traits_1<AlgKernel>	Traits_2;
typedef Traits_2::X_monotone_curve_2		Curve;
typedef Traits_2::Point_2			Point;
typedef CGAL::Arrangement_2<Traits_2>		Arrangement_2;

Polynomial parse_poly(AlgKernel ker,std::string tstr){
	CGAL::Polynomial_parser_1 parser;
	parser.parse(tstr);
	CGAL_assertion(parser.is_correct());
	std::vector<Coefficient> Coeff;
	// We get the coefficients from the parser.
	// Notice that we use our convertor. If we didn't supply a convertor
	// then the default convertor would be used (to ints)
	parser.result(std::back_inserter(Coeff),CGAL::The_Convert_to());
	return ker.construct_polynomial_1_object()(Coeff.begin(),Coeff.end());
}

int main(){
	AlgKernel ker;
	Arrangement_2 arr;
	std::string s;
	unsigned n;
	std::cout<<"number of polynomials in the arrangement: ";
	std::cin>>n;
	Polynomial p[n];
	Coefficient left[n],right[n];
	for(unsigned i=0;i<n;++i){
		std::cout<<"\npolynomial? ";
		std::cin>>s;
		p[i]=parse_poly(ker,s);
		std::cout<<"left bound? ";
		std::cin>>s;
		left[i]=s;
		std::cout<<"right bound? ";
		std::cin>>s;
		right[i]=s;
		CGAL_assertion(left[i]<right[i]);
	}
#ifdef __TEST_ARR
	clock_t start,end;
	start=clock();
#endif
	for(unsigned i=0;i<n;++i)
		insert_curve(arr,Curve(p[i],left[i],right[i]));
#ifdef __TEST_ARR
	end=clock();
	std::cout<<"\ntime: "<<((double)(end-start))/CLOCKS_PER_SEC<<" seconds, ";
#else
	print_arrangement(arr);
#endif
	std::cout<<"valid="<<arr.is_valid()<<std::endl;
	return (0);
}
