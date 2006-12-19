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

#include <CGAL/Arr_poly_traits_1.h>
#include <gmp.h>
#include <vector>

typedef CGAL::GBRS_algebraic_kernel<CGAL::Gmpz>	AlgKernel;
typedef AlgKernel::Coefficient			Coefficient;
typedef AlgKernel::Algebraic_real_1		Algebraic;
typedef AlgKernel::Polynomial_1			Polynomial;

typedef CGAL::GBRS_algebraic_kernel<CGAL::Gmpz>	AlgKernel;
typedef CGAL::Arr_poly_traits_1<AlgKernel>	Traits;
typedef Traits::X_monotone_curve_2		Curve;
typedef Traits::Point_2				Point;
typedef Traits::Multiplicity			Mult;

int main () {
	AlgKernel ker;
	// construct the polynomial x^2-4
	std::vector<Coefficient> coefsp;
	coefsp.push_back(Coefficient(-4));	// x^0
	coefsp.push_back(Coefficient(0));	// x^1
	coefsp.push_back(Coefficient(1));	// x^2
	Polynomial p=ker.construct_polynomial_1_object()
		(coefsp.begin(),coefsp.end());

	// now we create the polynomial q = x-1
	std::vector<Coefficient> coefsq;
	coefsq.push_back(Coefficient(-1));
	coefsq.push_back(Coefficient(1));
	Polynomial q=ker.construct_polynomial_1_object()
		(coefsq.begin(),coefsq.end());

	// we have to intersect them
	Traits t;
	std::vector<std::pair<Point,Mult> > points;
	Curve cv1(p,-4,5);
	Curve cv2(q,-5,4);
	std::cout<<"Curves are:\n"<<cv1<<"\n"<<cv2<<std::endl;
	t.intersect_2_object()(cv1,cv2,std::back_inserter(points));

	std::cout<<"\ncv1 is ";
	if(!t.is_vertical_2_object()(cv1))
		std::cout<<"not ";
	std::cout<<"vertical"<<std::endl;
	Point cv1_left(t.construct_min_vertex_2_object()(cv1));
	std::cout<<"left point of cv1:"<<cv1_left<<std::endl;
	std::cout<<"right point of cv1:"<<
		t.construct_max_vertex_2_object()(cv1)<<std::endl;

	std::cout<<"\nintersection points:\n";
	for(std::vector<std::pair<Point,Mult> >::iterator itv=points.begin();
			itv!=points.end();++itv)
		itv->first.show(std::cout)<<" mult "<<itv->second<<std::endl;

	std::cout<<"\ncomparison of the intersection points:"<<
		t.compare_xy_2_object()(points[0].first,points[1].first)<<
		std::endl;

	std::cout<<"\ncv1's left endpoint is ";
	switch(t.compare_y_at_x_2_object()(cv1_left,cv2)){
		case CGAL::EQUAL:std::cout<<"in";break;
		case CGAL::LARGER:std::cout<<"over";break;
		default:std::cout<<"below";break;
	}
	std::cout<<" cv2"<<std::endl;

	std::cout<<"\ncompare_y_at_left_2(cv1,cv2,p[0])="<<
		t.compare_y_at_x_left_2_object()(cv1,cv2,points[0].first)<<
		std::endl;

	std::cout<<"compare_y_at_right_2(cv1,cv2,p[0])="<<
		t.compare_y_at_x_right_2_object()(cv1,cv2,points[0].first)<<
		std::endl;

	std::cout<<"\np[0]==p[1]? "<<t.equal_2_object()(points[0].first,points[1].first)<<std::endl;
	std::cout<<"p[1]==p[0]? "<<t.equal_2_object()(points[1].first,points[0].first)<<std::endl;
	std::cout<<"p[0]==p[0]? "<<t.equal_2_object()(points[0].first,points[0].first)<<std::endl;
	std::cout<<"attention: EQUAL is "<<CGAL::EQUAL<<std::endl;

	Curve split1,split2;
	t.split_2_object()(cv1,points[0].first,split1,split2);
	std::cout<<"\nsplitting cv1 in the first point gives:\n"<<split1<<
		"\nand "<<split2<<std::endl;

	std::vector<Curve> monotones;
	t.make_x_monotone_2_object()(cv1,std::back_inserter(monotones));
	std::cout<<"splitting cv1 in monotone pieces gives:"<<std::endl;
	for(std::vector<Curve>::iterator it=monotones.begin();
			it!=monotones.end();++it)
		std::cout<<*it<<std::endl;

	return 0;
}
