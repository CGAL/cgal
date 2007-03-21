//! \file examples/Arrangement_2/ex_rational_functions.C
// Constructing an arrangement of arcs of rational functions.
#include <CGAL/basic.h>

#ifdef CGAL_USE_CGAL_CORE
#define CGAL_USE_CORE 1
#endif

#ifndef CGAL_USE_CORE
#include <iostream>
int main ()
{
	std::cout << "Sorry, this example needs CORE ..." << std::endl; 
	return (0);
}
#else

#ifdef __TEST_ARR
#include <ctime>
#endif

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_rational_arc_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
typedef Nt_traits::Rational                           Rational;
typedef Nt_traits::Algebraic                          Algebraic;
typedef CGAL::Cartesian<Algebraic>                    Alg_kernel;
typedef CGAL::Arr_rational_arc_traits_2<Alg_kernel,
	Nt_traits>    Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Curve_2                             Rational_arc_2;
typedef Traits_2::Rat_vector                          Rat_vector;
typedef std::list<Rational_arc_2>                     Rat_arcs_list;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;

int main ()
{
	// y = x^2-4 [-4, 5]:
	Rat_vector P1(51);
	for(int i=0;i<51;++i)
		P1[i]=(i%2)?-i:i;
	//--------------------------------------------------
	// Rat_vector        P1(3);  
	// P1[2] = 1; P1[1] = 0; P1[0] = -4;
	//-------------------------------------------------- 
	Rational_arc_2    a1 (P1, Algebraic(-4), Algebraic(5));

	// y = x-1 [-5, 4]:
	Rat_vector        P2(2);  
	P2[1] = 1; P2[0] = -1;
	Rational_arc_2    a2 (P2, Algebraic(-5), Algebraic(4));

	// y = x^3+1 [-6, 6]:
	Rat_vector        P3(4);  
	P3[3] = 1; P3[2] = 0; P3[1] = 0; P3[0] = 1;
	Rational_arc_2    a3 (P3, Algebraic(-6), Algebraic(6));

	// y = 2 [-5, 5]:
	Rat_vector        P4(2);  
	P4[0] = 2;
	Rational_arc_2    a4 (P4, Algebraic(-5), Algebraic(5));

	// Construct the arrangement of the four arcs.
	Arrangement_2              arr;
	std::list<Rational_arc_2>  arcs;

	arcs.push_back (a1);
	arcs.push_back (a2);
	arcs.push_back (a3);
	arcs.push_back (a4);
#ifdef __TEST_ARR
	clock_t start,end;
	start=clock();
#endif
	insert_curves (arr, arcs.begin(), arcs.end());
#ifdef __TEST_ARR
	end=clock();
	std::cout<<"time: "<<((double)(end-start))/CLOCKS_PER_SEC<<" seconds"<<std::endl;
#else
	// Print the arrangement size.
	std::cout << "The arrangement size:" << std::endl
		<< "   V = " << arr.number_of_vertices()
		<< ",  E = " << arr.number_of_edges() 
		<< ",  F = " << arr.number_of_faces() << std::endl;
#endif

	return 0;
}

#endif
