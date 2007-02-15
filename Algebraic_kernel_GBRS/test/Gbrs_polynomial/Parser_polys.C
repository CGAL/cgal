#include <CGAL/Arr_poly_traits_1.h>
#include <CGAL/Arrangement_2.h>
#include "../../../../CGAL-3.2.1/examples/Arrangement_2/arr_print.h"
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
	Polynomial p;
	if(parser.is_correct()){
		// The coefficients of the polynomial
		std::vector<Coefficient> Coeff;
		// We get the coefficients from the parser.
		// Notice that we use our convertor. If we didn't supply a convertor
		// then the default convertor would be used (to ints)
		parser.result(std::back_inserter(Coeff),CGAL::The_Convert_to());

		// Output the result
		//std::cout.precision(20);

		/*std::cout << "Coeff: ";
		std::copy(Coeff.begin(),Coeff.end(),std::ostream_iterator<Coefficient>(std::cout," "));
		std::cout << std::endl; 
	*/
	p=ker.construct_polynomial_1_object()
		(Coeff.begin(),Coeff.end());
	} else {
		// Something went wrong
		std::cout << "Failure..." << std::endl; 
		// The error was at...
		std::cout << "at: " << parser.error() << std::endl; 
	}
	return p;
}

int main()
{
	AlgKernel ker;
	Arrangement_2 arr;

	/*Polynomial p=parse_poly(ker,"x^2-4");
	Polynomial q=parse_poly(ker,"x-1");
	Polynomial r=parse_poly(ker,"x^3+1");
	Polynomial s=parse_poly(ker,"2");

	Curve cv0=Curve(p,-4,5);	// x^2-4 [-4,5]
	Curve cv1=Curve(q,-5,4);	// x-1 [-5,4]
	Curve cv2=Curve(r,-6,6);	// x^3+1 [-6,6]
	Curve cv3=Curve(s,-5,5);	// 2 [-5,5]

	insert_curve(arr,cv0);
	insert_curve(arr,cv1);
	insert_curve(arr,cv2);
	insert_curve(arr,cv3);*/

	Polynomial q=parse_poly(ker,"x^200-4");
	Polynomial r=parse_poly(ker,"x^2");
	Curve cvq=Curve(q,-10,10);
	Curve cvr=Curve(r,-10,10);
	insert_curve(arr,cvq);
	insert_curve(arr,cvr);

	print_arrangement(arr);
	std::cout<<"is valid? "<<arr.is_valid()<<std::endl;

	return (0);
}
