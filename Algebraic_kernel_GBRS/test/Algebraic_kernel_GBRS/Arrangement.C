#include <CGAL/Arr_poly_traits_1.h>
#include <CGAL/Arrangement_2.h>
#ifdef __TEST_ARR
#include <ctime>
#else
#include "../../../../CGAL-3.2.1/examples/Arrangement_2/arr_print.h"
#endif

typedef CGAL::GBRS_algebraic_kernel<CGAL::Gmpz>	AlgKernel;
typedef AlgKernel::Coefficient			Coefficient;
typedef AlgKernel::Polynomial_1			Polynomial;
typedef CGAL::Arr_poly_traits_1<AlgKernel>	Traits_2;
typedef Traits_2::X_monotone_curve_2		Curve;
typedef Traits_2::Point_2			Point;
typedef CGAL::Arrangement_2<Traits_2>		Arrangement_2;

int main()
{
	AlgKernel ker;
	Arrangement_2 arr;

	// construct the polynomial x^2-4
	std::vector<Coefficient>coefsp;
	for(int i=0;i<51;++i)
		coefsp.push_back(Coefficient((i%2)?(-i):i));	// x^i
	//--------------------------------------------------
	// coefsp.push_back(Coefficient(-4));	// x^0
	// coefsp.push_back(Coefficient(0));	// x^1
	// coefsp.push_back(Coefficient(1));	// x^2
	//-------------------------------------------------- 
	Polynomial p=ker.construct_polynomial_1_object()
		(coefsp.begin(),coefsp.end());
	// now we create the polynomial q = x-1
	std::vector<Coefficient> coefsq;
	coefsq.push_back(Coefficient(-1));
	coefsq.push_back(Coefficient(1));
	Polynomial q=ker.construct_polynomial_1_object()
		(coefsq.begin(),coefsq.end());

	// r(x)=x^3+1
	std::vector<Coefficient> coefsr;
	coefsr.push_back(Coefficient(1));
	coefsr.push_back(Coefficient(0));
	coefsr.push_back(Coefficient(0));
	coefsr.push_back(Coefficient(1));
	Polynomial r=ker.construct_polynomial_1_object()
		(coefsr.begin(),coefsr.end());
	// s(x)=2
	std::vector<Coefficient> coefss;
	coefss.push_back(Coefficient(2));
	Polynomial s=ker.construct_polynomial_1_object()
		(coefss.begin(),coefss.end());

	Curve cv0=Curve(p,-4,5,"cv0");	// x^2-4 [-4,5]
	Curve cv1=Curve(q,-5,4);	// x-1 [-5,4]
	Curve cv2=Curve(r,-6,6);	// x^3+1 [-6,6]
	Curve cv3=Curve(s,-5,5);	// 2 [-5,5]

#ifdef __TEST_ARR
	clock_t start,end;
	start=clock();
#endif
	insert_curve(arr,cv0);
	insert_curve(arr,cv1);
	insert_curve(arr,cv2);
	insert_curve(arr,cv3);

#ifdef __TEST_ARR
	end=clock();
	std::cout<<"time: "<<((double)(end-start))/CLOCKS_PER_SEC<<" seconds"<<std::endl;
#else
	print_arrangement(arr);
	std::cout<<"is valid? "<<arr.is_valid()<<std::endl;
#endif

	return (0);
}
