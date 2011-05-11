#include <typeinfo>
#include <CGAL/myeigen.h>
#include <CGAL/Kernel_d/Cartesian_base.h>
#include <CGAL/Kernel_d/Cartesian_filter_NT.h>
#include <CGAL/Kernel_d/Cartesian_filter_K.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Interval_nt.h>
#include <iostream>
typedef CGAL::Cartesian_base_d<double,CGAL::Dimension_tag<2> > K0;
#if 0
typedef CGAL::Cartesian_filter_NT<K0> K1;
#elif 1
typedef CGAL::Cartesian_base_d<CGAL::Interval_nt_advanced,CGAL::Dimension_tag<2> > KA;
typedef CGAL::Cartesian_base_d<CGAL::Gmpq,CGAL::Dimension_tag<2> > KE;
typedef CGAL::Cartesian_filter_K<K0,KA,KE> K1;
#endif
typedef K1::Point P;
typedef K1::Vector V;
typedef K1::Segment S;
//typedef K1::Construct_point CP;
typedef K1::Construct<CGAL::Construct_point_tag>::type CP;
typedef K1::Construct_vector CV;
typedef K1::Construct_cartesian_const_iterator CCI;
typedef K1::Predicate<CGAL::Orientation_tag>::type PO;
typedef K1::Cartesian_const_iterator CI;
typedef K1::Compute_cartesian_coordinate CC;

int main(){
	CP cp; CCI ci; CC cc; PO po;
	P a=cp(3,4);
	P b=cp(5,6,7);
	int rr[]={3,5,2};
	int* r=rr;
	P c=cp(r,r+2);
	P d=cp(r,r+3,CGAL::Homogeneous_tag());
	S s(c,d);
	std::cout << cc(a,1) << std::endl;
	std::cout << cc(b,1) << std::endl;
	std::cout << cc(s[0],1) << std::endl;
	std::cout << cc(s[1],1) << std::endl;
	for(CI i=ci.begin(a);i!=ci.end(a);++i)
		std::cout << *i << ' ';
	std::cout << '\n';
	P tab[]={a,b,c};
	std::cout << po(&tab[0],tab+3) << std::endl;
}
