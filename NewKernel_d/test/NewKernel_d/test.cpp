#include <typeinfo>
#include <CGAL/myeigen.h>
#include <CGAL/Kernel_d/Cartesian_base.h>
#include <CGAL/Kernel_d/Cartesian_static_filters.h>
#include <CGAL/Kernel_d/Cartesian_filter_NT.h>
#include <CGAL/Kernel_d/Cartesian_filter_K.h>
#include <CGAL/Kernel_d/Lazy_cartesian.h>
#include <CGAL/Kernel_d/Define_segment.h>
#include <CGAL/Kernel_d/Define_kernel_types.h>
#include <CGAL/Kernel_d/Wrapper/Cartesian_wrap.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Interval_nt.h>
#include <iostream>
//typedef CGAL::Cartesian_base_d<double,CGAL::Dimension_tag<2> > K0;
//typedef CGAL::Cartesian_base_d<CGAL::Interval_nt_advanced,CGAL::Dimension_tag<2> > KA;
struct KA : CGAL::Cartesian_static_filters<CGAL::Dimension_tag<2>, CGAL::Define_segment<CGAL::Cartesian_base_d<CGAL::Interval_nt_advanced,CGAL::Dimension_tag<2>,KA>, KA>, KA> {};
typedef CGAL::Define_segment<CGAL::Cartesian_base_d<CGAL::Gmpq,CGAL::Dimension_tag<2> > > KE;

struct RC: public
CGAL::Cartesian_static_filters<CGAL::Dimension_tag<2>, // Yes, it is silly to put it there.
CGAL::Cartesian_complete_predicates<
CGAL::Cartesian_complete_constructors<
CGAL::Cartesian_complete_computes<
CGAL::Cartesian_complete_types<
CGAL::Cartesian_refcount<
CGAL::Cartesian_LA_base_d<double,CGAL::Dimension_tag<2> >
>
>, false, RC
>, false, RC
>, false, RC
>, RC
>
{
	RC(){}
	RC(int){}
};

typedef CGAL::Define_segment<RC> K0;


#if 0
typedef K0 K2;
#elif 0
typedef CGAL::Cartesian_filter_NT<K0> K2;
#elif 1
typedef CGAL::Cartesian_filter_K<K0,KA,KE> K2;
#elif 1
typedef CGAL::Lazy_cartesian<KE,KA,CGAL::KernelD_converter<KE,KA> > K2;
#endif

#if 0
typedef K2 KK;
#elif 1
typedef CGAL::Cartesian_wrap<K2> KK;
#elif 1
typedef CGAL::Cartesian_wrap<K2> K3;
typedef CGAL::Cartesian_wrap<K3> KK;
#endif

#if 0
typedef KK K1;
typedef K1::Type<CGAL::Point_tag>::type P;
typedef K1::Iterator<CGAL::Point_cartesian_const_iterator_tag>::type CI;
typedef K1::Type<CGAL::Vector_tag>::type V;
typedef K1::Type<CGAL::Segment_tag>::type S;
#elif 1
//typedef CGAL::Define_kernel_types<KK> K1;
typedef KK K1;
typedef K1::Point P;
typedef K1::Point_cartesian_const_iterator CI;
typedef K1::Vector V;
typedef K1::Segment S;
#endif

//typedef K1::Construct_point CP;
typedef K1::Functor<CGAL::Construct_ttag<CGAL::Point_tag> >::type CP;
typedef K1::Functor<CGAL::Construct_ttag<CGAL::Vector_tag> >::type CV;
typedef K1::Functor<CGAL::Construct_ttag<CGAL::Segment_tag> >::type CS;
typedef K1::Functor<CGAL::Segment_extremity_tag>::type CSE;
typedef K1::Functor<CGAL::Construct_ttag<CGAL::Point_cartesian_const_iterator_tag> >::type CCI;
typedef K1::Functor<CGAL::Orientation_of_points_tag>::type PO;
typedef K1::Functor<CGAL::Side_of_oriented_sphere_tag>::type SOS;
//typedef K1::Point_cartesian_const_iterator CI;
typedef K1::Functor<CGAL::Compute_point_cartesian_coordinate_tag>::type CC;

#if 1
#define Kinit (k)
#else
#define Kinit
#endif

int main(){
	K1 k;
	CP cp Kinit;
	CCI ci Kinit;
	CC cc Kinit;
	PO po Kinit;
	CS cs Kinit;
	CSE cse Kinit;
	SOS sos Kinit;
	P a=cp(3,4);
	P b=cp(5,6,7);
	int rr[]={3,5,2};
	int* r=rr;
	P c=cp(r,r+2);
	P d=cp(r,r+3,CGAL::Homogeneous_tag());
	S s=cs(c,d);
	std::cout << cc(a,1) << std::endl;
	std::cout << cc(b,1) << std::endl;
	std::cout << cc(cse(s,0),1) << std::endl;
	std::cout << cc(cse(s,1),1) << std::endl;
	for(CI i=ci(a,CGAL::Begin_tag());i!=ci(a,CGAL::End_tag());++i)
		std::cout << *i << ' ';
	std::cout << '\n';
	P tab[]={a,b,c,d};
	std::cout << po (&tab[0],tab+3) << std::endl;
	std::cout << sos(&tab[0],tab+4) << std::endl;
}
