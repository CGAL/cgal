#define BOOST_RESULT_OF_USE_DECLTYPE 1
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
#include <CGAL/Kernel_d/Kernel_d_interface.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Interval_nt.h>
#include <iostream>

template<class>void marc_use(){}
#define USE_TYPE(T) marc_use<T>()

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

template<class Ker>
void test2(){
  typedef Ker K1;
  typedef typename K1::Point_d P;
  typedef typename K1::Cartesian_const_iterator_d CI;
  typedef typename K1::Vector_d V;
  typedef typename K1::Segment_d S;

#if 1
#define Kinit(f) =k.f()
#else
#define Kinit(f)
#endif

  //typedef K1::Construct_point CP;
  typedef typename K1::Construct_point_d CP;
  typedef typename K1::Construct_vector_d CV;
  typedef typename K1::Construct_segment_d CS;
  typedef typename K1::template Functor<CGAL::Segment_extremity_tag>::type CSE;
  typedef typename K1::Construct_cartesian_const_iterator_d CCI;
  typedef typename K1::Orientation_d PO;
  typedef typename K1::Side_of_oriented_sphere_d SOS;
  typedef typename K1::Compute_coordinate_d CC;

  USE_TYPE(V);
  USE_TYPE(CV);
  Ker k;
  CP cp Kinit(construct_point_d_object);
  CCI ci Kinit(construct_cartesian_const_iterator_d_object);
  CC cc Kinit(compute_coordinate_d_object);
  PO po Kinit(orientation_d_object);
  CS cs Kinit(construct_segment_d_object);
  CSE cse (k);
  SOS sos Kinit(side_of_oriented_sphere_d_object);
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
int main(){
  test2<CGAL::Kernel_d_interface<KK> >();
}
