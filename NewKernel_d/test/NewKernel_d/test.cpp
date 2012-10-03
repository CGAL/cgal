#define BOOST_RESULT_OF_USE_DECLTYPE 1
#include <CGAL/Epick_d.h>
#include <typeinfo>
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

#if 1
#define Kinit(f) =k.f()
#else
#define Kinit(f)
#endif

template<class Ker>
void test2(){
  typedef Ker K1;
  typedef typename K1::Point_d P;
  typedef typename K1::Cartesian_const_iterator_d CI;
  typedef typename K1::Vector_d V;
  typedef typename K1::Segment_d S;
  typedef typename K1::Flat_orientation_d FO;

  //typedef K1::Construct_point CP;
  typedef typename K1::Construct_point_d CP;
  typedef typename K1::Construct_vector_d CV;
  typedef typename K1::Construct_segment_d CS;
  typedef typename K1::template Functor<CGAL::Segment_extremity_tag>::type CSE;
  typedef typename K1::Construct_cartesian_const_iterator_d CCI;
  typedef typename K1::Orientation_d PO;
  typedef typename K1::Side_of_oriented_sphere_d SOS;
  typedef typename K1::Compute_coordinate_d CC;
  typedef typename K1::Construct_flat_orientation_d CFO;
  typedef typename K1::In_flat_orientation_d IFO;
  typedef typename K1::In_flat_side_of_oriented_sphere_d IFSOS;
  typedef typename K1::Contained_in_affine_hull_d CAH;
  typedef typename K1::Compare_lexicographically_d CL;

  USE_TYPE(V);
  USE_TYPE(CV);
  USE_TYPE(FO);
  USE_TYPE(CL);
  Ker k
#if 0
    (2)
#endif
    ;
  CP cp Kinit(construct_point_d_object);
  CCI ci Kinit(construct_cartesian_const_iterator_d_object);
  CC cc Kinit(compute_coordinate_d_object);
  PO po Kinit(orientation_d_object);
  CS cs Kinit(construct_segment_d_object);
  CSE cse (k);
  SOS sos Kinit(side_of_oriented_sphere_d_object);
  CFO cfo Kinit(construct_flat_orientation_d_object);
  IFO ifo Kinit(in_flat_orientation_d_object);
  IFSOS ifsos Kinit(in_flat_side_of_oriented_sphere_d_object);
  CAH cah Kinit(contained_in_affine_hull_d_object);
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
  for(CI i=ci(a);i!=ci(a,0);++i)
    std::cout << *i << ' ';
  std::cout << '\n';
  P tab[]={a,b,c,d};
  std::cout << po (&tab[0],tab+3) << std::endl;
  std::cout << sos(&tab[0],tab+4) << std::endl;
}

template<class CP> struct Construct_point3_helper {
  CP const& cp;
  Construct_point3_helper(CP const& x) : cp(x) {}
  template<class T1,class T2,class T3>
  typename CP::result_type operator()(T1 const&t1, T2 const&t2, T3 const&t3)const{
    double tab[]={(double)t1,(double)t2,(double)t3};
    return cp(tab+0,tab+3);
  }
  template<class T1,class T2,class T3,class T4>
  typename CP::result_type operator()(T1 const&t1, T2 const&t2, T3 const&t3, T4 const&t4)const{
    double tab[]={(double)t1,(double)t2,(double)t3};
    return cp(tab+0,tab+3,t4);
  }
};

template<class Ker>
void test3(){
  typedef Ker K1;
  typedef typename K1::Point_d P;
  typedef typename K1::Cartesian_const_iterator_d CI;
  typedef typename K1::Vector_d V;
  typedef typename K1::Segment_d S;
  typedef typename K1::Flat_orientation_d FO;

  //typedef K1::Construct_point CP;
  typedef typename K1::Construct_point_d CP_;
  typedef typename K1::Construct_vector_d CV;
  typedef typename K1::Construct_segment_d CS;
  typedef typename K1::template Functor<CGAL::Segment_extremity_tag>::type CSE;
  typedef typename K1::Construct_cartesian_const_iterator_d CCI;
  typedef typename K1::Orientation_d PO;
  typedef typename K1::Side_of_oriented_sphere_d SOS;
  typedef typename K1::Compute_coordinate_d CC;
  typedef typename K1::Construct_flat_orientation_d CFO;
  typedef typename K1::In_flat_orientation_d IFO;
  typedef typename K1::In_flat_side_of_oriented_sphere_d IFSOS;
  typedef typename K1::Contained_in_affine_hull_d CAH;
  typedef typename K1::Compare_lexicographically_d CL;

  USE_TYPE(V);
  USE_TYPE(CV);
  USE_TYPE(FO);
  Ker k
#if 1
    (3)
#endif
    ;
  CP_ cp_ Kinit(construct_point_d_object);
  typename boost::conditional<boost::is_same<typename Ker::Default_ambient_dimension,CGAL::Dynamic_dimension_tag>::value,Construct_point3_helper<CP_>,CP_>::type cp(cp_);
  CCI ci Kinit(construct_cartesian_const_iterator_d_object);
  CC cc Kinit(compute_coordinate_d_object);
  CL cl Kinit(compare_lexicographically_d_object);
  PO po Kinit(orientation_d_object);
  CS cs Kinit(construct_segment_d_object);
  CSE cse (k);
  SOS sos Kinit(side_of_oriented_sphere_d_object);
  CFO cfo Kinit(construct_flat_orientation_d_object);
  IFO ifo Kinit(in_flat_orientation_d_object);
  IFSOS ifsos Kinit(in_flat_side_of_oriented_sphere_d_object);
  CAH cah Kinit(contained_in_affine_hull_d_object);
  P a=cp(2,3,4);
  P b=cp(5,6,7,8);
  int rr[]={3,5,2,3};
  int* r=rr;
  P c=cp_(3,r,r+3);
  P d=cp_(r,r+4,CGAL::Homogeneous_tag());
  S s=cs(c,d);
  std::cout << cc(a,1) << std::endl;
  std::cout << cc(b,2) << std::endl;
  std::cout << cc(cse(s,0),1) << std::endl;
  std::cout << cc(cse(s,1),2) << std::endl;
  for(CI i=ci(a);i!=ci(a,0);++i)
    std::cout << *i << ' ';
  std::cout << '\n';
  P e=cp(-2,3,0);
  P tab[]={a,b,c,d,e};
  std::cout << po (&tab[0],tab+4) << std::endl;
  std::cout << sos(&tab[0],tab+5) << std::endl;
  FO fo=cfo(&tab[0],tab+3);
  std::cout << fo;
  P x[]={cp(2,2,3),cp(2,2,0),cp(1,2,1)};
  FO fo2=cfo(&x[0],x+3);
  std::cout << fo2;
  P y[]={cp(0,2,4),cp(3,1,2),cp(3,3,6),cp(0,4,8)};
  FO fo3=cfo(&y[0],y+3);
  assert(fo3.rest.size()==1 && fo3.rest[0]!=3);
  std::cout << fo3;
  CGAL::Orientation base=ifo(fo3,&y[0],y+3);
  assert(ifo(fo3,y+1,y+4)==base);
  P yy[]={y[1],y[3],y[0],y[2]};
  assert(ifo(fo3,yy+0,yy+3)==base);
  assert(ifo(fo3,yy+1,yy+4)==base);
  std::cout << ifsos(fo3,y+0,y+3,y[3]) << ' ';
  std::cout << ifsos(fo3,y+1,y+4,y[0]) << ' ';
  std::cout << ifsos(fo3,yy+0,yy+3,yy[3]) << ' ';
  std::cout << ifsos(fo3,yy+1,yy+4,yy[0]) << '\n';
  P buf[]={cp(100,900,0),y[0],y[1],y[2],y[3]};
  std::cout << sos(buf+0,buf+5) << ' ';
  buf[1]=y[1];buf[2]=y[2];buf[3]=y[3];buf[4]=y[0];
  std::cout << sos(buf+0,buf+5) << ' ';
  buf[1]=yy[0];buf[2]=yy[1];buf[3]=yy[2];buf[4]=yy[3];
  std::cout << sos(buf+0,buf+5) << ' ';
  buf[1]=yy[1];buf[2]=yy[2];buf[3]=yy[3];buf[4]=yy[0];
  std::cout << sos(buf+0,buf+5) << '\n';
  assert(cah(y+0,y+3,y[3]));
  assert(!cah(y+0,y+3,buf[0]));
  assert(cl(a,a)==CGAL::EQUAL);
  assert(cl(a,b)==CGAL::LARGER);
}

int main(){
  test2<CGAL::Kernel_d_interface<KK> >();
  test2<CGAL::Epick_d<2> >();
  test3<CGAL::Epick_d<3> >();
  test3<CGAL::Epick_d<CGAL::UNKNOWN_DIMENSION> >();
}
