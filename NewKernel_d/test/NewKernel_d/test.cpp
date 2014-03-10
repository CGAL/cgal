//#define BOOST_RESULT_OF_USE_DECLTYPE 1
#include <CGAL/Epick_d.h>
#include <typeinfo>
#include <CGAL/Kernel_d/Cartesian_base.h>
#include <CGAL/Kernel_d/Cartesian_static_filters.h>
#include <CGAL/Kernel_d/Cartesian_filter_NT.h>
#include <CGAL/Kernel_d/Cartesian_filter_K.h>
#include <CGAL/Kernel_d/Lazy_cartesian.h>
#include <CGAL/Kernel_d/Wrapper/Cartesian_wrap.h>
#include <CGAL/Kernel_d/Kernel_d_interface.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Interval_nt.h>
#include <iostream>

template<class>void marc_use(){}
#define USE_TYPE(T) marc_use<T>()

//typedef CGAL::Cartesian_base_d<double,CGAL::Dimension_tag<2> > K0;
//typedef CGAL::Cartesian_base_d<CGAL::Interval_nt_advanced,CGAL::Dimension_tag<2> > KA;
struct KA : CGAL::Cartesian_static_filters<CGAL::Dimension_tag<2>, CGAL::Cartesian_base_d<CGAL::Interval_nt_advanced,CGAL::Dimension_tag<2>, KA>, KA> {};
typedef CGAL::Cartesian_base_d<CGAL::Gmpq,CGAL::Dimension_tag<2> > KE;

struct RC: public
CGAL::Cartesian_static_filters<CGAL::Dimension_tag<2>, // Yes, it is silly to put it there.
 CGAL::Cartesian_refcount<
  CGAL::Cartesian_LA_base_d<double,CGAL::Dimension_tag<2> >
 >, RC
>
{
	RC(){}
	RC(int){}
};

struct K0 : RC {};


#if 0
typedef K0 K2;
#elif 0
typedef CGAL::Cartesian_filter_NT<K0> K2;
#elif 0
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
  typedef typename K1::Hyperplane_d H;
  typedef typename K1::Sphere_d Sp;
  typedef typename K1::Vector_d V;
  typedef typename K1::Segment_d S;
  typedef typename K1::Flat_orientation_d FO;

  //typedef K1::Construct_point CP;
  typedef typename K1::Construct_point_d CP;
  typedef typename K1::Construct_vector_d CV;
  typedef typename K1::Construct_segment_d CS;
  typedef typename K1::Construct_sphere_d CSp;
  typedef typename CGAL::Get_functor<K1, CGAL::Segment_extremity_tag>::type CSE;
  typedef typename K1::Construct_cartesian_const_iterator_d CCI;
  typedef typename K1::Linear_base_d LB;
  typedef typename K1::Orientation_d PO;
  typedef typename K1::Side_of_oriented_sphere_d SOS;
  typedef typename K1::Side_of_bounded_sphere_d SBS;
  typedef typename K1::Compute_coordinate_d CC;
  typedef typename K1::Construct_flat_orientation_d CFO;
  typedef typename K1::In_flat_orientation_d IFO;
  typedef typename K1::In_flat_side_of_oriented_sphere_d IFSOS;
  typedef typename K1::Contained_in_affine_hull_d CAH;
  typedef typename K1::Compare_lexicographically_d CL;
  typedef typename K1::Value_at_d VA;
  typedef typename K1::Construct_hyperplane_d CH;
  typedef typename K1::Center_of_sphere_d COS;
  typedef typename K1::Affine_rank_d AR;
  typedef typename K1::Linear_rank_d LR;
  typedef typename K1::Affinely_independent_d AI;
  typedef typename K1::Linearly_independent_d LI;
  typedef typename K1::Has_on_positive_side_d HOPS;
  typedef typename K1::Less_coordinate_d LC;
  typedef typename K1::Less_lexicographically_d LL;
  typedef typename K1::Less_or_equal_lexicographically_d LEL;
  typedef typename K1::Midpoint_d M;
  typedef typename K1::Oriented_side_d OS;
  typedef typename K1::Orthogonal_vector_d OV;
  typedef typename K1::Point_dimension_d PD;
  typedef typename K1::Point_of_sphere_d PS;
  typedef typename K1::Point_to_vector_d PV;
  typedef typename K1::Vector_to_point_d VP;
  typedef typename K1::Barycentric_coordinates_d BC;
  typedef typename K1::Construct_direction_d CD;
  typedef typename K1::Construct_line_d CLi;
  typedef typename K1::Construct_ray_d CR;
  typedef typename K1::Construct_iso_box_d CIB;
  typedef typename K1::Construct_aff_transformation_d CAT;
  typedef typename K1::Position_on_line_d PoL;

  // FIXME: really test everything at least once (clang warnings can list untested things).
  USE_TYPE(FO);
  Ker k
#if 0
    (2)
#endif
    ;
  CP cp Kinit(construct_point_d_object);
  CV cv Kinit(construct_vector_d_object);
  CCI ci Kinit(construct_cartesian_const_iterator_d_object);
  CC cc Kinit(compute_coordinate_d_object);
  PO po Kinit(orientation_d_object);
  CS cs Kinit(construct_segment_d_object);
  CSp csp Kinit(construct_sphere_d_object);
  VA va Kinit(value_at_d_object);
  CH ch Kinit(construct_hyperplane_d_object);
  CSE cse (k);
  SOS sos Kinit(side_of_oriented_sphere_d_object);
  SBS sbs Kinit(side_of_bounded_sphere_d_object);
  CFO cfo Kinit(construct_flat_orientation_d_object);
  IFO ifo Kinit(in_flat_orientation_d_object);
  IFSOS ifsos Kinit(in_flat_side_of_oriented_sphere_d_object);
  CAH cah Kinit(contained_in_affine_hull_d_object);
  LB lb Kinit(linear_base_d_object);
  COS cos Kinit(center_of_sphere_d_object);
  LR lr Kinit(linear_rank_d_object);
  AR ar Kinit(affine_rank_d_object);
  AI ai Kinit(affinely_independent_d_object);
  LI li Kinit(linearly_independent_d_object);
  HOPS hops Kinit(has_on_positive_side_d_object);
  LC lc Kinit(less_coordinate_d_object);
  CL cl Kinit(compare_lexicographically_d_object);
  LL ll Kinit(less_lexicographically_d_object);
  LEL lel Kinit(less_or_equal_lexicographically_d_object);
  M m Kinit(midpoint_d_object);
  OS os Kinit(oriented_side_d_object);
  OV ov Kinit(orthogonal_vector_d_object);
  PD pd Kinit(point_dimension_d_object);
  PS ps Kinit(point_of_sphere_d_object);
  PV pv Kinit(point_to_vector_d_object);
  VP vp Kinit(vector_to_point_d_object);
  BC bc Kinit(barycentric_coordinates_d_object);
  CD cd Kinit(construct_direction_d_object);
  CLi cli Kinit(construct_line_d_object);
  CR cr Kinit(construct_ray_d_object);
  CIB cib Kinit(construct_iso_box_d_object);
  CAT cat Kinit(construct_aff_transformation_d_object);
  PoL pol Kinit(position_on_line_d_object);

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
  std::cout << po (&tab[0],tab+3) << ' ';
  std::cout << sos(&tab[0],tab+4) << ' ';
  std::cout << sbs(&tab[0],tab+4) << std::endl;
  P x1=cp(0,1);
  P x2=cp(-1,-1);
  P x3=cp(1,-1);
  P x4=cp(0,0);
  P x5=cp(0,-1);
  P tab2[]={x1,x2,x3,x4};
  assert(po(tab2+0,tab2+3)==CGAL::COUNTERCLOCKWISE);
  assert(sos(tab2+0,tab2+3,x4)==CGAL::ON_POSITIVE_SIDE);
  assert(sbs(tab2+0,tab2+3,x4)==CGAL::ON_BOUNDED_SIDE);
  V y1=cp(1,-1);
  V y2=cp(3,-3);
  P tab3[]={y1,y2};
  std::vector<V> v;
  std::back_insert_iterator<std::vector<V> > bii(v);
  lb(tab3+0,tab3+2,bii);
  assert(v.size()==1);
  H h=ch(tab2+1,tab2+3);
  assert(fabs(va(h,x2)-1)<.0001);
  assert(fabs(va(h,x3)-1)<.0001);
  assert(fabs(va(h,x1)+1)<.0001);
#if 1
  // Doesn't compile with Lazy yet.
  FO fo=cfo(tab2+1,tab2+3);
  assert(ifo(fo,tab2+1,tab2+3)==CGAL::POSITIVE);
  assert(ifsos(fo,tab2+1,tab2+3,x5)==CGAL::ON_POSITIVE_SIDE);
#endif
  P z0=cp( 0+2,5-3);
  P z1=cp(-5+2,0-3);
  P z2=cp( 3+2,4-3);
  P tabz[]={z0,z1,z2};
  Sp sp = csp(tabz+0,tabz+3);
  P cent0=cos(sp);
  P cent1=cos(tabz+0,tabz+3);
  assert(fabs(cent0[0]-2)<.0001);
  assert(fabs(cent0[1]+3)<.0001);
  assert(fabs(cent1[0]-2)<.0001);
  assert(fabs(cent1[1]+3)<.0001);
  assert(fabs(sp.squared_radius()-25)<.0001);
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
  typedef typename CGAL::Get_functor<K1, CGAL::Segment_extremity_tag>::type CSE;
  typedef typename K1::Construct_cartesian_const_iterator_d CCI;
  typedef typename K1::Orientation_d PO;
  typedef typename K1::Side_of_oriented_sphere_d SOS;
  typedef typename K1::Side_of_bounded_sphere_d SBS;
  typedef typename K1::Compute_coordinate_d CC;
  typedef typename K1::Construct_flat_orientation_d CFO;
  typedef typename K1::In_flat_orientation_d IFO;
  typedef typename K1::In_flat_side_of_oriented_sphere_d IFSOS;
  typedef typename K1::Contained_in_affine_hull_d CAH;
  typedef typename K1::Contained_in_linear_hull_d CLH;
  typedef typename K1::Contained_in_simplex_d CiS;
  typedef typename K1::Compare_lexicographically_d CL;
  typedef typename K1::Component_accessor_d CA;
  typedef typename K1::Equal_d E;
  typedef typename K1::Squared_distance_d SD;

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
  SBS sbs Kinit(side_of_bounded_sphere_d_object);
  CFO cfo Kinit(construct_flat_orientation_d_object);
  IFO ifo Kinit(in_flat_orientation_d_object);
  IFSOS ifsos Kinit(in_flat_side_of_oriented_sphere_d_object);
  CAH cah Kinit(contained_in_affine_hull_d_object);
  CLH clh Kinit(contained_in_linear_hull_d_object);
  CiS cis Kinit(contained_in_simplex_d_object);
  CA ca Kinit(component_accessor_d_object);
  E ed Kinit(equal_d_object);
  SD sd Kinit(squared_distance_d_object);
  P a=cp(2,3,4);
  P b=cp(5,6,7,8);
  int rr[]={3,5,2,3};
  int* r=rr;
  P c=cp_(3,r,r+3);
  P d=cp_(r,r+4,CGAL::Homogeneous_tag());
  S s=cs(c,d);
  std::cout << cc(a,1) << std::endl;
  std::cout << cc(b,2) << std::endl;
  std::cout << cse(s,0)[1] << std::endl;
  std::cout << cc(cse(s,1),2) << std::endl;
  for(CI i=ci(a);i!=ci(a,0);++i)
    std::cout << *i << ' ';
  std::cout << '\n';
  P e=cp(-2,3,0);
  P tab[]={a,b,c,d,e};
  std::cout << po (&tab[0],tab+4) << ' ';
  std::cout << sos(&tab[0],tab+5) << ' ';
  std::cout << sbs(&tab[0],tab+5) << std::endl;
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
  P x1=cp(0,1,-1);
  P x2=cp(-1,-1,-1);
  P x3=cp(1,-1,-1);
  P x4=cp(0,0,1);
  P x5=cp(0,0,0);
  P x6=cp(0,0,-1);
  P tab2[]={x1,x2,x3,x4,x5};
  assert(po(tab2+0,tab2+4)==CGAL::POSITIVE);
  assert(sos(tab2+0,tab2+4,x5)==CGAL::ON_POSITIVE_SIDE);
  assert(sbs(tab2+0,tab2+4,x5)==CGAL::ON_BOUNDED_SIDE);
  FO fo4=cfo(tab2+0,tab2+3);
  assert(ifo(fo4,tab2+0,tab2+3)==CGAL::POSITIVE);
  assert(ifsos(fo4,tab2+0,tab2+3,x6)==CGAL::ON_POSITIVE_SIDE);
}
template struct CGAL::Epick_d<CGAL::Dimension_tag<2> >;
template struct CGAL::Epick_d<CGAL::Dimension_tag<3> >;
template struct CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
int main(){
  //Broken with Linear_base_d (output iterator)
  //test2<CGAL::Kernel_d_interface<KK> >();
  test2<CGAL::Epick_d<CGAL::Dimension_tag<2> > >();
  test3<CGAL::Epick_d<CGAL::Dimension_tag<3> > >();
  test3<CGAL::Epick_d<CGAL::Dynamic_dimension_tag> >();
}
