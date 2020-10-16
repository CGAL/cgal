#include <boost/config.hpp>
#if defined(BOOST_GCC) && (__GNUC__ <= 4) && (__GNUC_MINOR__ < 4)

#include <iostream>
int main()
{
  std::cerr << "NOTICE: This test requires G++ >= 4.4, and will not be compiled." << std::endl;
}

#else

//#define BOOST_RESULT_OF_USE_DECLTYPE 1
#include <CGAL/Epick_d.h>
#include <CGAL/Epeck_d.h>
#include <typeinfo>

#include <CGAL/NewKernel_d/Cartesian_base.h>
#include <CGAL/NewKernel_d/Cartesian_static_filters.h>
#include <CGAL/NewKernel_d/Cartesian_filter_NT.h>
#include <CGAL/NewKernel_d/Cartesian_filter_K.h>
#include <CGAL/NewKernel_d/Lazy_cartesian.h>
#include <CGAL/NewKernel_d/Wrapper/Cartesian_wrap.h>
#include <CGAL/NewKernel_d/Kernel_d_interface.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/use.h>
#include <iostream>
#include <sstream>
#include <CGAL/NewKernel_d/Types/Weighted_point.h>

#include <cmath>

//typedef CGAL::Cartesian_base_d<double,CGAL::Dimension_tag<2> > K0;
//typedef CGAL::Cartesian_base_d<CGAL::Interval_nt_advanced,CGAL::Dimension_tag<2> > KA;
struct KA : CGAL::Cartesian_static_filters<CGAL::Dimension_tag<2>, CGAL::Cartesian_base_d<CGAL::Interval_nt_advanced,CGAL::Dimension_tag<2>, KA>, KA> {};
typedef CGAL::Cartesian_base_d<CGAL::Exact_rational,CGAL::Dimension_tag<2> > KE;

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
#elif 1
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
  typedef typename K1::Aff_transformation_d AT;
  typedef typename K1::Direction_d D;
  typedef typename K1::Line_d L;
  typedef typename K1::Ray_d R;
  typedef typename K1::Iso_box_d IB;
  typedef typename K1::Flat_orientation_d FO;
  typedef typename K1::Weighted_point_d WP;

  //typedef K1::Construct_point CP;
  typedef typename K1::Construct_point_d CP;
  typedef typename K1::Construct_vector_d CV;
  typedef typename K1::Construct_segment_d CS;
  typedef typename K1::Construct_sphere_d CSp;
  typedef typename CGAL::Get_functor<K1, CGAL::Segment_extremity_tag>::type CSE;
  typedef typename K1::Construct_cartesian_const_iterator_d CCI;
  typedef typename K1::Construct_circumcenter_d CCc;
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
  typedef typename K1::Has_on_positive_side_d HOPS;
  typedef typename K1::Less_coordinate_d LC;
  typedef typename K1::Less_lexicographically_d LL;
  typedef typename K1::Less_or_equal_lexicographically_d LEL;
  typedef typename K1::Midpoint_d M;
  typedef typename K1::Oriented_side_d OS;
  typedef typename K1::Orthogonal_vector_d OV;
  typedef typename K1::Point_dimension_d PD;
  typedef typename K1::Point_to_vector_d PV;
  typedef typename K1::Vector_to_point_d VP;
  typedef typename K1::Barycentric_coordinates_d BC;
  typedef typename K1::Construct_direction_d CD;
  typedef typename K1::Construct_line_d CLi;
  typedef typename K1::Construct_ray_d CR;
  typedef typename K1::Construct_iso_box_d CIB;
  typedef typename K1::Construct_aff_transformation_d CAT;
  typedef typename K1::Position_on_line_d PoL;
  typedef typename K1::Equal_d E;
  typedef typename K1::Squared_distance_d SD;
  typedef typename K1::Squared_length_d SL;
  typedef typename K1::Scalar_product_d SP;
  typedef typename K1::Difference_of_vectors_d DV;
  typedef typename K1::Difference_of_points_d DP;
  typedef typename K1::Construct_min_vertex_d CmV;
  typedef typename K1::Construct_max_vertex_d CMV;
  typedef typename K1::Compute_squared_radius_d SR;
  typedef typename K1::Translated_point_d TP;
  typedef typename K1::Construct_power_sphere_d PC;
  typedef typename K1::Compute_power_product_d PoD;
  typedef typename K1::Construct_weighted_point_d CWP;
  typedef typename K1::Power_side_of_bounded_power_sphere_d PSBPS;
  typedef typename K1::Compute_squared_radius_smallest_orthogonal_sphere_d CSRSOS;
  //typedef typename K1::Point_drop_weight_d PDW;
  typedef CP PDW;
  typedef typename K1::Compute_weight_d PW;

  CGAL_USE_TYPE(AT);
  CGAL_USE_TYPE(D);
  CGAL_USE_TYPE(L);
  CGAL_USE_TYPE(R);
  CGAL_USE_TYPE(IB);

  Ker k
#if 0
    (2)
#endif
    ;
  CP cp Kinit(construct_point_d_object);
  CV cv Kinit(construct_vector_d_object);
  CCI ci Kinit(construct_cartesian_const_iterator_d_object);
  CC cc Kinit(compute_coordinate_d_object);
  CCc ccc Kinit(construct_circumcenter_d_object);
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
  HOPS hops Kinit(has_on_positive_side_d_object);
  LC lc Kinit(less_coordinate_d_object);
  CL cl Kinit(compare_lexicographically_d_object);
  LL ll Kinit(less_lexicographically_d_object);
  LEL lel Kinit(less_or_equal_lexicographically_d_object);
  M m Kinit(midpoint_d_object);
  OS os Kinit(oriented_side_d_object);
  OV ov Kinit(orthogonal_vector_d_object);
  PD pd Kinit(point_dimension_d_object);
  PV pv Kinit(point_to_vector_d_object);
  VP vp Kinit(vector_to_point_d_object);
  BC bc Kinit(barycentric_coordinates_d_object);
  CD cd Kinit(construct_direction_d_object);
  CLi cli Kinit(construct_line_d_object);
  CR cr Kinit(construct_ray_d_object);
  CIB cib Kinit(construct_iso_box_d_object);
  CAT cat Kinit(construct_aff_transformation_d_object);
  PoL pol Kinit(position_on_line_d_object);
  E ed Kinit(equal_d_object);
  SD sd Kinit(squared_distance_d_object);
  SL sl Kinit(squared_length_d_object);
  SP spr Kinit(scalar_product_d_object);
  DV dv Kinit(difference_of_vectors_d_object);
  DP dp Kinit(difference_of_points_d_object);
  CmV cmv Kinit(construct_min_vertex_d_object);
  CMV cMv Kinit(construct_max_vertex_d_object);
  SR sr Kinit(compute_squared_radius_d_object);
  TP tp Kinit(translated_point_d_object);
  PC pc Kinit(construct_power_sphere_d_object);
  CWP cwp Kinit(construct_weighted_point_d_object);
  //PDW pdw Kinit(point_drop_weight_d_object);
  PDW const& pdw = cp;
  PW pw Kinit(compute_weight_d_object);
  PoD pod Kinit(compute_power_product_d_object);
  PSBPS psbps Kinit(power_side_of_bounded_power_sphere_d_object);
  CSRSOS csrsos Kinit(compute_squared_radius_smallest_orthogonal_sphere_d_object);

  CGAL_USE(bc);
  CGAL_USE(pol);
  CGAL_USE(cat);
  CGAL_USE(cd);
  CGAL_USE(cli);
  CGAL_USE(cr);
  using std::abs;
  P a=cp(3,4);
  assert(pd(a)==2);
  assert(pv(a)[1]==4);
  P b=vp(cv(5,6,7));
  assert(abs(b[0]-5./7)<.0001);
  assert(lc(b,a,1));
  assert(!lc(a,b,0));
  int rr[]={3,5,2};
  int* r=rr;
  P c=cp(r,r+2);
  assert(!ll(a,a));
  assert(lel(a,a));
  assert(cl(a,a)==CGAL::EQUAL);
  assert(ll(a,c));
  assert(!lel(c,a));
  assert(cl(a,c)==CGAL::SMALLER);
  assert(ll(b,c));
  assert(cl(c,b)==CGAL::LARGER);
  assert(abs(m(a,c)[0]-3)<.0001);
  assert(abs(m(a,c)[1]-4.5)<.0001);
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
  assert(po (&tab[0],tab+3) == CGAL::CLOCKWISE);
  std::cout << sos(tab+1,tab+4,a) << ' ';
  std::cout << sbs(tab+1,tab+4,a) << std::endl;
  P tabp[]={P(0,0),P(1,0),P(0,1)};
  P tabn[]={P(0,0),P(0,1),P(1,0)};
  assert(po(tabp+0,tabp+3)==CGAL::POSITIVE);
  assert(po(tabn+0,tabn+3)==CGAL::NEGATIVE);
  assert(sos(tabp+0,tabp+3,P(3,3))==CGAL::ON_NEGATIVE_SIDE);
  assert(sos(tabn+0,tabn+3,P(3,3))==CGAL::ON_POSITIVE_SIDE);
  assert(sbs(tabp+0,tabp+3,P(3,3))==CGAL::ON_UNBOUNDED_SIDE);
  assert(sbs(tabn+0,tabn+3,P(3,3))==CGAL::ON_UNBOUNDED_SIDE);
  assert(sbs(tabp+1,tabp+3,P(1,1))==CGAL::ON_BOUNDARY);
  assert(ccc(tabp+1,tabp+2)==tabp[1]);
  assert(ccc(tabn+0,tabn+2)==P(0,.5));
  assert(sr(tabp+2,tabp+3)==0);
  assert(sr(tabp+1,tabp+3)==.5);
  assert(sbs(tabp+1,tabp+3,P(10,-1))==CGAL::ON_UNBOUNDED_SIDE);
  assert(sos(tabp+0,tabp+3,P(.5,.5))==CGAL::ON_POSITIVE_SIDE);
  assert(sos(tabn+0,tabn+3,P(.5,.5))==CGAL::ON_NEGATIVE_SIDE);
  assert(sbs(tabp+0,tabp+3,P(.5,.5))==CGAL::ON_BOUNDED_SIDE);
  assert(sbs(tabn+0,tabn+3,P(.5,.5))==CGAL::ON_BOUNDED_SIDE);
  P x1=cp(0,1);
  P x2=cp(-1,-1);
  P x3=cp(1,-1);
  P x4=cp(0,0);
  P x5=cp(0,-1);
  P tab2[]={x1,x2,x3,x4};
  assert(dp(x1,x2)[1]==2);
  assert(po(tab2+0,tab2+3)==CGAL::COUNTERCLOCKWISE);
  assert(sos(tab2+0,tab2+3,x4)==CGAL::ON_POSITIVE_SIDE);
  assert(sbs(tab2+0,tab2+3,x4)==CGAL::ON_BOUNDED_SIDE);
  V y1=cv(1,-1);
  assert(y1.squared_length()==2);
  assert(sl(y1)==2);
  V y2=cv(3,-3);
  assert(spr(y1,y2)==6);
  assert(dv(y2,y1)[0]==2);
  V tab3[]={y1,y2};
  std::vector<V> v;
  std::back_insert_iterator<std::vector<V> > bii(v);
  lb(tab3+0,tab3+2,bii);
  assert(v.size()==1);
  assert(lr(tab3+0,tab3+2)==1);
  H h=ch(tab2+1,tab2+3,tab2[0]);
  assert(abs(va(h,x2)-1)<.0001);
  assert(abs(va(h,x3)-1)<.0001);
  assert(abs(va(h,x1)+1)<.0001);
  H h2=ch(tab2+1,tab2+3,x1,CGAL::ON_POSITIVE_SIDE);
  assert(hops(h2,x1));
  assert(os(h2,x1)==CGAL::ON_POSITIVE_SIDE);
  H h3=ch(tab2+1,tab2+3,x1,CGAL::ON_NEGATIVE_SIDE);
  assert(!hops(h3,x1));
  P tab4[]={cp(-1,1),cp(1,-1),cp(-1,-1)};
  H h4=ch(tab4+0,tab4+2,tab4[2],CGAL::ON_POSITIVE_SIDE);
  assert(hops(h4,tab4[2]));
  assert(os(h4,tab4[1])==CGAL::ON_ORIENTED_BOUNDARY);
  V hv=ov(h); CGAL_USE(hv);
#if 1
  // Doesn't compile with Lazy yet.
  FO fo=cfo(tab2+1,tab2+3);
  assert(ifo(fo,tab2+1,tab2+3)==CGAL::POSITIVE);
  assert(ifsos(fo,tab2+1,tab2+3,x5)==CGAL::ON_POSITIVE_SIDE);
  P tab_h[]={P(0,42),P(1,42),P(4,42),P(2,42),P(3,42)};
  assert(cah(tab_h+0,tab_h+2,tab_h[4]));
  P py2=cp(3,-3);
  assert(!cah(tab_h+0,tab_h+2,py2));
  FO fo_hp = cfo (tab_h+0, tab_h+2);
  FO fo_hn = cfo (tab_h+2, tab_h+4);
  assert(ifo(fo_hp, tab_h+1, tab_h+3)==CGAL::POSITIVE);
  assert(ifo(fo_hn, tab_h+1, tab_h+3)==CGAL::NEGATIVE);
  assert(ifo(fo_hn, tab_h+2, tab_h+4)==CGAL::POSITIVE);
  assert(ifo(fo_hp, tab_h+2, tab_h+4)==CGAL::NEGATIVE);
  assert(ifsos(fo_hp, tab_h+1, tab_h+3, tab_h[3])==CGAL::ON_POSITIVE_SIDE);
  assert(ifsos(fo_hn, tab_h+1, tab_h+3, tab_h[3])==CGAL::ON_NEGATIVE_SIDE);
  assert(ifsos(fo_hp, tab_h+0, tab_h+2, tab_h[2])==CGAL::ON_NEGATIVE_SIDE);
  assert(ifsos(fo_hn, tab_h+0, tab_h+2, tab_h[2])==CGAL::ON_POSITIVE_SIDE);
  assert(ifsos(fo_hp, tab_h+2, tab_h+4, tab_h[1])==CGAL::ON_POSITIVE_SIDE);
  assert(ifsos(fo_hn, tab_h+2, tab_h+4, tab_h[1])==CGAL::ON_NEGATIVE_SIDE);
  assert(ifsos(fo_hp, tab_h+2, tab_h+4, tab_h[4])==CGAL::ON_NEGATIVE_SIDE);
  assert(ifsos(fo_hn, tab_h+2, tab_h+4, tab_h[4])==CGAL::ON_POSITIVE_SIDE);
  P tab_v[]={P(42,0),P(42,1),P(42,4),P(42,2),P(42,3)};
  assert(ar(tab_v+0,tab_v+5)==1);
  // FIXME: Triangulation says cah is only for independent range, but not Kernel_d
  // assert(cah(tab_v+0,tab_v+4,tab_v[4]));
  assert(cah(tab_v+0,tab_v+2,tab_v[4]));
  FO fo_vp = cfo (tab_v+0, tab_v+2);
  FO fo_vn = cfo (tab_v+2, tab_v+4);
  assert(ifo(fo_vp, tab_v+1, tab_v+3)==CGAL::POSITIVE);
  assert(ifo(fo_vn, tab_v+1, tab_v+3)==CGAL::NEGATIVE);
  assert(ifo(fo_vn, tab_v+2, tab_v+4)==CGAL::POSITIVE);
  assert(ifo(fo_vp, tab_v+2, tab_v+4)==CGAL::NEGATIVE);
  assert(ifsos(fo_vp, tab_v+1, tab_v+3, tab_v[3])==CGAL::ON_POSITIVE_SIDE);
  assert(ifsos(fo_vn, tab_v+1, tab_v+3, tab_v[3])==CGAL::ON_NEGATIVE_SIDE);
  assert(ifsos(fo_vp, tab_v+0, tab_v+2, tab_v[2])==CGAL::ON_NEGATIVE_SIDE);
  assert(ifsos(fo_vn, tab_v+0, tab_v+2, tab_v[2])==CGAL::ON_POSITIVE_SIDE);
  assert(ifsos(fo_vp, tab_v+2, tab_v+4, tab_v[1])==CGAL::ON_POSITIVE_SIDE);
  assert(ifsos(fo_vn, tab_v+2, tab_v+4, tab_v[1])==CGAL::ON_NEGATIVE_SIDE);
  assert(ifsos(fo_vp, tab_v+2, tab_v+4, tab_v[4])==CGAL::ON_NEGATIVE_SIDE);
  assert(ifsos(fo_vn, tab_v+2, tab_v+4, tab_v[4])==CGAL::ON_POSITIVE_SIDE);
#endif
  P z0=cp( 0+2,5-3);
  P z1=cp(-5+2,0-3);
  assert(abs(sd(z0,z1)-50)<.0001);
  assert(ed(z0,z0) && !ed(z0,z1));
  P z2=cp( 3+2,4-3);
  P tabz[]={z0,z1,z2};
  Sp sp = csp(tabz+0,tabz+3);
  P cent0=cos(sp);
  P cent1=cos(tabz+0,tabz+3);
  assert(abs(cent0[0]-2)<.0001);
  assert(abs(cent0[1]+3)<.0001);
  assert(abs(cent1[0]-2)<.0001);
  assert(abs(cent1[1]+3)<.0001);
  assert(abs(sp.squared_radius()-25)<.0001);
  P x2py1 = tp(x2,y1);
  assert(x2py1[1]==-2);
  WP tw[]={cwp(cp(5,0),1.5),cwp(cp(2,std::sqrt(3)),1),cwp(cp(2,-std::sqrt(3)),1)};
  WP xw=pc(tw+0,tw+3);
  assert(abs(pod(xw,tw[0]))<.0001);
  assert(abs(pod(xw,tw[1]))<.0001);
  assert(abs(pod(xw,tw[2]))<.0001);
  assert(pdw(xw)[0]<2.95);
  assert(pdw(xw)[0]>2.5);
  assert(pw(xw)<2.95);
  assert(pw(xw)>2.5);
  assert(psbps(tw+0,tw+3,cwp(cp(5,0),1.499)) == CGAL::ON_UNBOUNDED_SIDE);
  assert(psbps(tw+0,tw+3,cwp(cp(5,0),1.500)) == CGAL::ON_BOUNDARY);
  assert(psbps(tw+0,tw+3,cwp(cp(5,0),1.501)) == CGAL::ON_BOUNDED_SIDE);
  WP tw2[]={cwp(cp(-3,2),1),cwp(cp(5,2),2),cwp(cp(1,6),1)};
  assert(psbps(tw2+0,tw2+2,tw2[2]) == CGAL::ON_UNBOUNDED_SIDE);
  assert(abs(csrsos(tw2+0,tw2+2)-14.5039)<.0001);
  assert(abs(csrsos(tw2+0,tw2+1)+1)<.0001);

  P tl=cp(2,5);
  P br=cp(4,-1);
  IB ib=cib(tl,br);
  P bl=cmv(ib);
  P tr=cMv(ib);
  assert(cc(bl,0)==2);
  assert(cc(bl,1)==-1);
  assert(cc(tr,0)==4);
  assert(cc(tr,1)==5);

  Sp un1; CGAL_USE(un1);
  H un2; CGAL_USE(un2);
  S un3; CGAL_USE(un3);
  P un4; CGAL_USE(un4);
  V un5; CGAL_USE(un5);
  CI un6; CGAL_USE(un6);
  FO un7; CGAL_USE(un7);
  L un8; CGAL_USE(un8);
  R un9; CGAL_USE(un9);
  D un10; CGAL_USE(un10);
}

// Fails for an exact kernel, so I split it here
template<class Ker>
void test2i(){
  typedef Ker K1;
  typedef typename K1::Point_d P;
  typedef typename K1::Sphere_d Sp;
  typedef typename K1::Point_of_sphere_d PS;
  typedef typename K1::Construct_point_d CP;
  typedef typename K1::Construct_sphere_d CSp;
  typedef typename K1::Equal_d E;
  typedef typename K1::Squared_distance_d SD;
  typedef typename K1::Center_of_sphere_d COS;
  Ker k
#if 0
    (2)
#endif
    ;
  CP cp Kinit(construct_point_d_object);
  PS ps Kinit(point_of_sphere_d_object);
  CSp csp Kinit(construct_sphere_d_object);
  E ed Kinit(equal_d_object);
  SD sd Kinit(squared_distance_d_object);
  COS cos Kinit(center_of_sphere_d_object);
  P z0=cp( 0+2,5-3);
  P z1=cp(-5+2,0-3);
  P z2=cp( 3+2,4-3);
  P tabz[]={z0,z1,z2};
  Sp sp = csp(tabz+0,tabz+3);
  P cent0=cos(sp);
  P psp0=ps(sp,0);
  P psp1=ps(sp,1);
  P psp2=ps(sp,2);
  assert(!ed(psp0,psp1));
  assert(!ed(psp0,psp2));
  assert(!ed(psp2,psp1));
  using std::abs;
  assert(abs(sd(cent0,psp0)-25)<.0001);
  assert(abs(sd(cent0,psp1)-25)<.0001);
  assert(abs(sd(cent0,psp2)-25)<.0001);
}

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable: 4512)
#endif

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

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

template<class Ker>
void test3(){
  typedef Ker K1;
  //typedef typename K1::FT FT;
  typedef typename K1::Point_d P;
  typedef typename K1::Cartesian_const_iterator_d CI;
  typedef typename K1::Vector_d V;
  typedef typename K1::Segment_d S;
  typedef typename K1::Flat_orientation_d FO;

  //typedef K1::Construct_point CP;
  typedef typename K1::Construct_point_d CP_;
  typedef typename K1::Construct_vector_d CV_;
  typedef typename K1::Construct_segment_d CS;
  typedef typename CGAL::Get_functor<K1, CGAL::Segment_extremity_tag>::type CSE;
  typedef typename K1::Construct_cartesian_const_iterator_d CCI;
  typedef typename K1::Orientation_d PO;
  typedef typename K1::Linearly_independent_d LI;
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
  typedef typename K1::Point_dimension_d PD;
  typedef typename K1::Affinely_independent_d AI;
  typedef typename K1::Scaled_vector_d SV;
  typedef typename K1::Side_of_bounded_sphere_d SBDS;

  Ker k
#if 1
    (3)
#endif
    ;
  CP_ cp_ Kinit(construct_point_d_object);
  CV_ cv_ Kinit(construct_vector_d_object);
  typename boost::mpl::if_<boost::is_same<typename Ker::Default_ambient_dimension,CGAL::Dynamic_dimension_tag>,Construct_point3_helper<CP_>,CP_>::type cp(cp_);
  typename boost::mpl::if_<boost::is_same<typename Ker::Default_ambient_dimension,CGAL::Dynamic_dimension_tag>,Construct_point3_helper<CV_>,CV_>::type cv(cv_);
  CCI ci Kinit(construct_cartesian_const_iterator_d_object);
  CC cc Kinit(compute_coordinate_d_object);
  CL cl Kinit(compare_lexicographically_d_object);
  PO po Kinit(orientation_d_object);
  CS cs Kinit(construct_segment_d_object);
  CSE cse (k);
  SV sv Kinit(scaled_vector_d_object);
  LI li Kinit(linearly_independent_d_object);
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
  PD pd Kinit(point_dimension_d_object);
  AI ai Kinit(affinely_independent_d_object);
  SBDS sbds Kinit(side_of_bounded_sphere_d_object);
  using std::abs;
  P a; // Triangulation needs this :-(
  a=cp(2,3,4);
  assert(pd(a)==3);
  P b=cp(5,6,7,8);
  assert(ed(a,a));
  assert(!ed(a,b));
  assert(ca.dimension(a)==3);
  assert(ca.cartesian(a,1)==3);
  assert(ca.homogeneous(a,1)==3);
  assert(ca.homogeneous(a,3)==1);
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
  assert(abs(sd(e,a)-32)<.0001);
  P tab[]={a,b,c,d,e};
  std::cout << po (&tab[0],tab+4) << ' ';
  std::cout << sos(&tab[1],tab+5,tab[0]) << ' ';
  std::cout << sbs(tab+1,tab+5,tab[0]) << std::endl;
  FO fo=cfo(&tab[0],tab+3);
  std::cout << fo;
  P x[]={cp(2,2,3),cp(2,2,0),cp(1,2,1)};
  FO fo2=cfo(&x[0],x+3);
  std::cout << fo2;
  P y[]={cp(0,2,4),cp(3,1,2),cp(3,3,6),cp(0,4,8)};
  assert(!cis(x+0,x+3,y[0]));
  V yv[]={cv(0,2,4),cv(3,1,2),cv(3,3,6),cv(0,4,8)};
  assert( clh(yv+0,yv+1,yv[3]));
  assert( clh(yv+0,yv+2,yv[2]));
  assert(!clh(yv+0,yv+1,yv[2]));
  assert( li(yv+0,yv+2));
  assert(!li(yv+0,yv+3));
  assert( cah(y+0,y+3,y[3]));
  assert(!cah(y+0,y+2,y[2]));
  assert( ai(y+0,y+3));
  assert(!ai(y+0,y+4));
  assert(sv(yv[0],3)[1]==6);
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
  std::cout << sos(buf+1,buf+5,buf[0]) << ' ';
  buf[1]=y[1];buf[2]=y[2];buf[3]=y[3];buf[4]=y[0];
  std::cout << sos(buf+1,buf+5,buf[0]) << ' ';
  buf[1]=yy[0];buf[2]=yy[1];buf[3]=yy[2];buf[4]=yy[3];
  std::cout << sos(buf+1,buf+5,buf[0]) << ' ';
  buf[1]=yy[1];buf[2]=yy[2];buf[3]=yy[3];buf[4]=yy[0];
  std::cout << sos(buf+1,buf+5,buf[0]) << '\n';
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
  assert(!ed(x1,x2));
  P tab2[]={x1,x2,x3,x4,x5};
  assert(cis(tab2+0,tab2+4,x5));
  assert(po(tab2+0,tab2+4)==CGAL::POSITIVE);
  assert(sos(tab2+0,tab2+4,x5)==CGAL::ON_POSITIVE_SIDE);
  assert(sbs(tab2+0,tab2+4,x5)==CGAL::ON_BOUNDED_SIDE);
  FO fo4=cfo(tab2+0,tab2+3);
  assert(ifo(fo4,tab2+0,tab2+3)==CGAL::POSITIVE);
  assert(ifsos(fo4,tab2+0,tab2+3,x6)==CGAL::ON_POSITIVE_SIDE);
  P tx[]={cp(1,1,42),cp(3,1,42),cp(1,3,42),cp(3,3,42),cp(2,2,42),cp(4,5,42)};
  FO foxp=cfo(tx+0,tx+3);
  FO foxn=cfo(tx+1,tx+4);
  assert(ifo(foxp, tx+0, tx+3) == CGAL::POSITIVE);
  assert(ifo(foxn, tx+0, tx+3) == CGAL::NEGATIVE);
  assert(ifo(foxp, tx+1, tx+4) == CGAL::NEGATIVE);
  assert(ifo(foxn, tx+1, tx+4) == CGAL::POSITIVE);
  assert(ifsos(foxp, tx+0, tx+3, tx[3]) == CGAL::ON_ORIENTED_BOUNDARY);
  assert(ifsos(foxn, tx+0, tx+3, tx[3]) == CGAL::ON_ORIENTED_BOUNDARY);
  assert(ifsos(foxp, tx+0, tx+3, tx[4]) == CGAL::ON_POSITIVE_SIDE);
  assert(ifsos(foxn, tx+0, tx+3, tx[4]) == CGAL::ON_NEGATIVE_SIDE);
  assert(ifsos(foxp, tx+0, tx+3, tx[5]) == CGAL::ON_NEGATIVE_SIDE);
  assert(ifsos(foxn, tx+0, tx+3, tx[5]) == CGAL::ON_POSITIVE_SIDE);
  P ty[]={cp(1,42,1),cp(3,42,1),cp(1,42,3),cp(3,42,3),cp(2,42,2),cp(4,42,5)};
  FO foyp=cfo(ty+0,ty+3);
  FO foyn=cfo(ty+1,ty+4);
  assert(ifo(foyp, ty+0, ty+3) == CGAL::POSITIVE);
  assert(ifo(foyn, ty+0, ty+3) == CGAL::NEGATIVE);
  assert(ifo(foyp, ty+1, ty+4) == CGAL::NEGATIVE);
  assert(ifo(foyn, ty+1, ty+4) == CGAL::POSITIVE);
  assert(ifsos(foyp, ty+0, ty+3, ty[3]) == CGAL::ON_ORIENTED_BOUNDARY);
  assert(ifsos(foyn, ty+0, ty+3, ty[3]) == CGAL::ON_ORIENTED_BOUNDARY);
  assert(ifsos(foyp, ty+0, ty+3, ty[4]) == CGAL::ON_POSITIVE_SIDE);
  assert(ifsos(foyn, ty+0, ty+3, ty[4]) == CGAL::ON_NEGATIVE_SIDE);
  assert(ifsos(foyp, ty+0, ty+3, ty[5]) == CGAL::ON_NEGATIVE_SIDE);
  assert(ifsos(foyn, ty+0, ty+3, ty[5]) == CGAL::ON_POSITIVE_SIDE);
  P tz[]={cp(1,1,42),cp(3,1,42),cp(1,3,42),cp(3,3,42),cp(2,2,42),cp(4,5,42)};
  FO fozp=cfo(tz+0,tz+3);
  FO fozn=cfo(tz+1,tz+4);
  assert(ifo(fozp, tz+0, tz+3) == CGAL::POSITIVE);
  assert(ifo(fozn, tz+0, tz+3) == CGAL::NEGATIVE);
  assert(ifo(fozp, tz+1, tz+4) == CGAL::NEGATIVE);
  assert(ifo(fozn, tz+1, tz+4) == CGAL::POSITIVE);
  assert(ifsos(fozp, tz+0, tz+3, tz[3]) == CGAL::ON_ORIENTED_BOUNDARY);
  assert(ifsos(fozn, tz+0, tz+3, tz[3]) == CGAL::ON_ORIENTED_BOUNDARY);
  assert(ifsos(fozp, tz+0, tz+3, tz[4]) == CGAL::ON_POSITIVE_SIDE);
  assert(ifsos(fozn, tz+0, tz+3, tz[4]) == CGAL::ON_NEGATIVE_SIDE);
  assert(ifsos(fozp, tz+0, tz+3, tz[5]) == CGAL::ON_NEGATIVE_SIDE);
  assert(ifsos(fozn, tz+0, tz+3, tz[5]) == CGAL::ON_POSITIVE_SIDE);
  P showit=cp(1,2,4);
  std::ostringstream output;
  output << showit;
  assert(output.str()=="3 1 2 4");
  std::istringstream input("3 5 6 9");
  input >> showit;
  assert(ed(showit,cp(5,6,9)));
  P t1[]={cp(1,2,3),cp(3,2,1),cp(2,4,2)};
  assert(sbds(t1+0,t1+2,cp(2,2,3.414)) == CGAL::ON_BOUNDED_SIDE);
  assert(sbds(t1+0,t1+2,cp(1,2,3)) == CGAL::ON_BOUNDARY);
  assert(sbds(t1+0,t1+2,cp(2,2,3.415)) == CGAL::ON_UNBOUNDED_SIDE);
  assert(sbds(t1+0,t1+3,cp(2.1,3.5,1.9)) == CGAL::ON_BOUNDED_SIDE);
  assert(sbds(t1+0,t1+3,cp(10,10,10)) == CGAL::ON_UNBOUNDED_SIDE);

  typedef typename K1::Weighted_point_d WP;
  typedef typename K1::Construct_weighted_point_d CWP;
  //typedef typename K1::Point_drop_weight_d PDW;
  typedef CP_ PDW;
  typedef typename K1::Compute_weight_d PW;
  typedef typename K1::Power_side_of_power_sphere_d PT;
  typedef typename K1::In_flat_power_side_of_power_sphere_d IFPT;
  CWP cwp Kinit(construct_weighted_point_d_object);
  //PDW pdw Kinit(point_drop_weight_d_object);
  PDW const& pdw = cp_;
  PW pw Kinit(compute_weight_d_object);
  PT pt Kinit(power_side_of_power_sphere_d_object);
  IFPT ifpt Kinit(in_flat_power_side_of_power_sphere_d_object);
  WP wp;
  wp = cwp (x1, 2);
  WP xw6 = cwp (x6, 0);
  assert (pw(wp) == 2);
  assert (ed(pdw(wp), x1));
  WP tabw[]={cwp(x1,0),cwp(x2,0),cwp(x3,0),cwp(x4,0),cwp(x5,0)};
  assert(pt(tabw+0,tabw+4,tabw[4])==CGAL::ON_POSITIVE_SIDE);
  assert(ifpt(fo4,tabw+0,tabw+3,xw6)==CGAL::ON_POSITIVE_SIDE);
  std::ostringstream swp; swp << wp; assert(swp.str()=="3 0 1 -1 2");
  std::istringstream swp2("3 4 5 6 7");
  swp2 >> wp;
  assert(ed(wp.point(),cp(4,5,6)) && wp.weight()==7);

  V v1=cv(3,2,1);
  std::ostringstream sv1; sv1 << v1; assert(sv1.str()=="3 3 2 1");
  std::istringstream sv2("3 4 5 6"); sv2 >> v1; assert(v1[0]==4&&v1[1]==5);
}
template<class Ker>
void test4(){
  typedef typename Ker::Point_d P;
  typedef typename Ker::Weighted_point_d WP;
  typedef typename Ker::Construct_circumcenter_d CCc;
  typedef typename Ker::Equal_d E;
  typedef typename Ker::Construct_power_sphere_d PC;
  typedef typename Ker::Compute_power_product_d PoD;
  typedef typename Ker::Affine_rank_d AR;
  Ker k(4);
  CCc ccc Kinit(construct_circumcenter_d_object);
  E ed Kinit(equal_d_object);
  PC pc Kinit(construct_power_sphere_d_object);
  PoD pod Kinit(compute_power_product_d_object);
  AR ar Kinit(affine_rank_d_object);
  auto mkpt=[](auto...x){double l[]{(double)x...};return P(std::begin(l), std::end(l));};
  P tab1[]={mkpt(15,20,40,80),mkpt(10,23,36,80),mkpt(10,20,40,85),mkpt(10,15,40,80),mkpt(13,20,40,76)};
  assert(ed(ccc(tab1+0, tab1+5),mkpt(10,20,40,80)));
  P tab2[]={mkpt(15,20,40,80),mkpt(13,24,40,80),mkpt(10,25,40,80),mkpt(10,20,43,84)};
  assert(ed(ccc(tab2+0, tab2+4),mkpt(10,20,40,80)));
  P tab3[]={mkpt(15,20,35,80),mkpt(10,25,40,75),mkpt(13,24,37,76)};
  assert(ed(ccc(tab3+0, tab3+3),mkpt(10,20,40,80)));
  auto mkwpt=[](auto...x){double l[]{(double)x...};auto last=std::prev(std::end(l));return WP(P(std::begin(l), last),*last);};
  WP tab4[]={mkwpt(89,17,29,97,14),mkwpt(86,99,64,26,44),mkwpt(40,9,13,91,20),mkwpt(41,30,93,13,10),mkwpt(45,6,98,9,0),mkwpt(0,0,0,0,0)};
  for(int i=5;i>=1;--i){
    tab4[i]=pc(tab4+0, tab4+i);
    for(int j=0;j<i;++j)
      assert(pod(tab4[i],tab4[j])==0);
    auto drop=[](WP const&x){return x.point();};
    assert(ar(CGAL::make_transforming_iterator(tab4+0,drop), CGAL::make_transforming_iterator(tab4+i+1,drop))==i-1);
  }
}
template struct CGAL::Epick_d<CGAL::Dimension_tag<2> >;
template struct CGAL::Epick_d<CGAL::Dimension_tag<3> >;
template struct CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
typedef CGAL::Epick_d<CGAL::Dimension_tag<2> > Ker2;
typedef CGAL::Epick_d<CGAL::Dimension_tag<3> > Ker3;
typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kerd;
CGAL_static_assertion((boost::is_same<CGAL::Dimension_tag<2>,Ker2::Dimension>::value));
CGAL_static_assertion((boost::is_same<CGAL::Dimension_tag<3>,Ker3::Dimension>::value));
CGAL_static_assertion((boost::is_same<CGAL::Dynamic_dimension_tag,Kerd::Dimension>::value));
CGAL_static_assertion((boost::is_same<CGAL::Dimension_tag<2>,CGAL::Ambient_dimension<Ker2::Point_d>::type>::value));
CGAL_static_assertion((boost::is_same<CGAL::Dimension_tag<3>,CGAL::Ambient_dimension<Ker3::Point_d,Ker3>::type>::value));
int main(){
  //Broken with Linear_base_d (output iterator)
  //test2<CGAL::Kernel_d_interface<KK> >();
  test2<Ker2>(); test2i<Ker2>();
  test3<Ker3>();
  test3<Kerd>();
#if !defined _MSC_VER || _MSC_VER >= 1910
  test2<CGAL::Epeck_d<CGAL::Dimension_tag<2>>>();
  test3<CGAL::Epeck_d<CGAL::Dimension_tag<3>>>();
  test3<CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>>();
  test4<CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>>();
#endif
}

#endif
