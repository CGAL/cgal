#line 1433 "Simple_extended_kernel.lw"
#define RPOLYNOMIAL_EXPLICIT_OUTPUT
#include <CGAL/Cartesian.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Extended_cartesian.h>
#include <CGAL/Filtered_extended_homogeneous.h>
#include <CGAL/test_macros.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer Integer;
template <>
struct ring_or_field<leda_integer> {
  typedef ring_with_gcd kind;
  typedef leda_integer RT;
  static RT gcd(const RT& r1, const RT& r2) 
  { return ::gcd(r1,r2); }
};
#include <CGAL/leda_real.h>
typedef leda_real Real;
template <>
struct ring_or_field<leda_real> {
  typedef field_with_div kind;
};
#else
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#include <CGAL/Quotient.h>
typedef CGAL::Gmpz Integer;
template <>
struct ring_or_field<CGAL::Gmpz> {
  typedef ring_with_gcd kind;
  typedef CGAL::Gmpz RT;
  static RT gcd(const RT& r1, const RT& r2) 
  { return CGAL::gcd(r1,r2); }
};
typedef CGAL::Quotient<Integer> Real;
template <>
struct ring_or_field<Real> {
  typedef field_with_div kind;
};
#else
typedef long   Integer;
typedef double Real;
#endif
#endif

#include <CGAL/intersection_2.h>
using namespace CGAL;

int main()
{
  SETDTHREAD(41);
  CGAL_TEST_START;
{
  typedef CGAL::Extended_homogeneous<Integer> EDec;
  
#line 1507 "Simple_extended_kernel.lw"
typedef EDec::Point_2     EP;
typedef EDec::Segment_2   ES;
typedef EDec::Direction_2 ED;
typedef EDec::Standard_kernel::Point_2  Point;
typedef EDec::Standard_kernel::Line_2   Line;
typedef EDec::Standard_RT RT;

EDec D;
CGAL::set_pretty_mode ( std::cerr );
Point ps1(0,0), ps2(1,1), ps3(1,0), ps4(0,1), ps5(1,1,2);
EDec::Point_type t1,t2,t3;
EP eps1 = D.construct_point(ps1);
EP eps2 = D.construct_point(ps2); 
EP eps3 = D.construct_point(ps3);
EP eps4 = D.construct_point(ps4);
EP eps5 = D.construct_point(ps5);
EP epn1 = D.construct_point(ps4,ps1,t1); // vertical down ray
EP epn2 = D.construct_point(ps1,ps4,t2); // vertical up ray
EP epn3 = D.construct_point(ps1,ps3,t3); // horizontal right ray
EP epn4 = D.construct_point(Line(2,3,4));
ES el1 =  D.construct_segment(D.construct_point(Line(ps1,ps4)),
                              D.construct_point(Line(ps4,ps1)));
ES el2 =  D.construct_segment(D.NW(),D.NE());

CGAL_TEST(D.type(D.SW())==EDec::SWCORNER);
CGAL_TEST(D.type(D.NW())==EDec::NWCORNER);
CGAL_TEST(D.type(D.SE())==EDec::SECORNER);
CGAL_TEST(D.type(D.NE())==EDec::NECORNER);
CGAL_TEST(D.type(epn1)==EDec::BOTTOMFRAME);
CGAL_TEST(D.type(epn2)==EDec::TOPFRAME);
CGAL_TEST(D.type(epn3)==EDec::RIGHTFRAME);
CGAL_TEST(D.type(eps1)==EDec::STANDARD);

CGAL_TEST(D.standard_line(epn1) == Line(ps4,ps1));
CGAL_TEST(D.standard_point(eps1) == ps1);

ES es1 = D.construct_segment(epn1,epn2);
ES es2 = D.construct_segment(eps1,eps5);
CGAL_TEST(D.source(es1) == epn1);
CGAL_TEST(D.target(es1) == epn2);
CGAL_TEST(D.orientation(es1,D.construct_point(Point(-1,-2))) > 0 );
CGAL_TEST(D.is_degenerate(D.construct_segment(epn1,epn1)));
CGAL_TEST(D.compare_xy(eps1,eps5)<0);
CGAL_TEST(D.compare_xy(eps1,epn2)<0);
CGAL_TEST(D.compare_xy(epn1,eps1)<0);
CGAL_TEST(D.intersection(es1,es2) == eps1);

CGAL_TEST(D.compare_xy(eps1,eps2)<0);
CGAL_TEST(D.compare_xy(eps4,eps1)>0);
CGAL_TEST(D.compare_xy(eps1,eps1)==0);

CGAL_TEST(D.compare_xy(D.NW(),eps2)<0);
CGAL_TEST(D.compare_xy(eps1,D.NE())<0);
CGAL_TEST(D.compare_xy(D.SW(),D.SE())<0);
CGAL_TEST(D.compare_xy(epn1,eps1)<0);
CGAL_TEST(D.compare_xy(eps1,epn2)<0);

CGAL_TEST(D.orientation(eps1,eps2,eps3)<0);
CGAL_TEST(D.orientation(eps1,eps3,eps2)>0);
CGAL_TEST(D.orientation(eps1,eps2,D.construct_point(Point(2,2)))==0);

CGAL_TEST(D.orientation(eps1,eps2,D.construct_point(ps1,ps2))==0);
CGAL_TEST(D.orientation(eps1,eps2,epn3)<0);
CGAL_TEST(D.orientation(eps1,eps2,epn2)>0);

CGAL_TEST(D.orientation(D.NW(),D.NE(),eps1)<0);
CGAL_TEST(D.orientation(D.NE(),D.NW(),eps1)>0);
CGAL_TEST(D.orientation(D.SW(),D.NE(),eps1)==0);
CGAL_TEST(D.orientation(epn1,epn2,eps1)==0);
CGAL_TEST(D.orientation(epn1,epn2,eps4)==0);
CGAL_TEST(D.orientation(epn1,epn2,eps3)<0);
CGAL_TEST(D.orientation(epn2,epn1,eps3)>0);

CGAL_TEST(D.first_pair_closer_than_second(eps1,eps5,eps1,eps2));
CGAL_TEST(!D.first_pair_closer_than_second(eps1,eps3,eps1,eps4));
CGAL_TEST(D.first_pair_closer_than_second(eps1,eps3,eps2,
          D.construct_point(Point(2,2,1))));
CGAL_TEST(D.first_pair_closer_than_second(eps1,eps3,eps2,D.NW()));
CGAL_TEST(D.first_pair_closer_than_second(eps1,D.SE(),D.SW(),D.NE()));
CGAL_TEST(!D.first_pair_closer_than_second(eps1,D.SE(),eps5,D.NE()));
CGAL_TEST(D.first_pair_closer_than_second(eps5,D.NE(),eps1,D.SE()));
CGAL_TEST(!D.first_pair_closer_than_second(eps1,D.SE(),eps1,D.NE()));
CGAL_TEST(D.first_pair_closer_than_second(D.SE(),D.NE(),D.NE(),D.SW()));
CGAL_TEST(!D.first_pair_closer_than_second(D.SE(),D.NE(),D.NW(),D.SW()));

CGAL_TEST(D.construct_direction(D.NW(),D.NE())==ED(RT(1),RT(0)));
CGAL_TEST(D.construct_direction(D.NE(),D.NW())==ED(RT(-1),RT(0)));
CGAL_TEST(D.construct_direction(D.SW(),D.NE())==ED(RT(1),RT(1)));
CGAL_TEST(D.construct_direction(D.NW(),D.SE())==ED(RT(1),RT(-1)));
CGAL_TEST(D.construct_direction(D.NW(),D.SW())==ED(RT(0),RT(-1)));
CGAL_TEST(D.construct_direction(D.SW(),D.NW())==ED(RT(0),RT(1)));

CGAL_TEST(D.construct_direction(eps5,D.NE())==ED(RT(1),RT(1)));
CGAL_TEST(D.construct_direction(eps5,D.SW())==ED(RT(-1),RT(-1)));

ES upper = D.construct_segment(D.NW(),D.NE()); 
ES left  = D.construct_segment(D.SW(),D.NW()); 
EP ep_res = D.intersection(el1,upper);
CGAL_TEST(ep_res==epn2);
ep_res = D.intersection(left, upper);
CGAL_TEST(ep_res == D.NW());

#line 1489 "Simple_extended_kernel.lw"
  
#line 1610 "Simple_extended_kernel.lw"
CGAL_IO_TEST(eps3,eps1);
CGAL_IO_TEST(epn2,epn1);


#line 1490 "Simple_extended_kernel.lw"
}
{  
  typedef CGAL::Extended_cartesian<Real> EDec;
  
#line 1507 "Simple_extended_kernel.lw"
typedef EDec::Point_2     EP;
typedef EDec::Segment_2   ES;
typedef EDec::Direction_2 ED;
typedef EDec::Standard_kernel::Point_2  Point;
typedef EDec::Standard_kernel::Line_2   Line;
typedef EDec::Standard_RT RT;

EDec D;
CGAL::set_pretty_mode ( std::cerr );
Point ps1(0,0), ps2(1,1), ps3(1,0), ps4(0,1), ps5(1,1,2);
EDec::Point_type t1,t2,t3;
EP eps1 = D.construct_point(ps1);
EP eps2 = D.construct_point(ps2); 
EP eps3 = D.construct_point(ps3);
EP eps4 = D.construct_point(ps4);
EP eps5 = D.construct_point(ps5);
EP epn1 = D.construct_point(ps4,ps1,t1); // vertical down ray
EP epn2 = D.construct_point(ps1,ps4,t2); // vertical up ray
EP epn3 = D.construct_point(ps1,ps3,t3); // horizontal right ray
EP epn4 = D.construct_point(Line(2,3,4));
ES el1 =  D.construct_segment(D.construct_point(Line(ps1,ps4)),
                              D.construct_point(Line(ps4,ps1)));
ES el2 =  D.construct_segment(D.NW(),D.NE());

CGAL_TEST(D.type(D.SW())==EDec::SWCORNER);
CGAL_TEST(D.type(D.NW())==EDec::NWCORNER);
CGAL_TEST(D.type(D.SE())==EDec::SECORNER);
CGAL_TEST(D.type(D.NE())==EDec::NECORNER);
CGAL_TEST(D.type(epn1)==EDec::BOTTOMFRAME);
CGAL_TEST(D.type(epn2)==EDec::TOPFRAME);
CGAL_TEST(D.type(epn3)==EDec::RIGHTFRAME);
CGAL_TEST(D.type(eps1)==EDec::STANDARD);

CGAL_TEST(D.standard_line(epn1) == Line(ps4,ps1));
CGAL_TEST(D.standard_point(eps1) == ps1);

ES es1 = D.construct_segment(epn1,epn2);
ES es2 = D.construct_segment(eps1,eps5);
CGAL_TEST(D.source(es1) == epn1);
CGAL_TEST(D.target(es1) == epn2);
CGAL_TEST(D.orientation(es1,D.construct_point(Point(-1,-2))) > 0 );
CGAL_TEST(D.is_degenerate(D.construct_segment(epn1,epn1)));
CGAL_TEST(D.compare_xy(eps1,eps5)<0);
CGAL_TEST(D.compare_xy(eps1,epn2)<0);
CGAL_TEST(D.compare_xy(epn1,eps1)<0);
CGAL_TEST(D.intersection(es1,es2) == eps1);

CGAL_TEST(D.compare_xy(eps1,eps2)<0);
CGAL_TEST(D.compare_xy(eps4,eps1)>0);
CGAL_TEST(D.compare_xy(eps1,eps1)==0);

CGAL_TEST(D.compare_xy(D.NW(),eps2)<0);
CGAL_TEST(D.compare_xy(eps1,D.NE())<0);
CGAL_TEST(D.compare_xy(D.SW(),D.SE())<0);
CGAL_TEST(D.compare_xy(epn1,eps1)<0);
CGAL_TEST(D.compare_xy(eps1,epn2)<0);

CGAL_TEST(D.orientation(eps1,eps2,eps3)<0);
CGAL_TEST(D.orientation(eps1,eps3,eps2)>0);
CGAL_TEST(D.orientation(eps1,eps2,D.construct_point(Point(2,2)))==0);

CGAL_TEST(D.orientation(eps1,eps2,D.construct_point(ps1,ps2))==0);
CGAL_TEST(D.orientation(eps1,eps2,epn3)<0);
CGAL_TEST(D.orientation(eps1,eps2,epn2)>0);

CGAL_TEST(D.orientation(D.NW(),D.NE(),eps1)<0);
CGAL_TEST(D.orientation(D.NE(),D.NW(),eps1)>0);
CGAL_TEST(D.orientation(D.SW(),D.NE(),eps1)==0);
CGAL_TEST(D.orientation(epn1,epn2,eps1)==0);
CGAL_TEST(D.orientation(epn1,epn2,eps4)==0);
CGAL_TEST(D.orientation(epn1,epn2,eps3)<0);
CGAL_TEST(D.orientation(epn2,epn1,eps3)>0);

CGAL_TEST(D.first_pair_closer_than_second(eps1,eps5,eps1,eps2));
CGAL_TEST(!D.first_pair_closer_than_second(eps1,eps3,eps1,eps4));
CGAL_TEST(D.first_pair_closer_than_second(eps1,eps3,eps2,
          D.construct_point(Point(2,2,1))));
CGAL_TEST(D.first_pair_closer_than_second(eps1,eps3,eps2,D.NW()));
CGAL_TEST(D.first_pair_closer_than_second(eps1,D.SE(),D.SW(),D.NE()));
CGAL_TEST(!D.first_pair_closer_than_second(eps1,D.SE(),eps5,D.NE()));
CGAL_TEST(D.first_pair_closer_than_second(eps5,D.NE(),eps1,D.SE()));
CGAL_TEST(!D.first_pair_closer_than_second(eps1,D.SE(),eps1,D.NE()));
CGAL_TEST(D.first_pair_closer_than_second(D.SE(),D.NE(),D.NE(),D.SW()));
CGAL_TEST(!D.first_pair_closer_than_second(D.SE(),D.NE(),D.NW(),D.SW()));

CGAL_TEST(D.construct_direction(D.NW(),D.NE())==ED(RT(1),RT(0)));
CGAL_TEST(D.construct_direction(D.NE(),D.NW())==ED(RT(-1),RT(0)));
CGAL_TEST(D.construct_direction(D.SW(),D.NE())==ED(RT(1),RT(1)));
CGAL_TEST(D.construct_direction(D.NW(),D.SE())==ED(RT(1),RT(-1)));
CGAL_TEST(D.construct_direction(D.NW(),D.SW())==ED(RT(0),RT(-1)));
CGAL_TEST(D.construct_direction(D.SW(),D.NW())==ED(RT(0),RT(1)));

CGAL_TEST(D.construct_direction(eps5,D.NE())==ED(RT(1),RT(1)));
CGAL_TEST(D.construct_direction(eps5,D.SW())==ED(RT(-1),RT(-1)));

ES upper = D.construct_segment(D.NW(),D.NE()); 
ES left  = D.construct_segment(D.SW(),D.NW()); 
EP ep_res = D.intersection(el1,upper);
CGAL_TEST(ep_res==epn2);
ep_res = D.intersection(left, upper);
CGAL_TEST(ep_res == D.NW());

#line 1494 "Simple_extended_kernel.lw"
  //IO does not work for LEDA reals
}
{
  typedef CGAL::Filtered_extended_homogeneous<Integer> EDec;
  
#line 1507 "Simple_extended_kernel.lw"
typedef EDec::Point_2     EP;
typedef EDec::Segment_2   ES;
typedef EDec::Direction_2 ED;
typedef EDec::Standard_kernel::Point_2  Point;
typedef EDec::Standard_kernel::Line_2   Line;
typedef EDec::Standard_RT RT;

EDec D;
CGAL::set_pretty_mode ( std::cerr );
Point ps1(0,0), ps2(1,1), ps3(1,0), ps4(0,1), ps5(1,1,2);
EDec::Point_type t1,t2,t3;
EP eps1 = D.construct_point(ps1);
EP eps2 = D.construct_point(ps2); 
EP eps3 = D.construct_point(ps3);
EP eps4 = D.construct_point(ps4);
EP eps5 = D.construct_point(ps5);
EP epn1 = D.construct_point(ps4,ps1,t1); // vertical down ray
EP epn2 = D.construct_point(ps1,ps4,t2); // vertical up ray
EP epn3 = D.construct_point(ps1,ps3,t3); // horizontal right ray
EP epn4 = D.construct_point(Line(2,3,4));
ES el1 =  D.construct_segment(D.construct_point(Line(ps1,ps4)),
                              D.construct_point(Line(ps4,ps1)));
ES el2 =  D.construct_segment(D.NW(),D.NE());

CGAL_TEST(D.type(D.SW())==EDec::SWCORNER);
CGAL_TEST(D.type(D.NW())==EDec::NWCORNER);
CGAL_TEST(D.type(D.SE())==EDec::SECORNER);
CGAL_TEST(D.type(D.NE())==EDec::NECORNER);
CGAL_TEST(D.type(epn1)==EDec::BOTTOMFRAME);
CGAL_TEST(D.type(epn2)==EDec::TOPFRAME);
CGAL_TEST(D.type(epn3)==EDec::RIGHTFRAME);
CGAL_TEST(D.type(eps1)==EDec::STANDARD);

CGAL_TEST(D.standard_line(epn1) == Line(ps4,ps1));
CGAL_TEST(D.standard_point(eps1) == ps1);

ES es1 = D.construct_segment(epn1,epn2);
ES es2 = D.construct_segment(eps1,eps5);
CGAL_TEST(D.source(es1) == epn1);
CGAL_TEST(D.target(es1) == epn2);
CGAL_TEST(D.orientation(es1,D.construct_point(Point(-1,-2))) > 0 );
CGAL_TEST(D.is_degenerate(D.construct_segment(epn1,epn1)));
CGAL_TEST(D.compare_xy(eps1,eps5)<0);
CGAL_TEST(D.compare_xy(eps1,epn2)<0);
CGAL_TEST(D.compare_xy(epn1,eps1)<0);
CGAL_TEST(D.intersection(es1,es2) == eps1);

CGAL_TEST(D.compare_xy(eps1,eps2)<0);
CGAL_TEST(D.compare_xy(eps4,eps1)>0);
CGAL_TEST(D.compare_xy(eps1,eps1)==0);

CGAL_TEST(D.compare_xy(D.NW(),eps2)<0);
CGAL_TEST(D.compare_xy(eps1,D.NE())<0);
CGAL_TEST(D.compare_xy(D.SW(),D.SE())<0);
CGAL_TEST(D.compare_xy(epn1,eps1)<0);
CGAL_TEST(D.compare_xy(eps1,epn2)<0);

CGAL_TEST(D.orientation(eps1,eps2,eps3)<0);
CGAL_TEST(D.orientation(eps1,eps3,eps2)>0);
CGAL_TEST(D.orientation(eps1,eps2,D.construct_point(Point(2,2)))==0);

CGAL_TEST(D.orientation(eps1,eps2,D.construct_point(ps1,ps2))==0);
CGAL_TEST(D.orientation(eps1,eps2,epn3)<0);
CGAL_TEST(D.orientation(eps1,eps2,epn2)>0);

CGAL_TEST(D.orientation(D.NW(),D.NE(),eps1)<0);
CGAL_TEST(D.orientation(D.NE(),D.NW(),eps1)>0);
CGAL_TEST(D.orientation(D.SW(),D.NE(),eps1)==0);
CGAL_TEST(D.orientation(epn1,epn2,eps1)==0);
CGAL_TEST(D.orientation(epn1,epn2,eps4)==0);
CGAL_TEST(D.orientation(epn1,epn2,eps3)<0);
CGAL_TEST(D.orientation(epn2,epn1,eps3)>0);

CGAL_TEST(D.first_pair_closer_than_second(eps1,eps5,eps1,eps2));
CGAL_TEST(!D.first_pair_closer_than_second(eps1,eps3,eps1,eps4));
CGAL_TEST(D.first_pair_closer_than_second(eps1,eps3,eps2,
          D.construct_point(Point(2,2,1))));
CGAL_TEST(D.first_pair_closer_than_second(eps1,eps3,eps2,D.NW()));
CGAL_TEST(D.first_pair_closer_than_second(eps1,D.SE(),D.SW(),D.NE()));
CGAL_TEST(!D.first_pair_closer_than_second(eps1,D.SE(),eps5,D.NE()));
CGAL_TEST(D.first_pair_closer_than_second(eps5,D.NE(),eps1,D.SE()));
CGAL_TEST(!D.first_pair_closer_than_second(eps1,D.SE(),eps1,D.NE()));
CGAL_TEST(D.first_pair_closer_than_second(D.SE(),D.NE(),D.NE(),D.SW()));
CGAL_TEST(!D.first_pair_closer_than_second(D.SE(),D.NE(),D.NW(),D.SW()));

CGAL_TEST(D.construct_direction(D.NW(),D.NE())==ED(RT(1),RT(0)));
CGAL_TEST(D.construct_direction(D.NE(),D.NW())==ED(RT(-1),RT(0)));
CGAL_TEST(D.construct_direction(D.SW(),D.NE())==ED(RT(1),RT(1)));
CGAL_TEST(D.construct_direction(D.NW(),D.SE())==ED(RT(1),RT(-1)));
CGAL_TEST(D.construct_direction(D.NW(),D.SW())==ED(RT(0),RT(-1)));
CGAL_TEST(D.construct_direction(D.SW(),D.NW())==ED(RT(0),RT(1)));

CGAL_TEST(D.construct_direction(eps5,D.NE())==ED(RT(1),RT(1)));
CGAL_TEST(D.construct_direction(eps5,D.SW())==ED(RT(-1),RT(-1)));

ES upper = D.construct_segment(D.NW(),D.NE()); 
ES left  = D.construct_segment(D.SW(),D.NW()); 
EP ep_res = D.intersection(el1,upper);
CGAL_TEST(ep_res==epn2);
ep_res = D.intersection(left, upper);
CGAL_TEST(ep_res == D.NW());

#line 1499 "Simple_extended_kernel.lw"
  
#line 1610 "Simple_extended_kernel.lw"
CGAL_IO_TEST(eps3,eps1);
CGAL_IO_TEST(epn2,epn1);


#line 1500 "Simple_extended_kernel.lw"
  D.print_statistics();
}
  CGAL_TEST_END;
}


