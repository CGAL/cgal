#include <CGAL/basic.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <strstream>

#include <CGAL/Fixed_precision_nt.h>
#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpz.h>
#else
#  include <CGAL/MP_Float.h>
#endif
#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/point_generators_3.h>

typedef CGAL::Cartesian<CGAL::Fixed_precision_nt>::Point_2 Point2;
typedef CGAL::Cartesian<CGAL::Fixed_precision_nt>::Point_3 Point3;
#ifdef CGAL_USE_GMP
typedef CGAL::Cartesian<CGAL::Gmpz>::Point_2               Point2gmp;
typedef CGAL::Cartesian<CGAL::Gmpz>::Point_3               Point3gmp;
#else
typedef CGAL::Cartesian<CGAL::MP_Float>::Point_2               Point2gmp;
typedef CGAL::Cartesian<CGAL::MP_Float>::Point_3               Point3gmp;
#endif

Point2gmp convert(Point2 p)
{return Point2gmp(
    (int)(CGAL::to_double(p.x())/CGAL::Fixed_precision_nt::unit_value()),
    (int)(CGAL::to_double(p.y())/CGAL::Fixed_precision_nt::unit_value()));}

Point3gmp convert(Point3 p)
{return Point3gmp(
    (int)(CGAL::to_double(p.x())/CGAL::Fixed_precision_nt::unit_value()),
    (int)(CGAL::to_double(p.y())/CGAL::Fixed_precision_nt::unit_value()),
    (int)(CGAL::to_double(p.z())/CGAL::Fixed_precision_nt::unit_value()));}

template < class InputIterator >
bool test_orient2(InputIterator first, int N, int& p, int& z, int& n)
{
  bool result = true;  
  Point2    p1,p2,p3;
  Point2gmp P1,P2,P3;
  z=0,n=0,p=0;

  while(N--) {
    p1 = *first; P1 = convert(p1); ++first;
    p2 = *first; P2 = convert(p2); ++first;
    p3 = *first; P3 = convert(p3); ++first;
    CGAL::Orientation oF( CGAL::orientation(p1,p2,p3));
    CGAL::Orientation oG( CGAL::orientation(P1,P2,P3));
    if (oF != oG) {
      std::cout << "---------- Problem"<<p1<<p2<<p3<<(int)oF<<std::endl
	   << " diff from        "<<P1<<P2<<P3<<(int)oG<<std::endl;
      result = false;
    }
    switch(oF){
    case CGAL::POSITIVE : p++; break;
    case CGAL::ZERO     : z++; break;
    case CGAL::NEGATIVE : n++; break;
    }
  }
  std::cout << "answers : "<<p << " +;  "<<z<<" 0;  "<<n<<" -;"<<std::endl;
  return result;
}

template < class InputIterator >
bool test_orient3(InputIterator first, int N, int& p, int& z, int& n)
{
  bool result = true;
  Point3    p1,p2,p3,p4;
  Point3gmp P1,P2,P3,P4;
  z=0,n=0,p=0;

  while(N--) {
    p1 = *first; P1 = convert(p1); ++first;
    p2 = *first; P2 = convert(p2); ++first;
    p3 = *first; P3 = convert(p3); ++first;
    p4 = *first; P4 = convert(p4); ++first;
    CGAL::Orientation oF( CGAL::orientation(p1,p2,p3,p4));
    CGAL::Orientation oG( CGAL::orientation(P1,P2,P3,P4));
    if (oF != oG) {
      std::cout << "---------- Problem"<<p1<<p2<<p3<<p4<<(int)oF<<std::endl
	   << " diff from        "<<P1<<P2<<P3<<P4<<(int)oG<<std::endl;
      result = false;
    }
    switch(oF){
    case CGAL::POSITIVE : p++; break;
    case CGAL::ZERO     : z++; break;
    case CGAL::NEGATIVE : n++; break;
    }
  }
  std::cout << "answers : "<<p << " +;  "<<z<<" 0;  "<<n<<" -;"<<std::endl;
  return result;
}

template < class InputIterator >
bool test_sphere2(InputIterator first, int N, int& p, int& z, int& n)
{
  bool result = true;
  Point2    p1,p2,p3,p4;
  Point2gmp P1,P2,P3,P4;
  z=0,n=0,p=0;

  while(N--) {
    p1 = *first; P1 = convert(p1); ++first;
    p2 = *first; P2 = convert(p2); ++first;
    p3 = *first; P3 = convert(p3); ++first;
    p4 = *first; P4 = convert(p4); ++first;
    CGAL::Oriented_side oF( CGAL::side_of_oriented_circle(p1,p2,p3,p4));
    CGAL::Oriented_side oG( CGAL::side_of_oriented_circle(P1,P2,P3,P4));
    if ( oG
	 && (! CGAL::Fixed_precision_nt::is_perturbed_incircle()) 
	 && (oF != oG) ) {
      std::cout << "---------- Problem"<<p1<<p2<<p3<<p4<<(int)oF<<std::endl
	   << " diff from        "<<P1<<P2<<P3<<P4<<(int)oG<<std::endl;
      result = false;
    }
    switch(oF){
    case CGAL::POSITIVE : p++; break;
    case CGAL::ZERO     : z++; break;
    case CGAL::NEGATIVE : n++; break;
    }
  }
  std::cout << "answers : "<<p << " +;  "<<z<<" 0;  "<<n<<" -;"<<std::endl;
  return result;
}

template < class InputIterator >
bool test_sphere3(InputIterator first, int N, int& p, int& z, int& n)
{
  bool result = true;
  Point3    p1,p2,p3,p4,p5;
  Point3gmp P1,P2,P3,P4,P5;
  z=0,n=0,p=0;

  while(N--) {
    p1 = *first; P1 = convert(p1); ++first;
    p2 = *first; P2 = convert(p2); ++first;
    p3 = *first; P3 = convert(p3); ++first;
    p4 = *first; P4 = convert(p4); ++first;
    p5 = *first; P5 = convert(p5); ++first;
    CGAL::Oriented_side oF( CGAL::side_of_oriented_sphere(p1,p2,p3,p4,p5));
    CGAL::Oriented_side oG( CGAL::side_of_oriented_sphere(P1,P2,P3,P4,P5));
    if ( oG
	 && (! CGAL::Fixed_precision_nt::is_perturbed_insphere())
	 && (oF != oG)) {
      std::cout << "---------- Problem"<<p1<<p2<<p3<<p4<<p5<<(int)oF<<std::endl
	   << " diff from        "<<P1<<P2<<P3<<P4<<P5<<(int)oG<<std::endl;
      result = false;
    }
    switch(oF){
    case CGAL::POSITIVE : p++; break;
    case CGAL::ZERO     : z++; break;
    case CGAL::NEGATIVE : n++; break;
    }
  }
  std::cout << "answers : "<<p << " +;  "<<z<<" 0;  "<<n<<" -;"<<std::endl;
  return result;
}



int main()
{
  // The following is necessary on Intel platforms.
  CGAL::force_ieee_double_precision();

  std::cout <<"Test program for class CGAL::Fixed_precision_nt by"
       << " comparison to Gmpz or MP_Float"
       <<std::endl<<std::endl;

#ifdef FIXED_ARE_INTEGERS
    std::cout <<"Fixed_precision_nt are integers "<<std::endl;
  float MAX=16000000.0;
  // max is about 2^24, thus a fixed is an integer.
#else
    std::cout <<"Fixed_precision_nt are not integers"<<std::endl;
  float MAX=1020.0;
  //something else bigger than fixed sized data below
#endif
  
  bool test_result =  CGAL::Fixed_precision_nt::init(MAX);
  bool test_all =  true;
  int p,z,n;

#ifdef FIXED_ARE_INTEGERS
  test_result &= (CGAL::Fixed_precision_nt(12345.6) 
		  == CGAL::Fixed_precision_nt(12346.4) );
  std::cout << "Verifying rounding on Fixed_precision_nt     "<<std::endl;
  std::cout << ((test_result) ? " OK " : "        FAILED --------") << std::endl;
  test_all &= test_result; test_result = true;
#endif


  CGAL::Fixed_precision_nt f1(1),f2(2),f3(3),f4(4);
  //  CGAL::Fixed_precision_nt fi(2*MAX); // overflow causes warning
  //  in test suite, thus this line is removed

  std::cout << "Verifying basic arithemtic on Fixed_precision_nt     "<<std::endl;
  test_result &= (f1+f2==f3);
  test_result &= (f4/f2==f2);
  test_result &= (f4-f1==f3);
  test_result &= (f2*f2==f4);
  //  test_result &= ! CGAL::is_valid(fi);
  test_result &= CGAL::is_finite(f1);
  test_result &= ((-f2) != f2);
  test_result &= (f1 <= f2);
  test_result &= (f1 <  f2);
  test_result &= (f2 >= f1);
  test_result &= (f2 >  f1);
  f2 = f4;
  f2 += f2;
  f2 *= f2;
  f2 /= f1;
  f2 -= f4;
  std::cout <<"60="<< f2<<std::endl;

  std::cout << ((test_result) ? " OK " : "        FAILED --------") << std::endl;
  test_all &= test_result; test_result = true;

  CGAL::Random_points_in_square_2<Point2,
                                CGAL::Creator_uniform_2<double,Point2> >
    Rnd2 ( CGAL::Fixed_precision_nt::upper_bound() );

  CGAL::Random_points_in_cube_3<Point3,
                                CGAL::Creator_uniform_3<double,Point3> >
    Rnd3 ( CGAL::Fixed_precision_nt::upper_bound() );

  int N=1000;
  std::cout << "Verifying 2D orientation test on "<<N<<" random tests"<<std::endl;
  test_result &= test_orient2(Rnd2,N,p,z,n);
  std::cout << ((test_result) ? " OK " : "        FAILED --------") << std::endl;
  test_all &= test_result; test_result = true;


  const int colinN=6;
  // All these triples are collinear
  Point2 colin[colinN*3] = {
    Point2(32,56), Point2(32,56), Point2(32,56),    // 3 equal points
    Point2(32,56), Point2(32,56), Point2(32,0),     // 2 equal points
    Point2(16,28), Point2(32,56), Point2(64,112),   // collinear with 0,0
    Point2(16,32), Point2(32,60), Point2(64,116),   // translated
    Point2(32,16), Point2(60,32), Point2(116,64),   // swapped
    Point2(32*CGAL::Fixed_precision_nt::unit_value(),
	   16*CGAL::Fixed_precision_nt::unit_value()),
    Point2(60*CGAL::Fixed_precision_nt::unit_value(), 
	   32*CGAL::Fixed_precision_nt::unit_value()),
    Point2(116*CGAL::Fixed_precision_nt::unit_value(),
	   64.4*CGAL::Fixed_precision_nt::unit_value()),
                                                    // fixed rounding
  };  
  std::cout << "Verifying 2D orientation test on collinear points"<<std::endl;
  test_result &= test_orient2(colin,colinN,p,z,n);
  test_result &= (z==colinN);
  std::cout << ((test_result) ? " OK " : "        FAILED --------") << std::endl;
  test_all &= test_result; test_result = true;


  // All these triples are positive
  const int posN=4;
  Point2 pos[posN*3] = {
    Point2(0,0),   Point2(1,0),   Point2(0,1),      // base points
    Point2(-MAX,-MAX), Point2(MAX+1,MAX-1), Point2(-MAX+1,-MAX+1),
                                                    // almost collinear
    Point2(-MAX+1,-MAX+1), Point2(-MAX,-MAX), Point2(MAX+1,MAX-1), 
                                                    // permuted
    Point2(32*CGAL::Fixed_precision_nt::unit_value(),
	   16*CGAL::Fixed_precision_nt::unit_value()),
    Point2(60*CGAL::Fixed_precision_nt::unit_value(),
	   32*CGAL::Fixed_precision_nt::unit_value()),
    Point2(116*CGAL::Fixed_precision_nt::unit_value(),
	   64.6*CGAL::Fixed_precision_nt::unit_value()),
                                                    // fixed rounding
  };
  std::cout << "Verifying 2D orientation test on positive triples"<<std::endl;
  test_result &= test_orient2(pos,posN,p,z,n);
  test_result &= (p==posN);
  std::cout << ((test_result) ? " OK " : "        FAILED --------") << std::endl;
  test_all &= test_result; test_result = true;



  std::cout << "Verifying incircle test on "<<N<<" random tests"<<std::endl;
  test_result &= test_sphere2(Rnd2,N,p,z,n);
  std::cout << ((test_result) ? " OK " : "        FAILED --------") << std::endl;
  test_all &= test_result; test_result = true;


  const int cocircN = 20;
  const int cocirccolinN = 2;
  Point2 cocirc[cocircN*4] = {
  // All these points are on a circle of radius 325
   Point2(36,323),   Point2(80,315),   Point2(91,-312),  Point2(125,300), 
   Point2(165,280),  Point2(195,-260), Point2(204,253),  Point2(36,-323),
   Point2(80,-315),  Point2(91,-312),  Point2(125,-300), Point2(165,-280),
   Point2(195,-260), Point2(204,-253), Point2(-36,-323), Point2(-80,-315), 
   Point2(-91,-312), Point2(-125,-300),Point2(-165,-280),Point2(-195,-260),
   Point2(-204,-253),Point2(-36,323),  Point2(-80,315),  Point2(-91,312),
   Point2(-125,300), Point2(-165,280), Point2(-195,260), Point2(-204,253),
   Point2(80,315),   Point2(-80,-315), Point2(-80,315),  Point2(80,-315), 
   Point2(195,-260), Point2(204,253),  Point2(36,-323),  Point2(80,-315), 
   Point2(91,-312),  Point2(125,-300), Point2(165,-280), Point2(195,-260),
   Point2(204,-253), Point2(-36,-323), Point2(-80,-315), Point2(-91,-312),
   Point2(-125,-300),Point2(-165,-280),Point2(-195,-260),Point2(-204,-253), 
   Point2(-36,323),  Point2(-80,315),  Point2(-91,312),  Point2(-125,300), 
   Point2(-165,280), Point2(-195,260), Point2(-204,253), Point2(36,323),
   // All these points are on a circle and an hyperbola
   Point2(-1,-2),    Point2(-2,-1),    Point2(1,2),      Point2(2,1),
   Point2(-1,-899),  Point2(-899,-1),  Point2(1,899),    Point2(899,1),
   Point2(-510,1020),Point2(-1020,510),Point2(1020,-510),Point2(510,-1020),
   Point2(-977,979), Point2(-979,977), Point2(977,-979), Point2(979,-977),
   // these points are collinear
   Point2(500,-501), Point2(501,-171), Point2(502,159),  Point2(503,489),
   Point2(-977,-979),Point2(-378,-382),Point2(221,215),  Point2(820,812)
  };


  std::cout << "Verifying incircle test on cocircular points"<<std::endl;
  CGAL::Fixed_precision_nt::unperturb_incircle();
  test_result &= ! CGAL::Fixed_precision_nt::is_perturbed_incircle();
  test_result &= test_sphere2(cocirc,cocircN,p,z,n);
  test_result &= (z==cocircN);
  std::cout << ((test_result) ? " OK " : "        FAILED --------") << std::endl;
  test_all &= test_result; test_result = true;

  std::cout<<"Verifying perturbed incircle test on cocircular points"<<std::endl;
  CGAL::Fixed_precision_nt::perturb_incircle();
  test_result &= CGAL::Fixed_precision_nt::is_perturbed_incircle();
  (void) test_sphere2(cocirc,cocircN,p,z,n);
  test_result &= ((z==cocirccolinN)&&(n+z+p==cocircN));
  std::cout << ((test_result) ? " OK " : "        FAILED --------") << std::endl;
  test_all &= test_result; test_result = true;


  std::cout <<"Verifying 3D orientation test on "<<N<<" random tests"<<std::endl;
  test_result &= test_orient3(Rnd3,N,p,z,n);
  std::cout << "Verifying insphere test on "<<N<<" random tests"<<std::endl;
  test_result &= test_sphere3(Rnd3,N,p,z,n);
  CGAL::Fixed_precision_nt::perturb_insphere();
  std::cout << ((test_result) ? " OK " : "        FAILED --------") << std::endl;

  std::cout <<std::endl<<std::endl
       << ((test_all) ? "test succeeds" 
                         : "TEST FAILED   TEST FAILED   TEST FAILED")
       <<std::endl;
  return test_all ? 0 : -1;
}

