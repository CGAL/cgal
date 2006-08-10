#ifdef CGAL_NEF3_USE_LEDA_INTEGER
#include <CGAL/leda_integer.h>
typedef leda_integer NT;
#endif

#ifdef CGAL_NEF3_USE_GMPZ
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz NT;
#endif

#ifdef CGAL_NEF3_USE_GMPQ
#include <CGAL/Gmpq.h>
typedef CGAL::Gmpq NT;
#endif

#ifdef CGAL_NEF3_USE_LAZY_EXACT_NT
#include <CGAL/Lazy_nt.h>
typedef CGAL::Lazy_nt<NT> RT;
#else
typedef NT RT;
#endif

#ifdef CGAL_NEF3_USE_SIMPLE_HOMOGENEOUS
#include <CGAL/Simple_homogeneous.h>
typedef CGAL::Simple_homogeneous<RT> Kernel;
#endif

#ifdef CGAL_NEF3_USE_HOMOGENEOUS
#include <CGAL/Homogeneous.h>
typedef CGAL::Homogeneous<RT> Kernel;
#endif

#ifdef CGAL_NEF3_USE_EXTENDED_HOMOGENEOUS
#include <CGAL/Extended_homogeneous.h>
typedef CGAL::Extended_homogeneous<RT> Kernel;
#endif

#ifdef CGAL_NEF3_USE_LEDA_RATIONAL
#include <CGAL/leda_rational.h>
typedef leda_rational NT;
#endif

#ifdef CGAL_NEF3_USE_GMPQ
#include <CGAL/Gmpq.h>
typedef CGAL::Gmpq NT;
#endif

#ifdef CGAL_NEF3_USE_SIMPLE_CARTESIAN
#include <CGAL/Simple_cartesian.h>
typedef CGAL::Simple_cartesian<NT> Kernel;
#endif

#ifdef CGAL_NEF3_USE_CARTESIAN
#include <CGAL/Cartesian.h>
typedef CGAL::Cartesian<NT> Kernel;
#endif

#ifdef CGAL_NEF3_USE_EXTENDED_CARTESIAN
#include <CGAL/Extended_cartesian.h>
typdef CGAL::Extended_cartesian<NT> Kernel;
#endif

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Timer.h>
#include <sstream>
#include <fstream>
#include "tetrahedron_generator.h"
#include "grid_generator.h"

typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Point_3 Point_3;
typedef Nef_polyhedron::Vector_3 Vector_3;
typedef Nef_polyhedron::Aff_transformation_3 Aff_transformation_3;
typedef tetrahedron_generator<Kernel> tgen;
typedef CGAL::grid_generator<Nef_polyhedron> ggen;
typedef CGAL::Random_points_in_cube_3<Point_3> Point_source;

bool cgal_nef3_timer_on = false;

int main(int argc, char* argv[]) {

  assert(argc < 8);
  
  int nx = argc>1 ? std::atoi(argv[1]) : 2;
  int ny = argc>2 ? std::atoi(argv[2]) : 2;
  int nz = argc>3 ? std::atoi(argv[3]) : 2;
  int n  = argc>4 ? std::atoi(argv[4]) : 2;
  int s  = argc>6 ? std::atoi(argv[6]) : 100;

  std::ostringstream out2;
  if(argc>5) {
    tgen t2(out2,s,std::atoi(argv[5]));
    t2.create_tetrahedra(nx,ny,nz);
  } else {
    tgen t2(out2,s);
    t2.create_tetrahedra(nx,ny,nz);    
  }
  std::istringstream in2(out2.str());
  Nef_polyhedron N2;
  in2 >> N2;
  CGAL_assertion(N2.is_valid());

  RT b=s*nx;
  N2.transform(Aff_transformation_3(CGAL::TRANSLATION,Vector_3(-b,-b,-b,2)));

  CGAL::Timer pl;
  pl.start();

  Point_source P(CGAL::to_double(b)/2);
  for(int i=0;i<n;++i) {
    N2.locate(*P++);
  }
  pl.stop();

  std::cout << "Input_size: " << N2.number_of_vertices() << std::endl;
  std::cout << "Number_of_point_location_queries: " << n << std::endl;
  std::cout << "Total_runtime: " << pl.time() << std::endl;
  std::cout << "Runtime_per_query: " << pl.time()/n << std::endl;
}
