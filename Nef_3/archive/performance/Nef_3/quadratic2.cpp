#ifdef CGAL_NEF3_USE_LEDA_INTEGER
#include <CGAL/leda_integer.h>
typedef leda_integer NT;
#endif

#ifdef CGAL_NEF3_USE_GMPZ
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz NT;
#endif

#ifdef CGAL_NEF3_USE_GMPZ
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz NT;
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

#ifdef CGAL_NEF3_USE_FILTERED_HOMOGENEOUS_3
#include <CGAL/Filtered_homogeneous_3.h>
typedef CGAL::Filtered_homogeneous_3<RT> Kernel;
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
#include "grid_generator.h"

typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Vector_3 Vector_3;
typedef Nef_polyhedron::Aff_transformation_3 Aff_transformation_3;
typedef CGAL::grid_generator<Nef_polyhedron> ggen;

bool cgal_nef3_timer_on = false;

int main(int argc, char* argv[]) {

  assert(argc>1 && argc < 5);

  int nx = argc>2 ? std::atoi(argv[2]) : 2;
  int ny = argc>3 ? std::atoi(argv[3]) : nx;

  std::ifstream in(argv[1]);
  Nef_polyhedron Nin;
  in >> Nin;

  std::ostringstream out1;
  ggen g(out1, Nin,false);
  g.print(nx,1,1);
  std::istringstream in1(out1.str());
  Nef_polyhedron N1;
  in1 >> N1;
  CGAL_assertion(N1.is_valid());

  Nin.transform(Aff_transformation_3(0,-1,0,
                                     1,0,0,
                                     0,0,1,1));

  Nin.transform(Aff_transformation_3(CGAL::TRANSLATION, Vector_3(0,0,1)));

  std::ostringstream out2;
  ggen g2(out2, Nin,false);
  g2.print(1,ny,1);
  std::istringstream in2(out2.str());
  Nef_polyhedron N2;
  in2 >> N2;
  CGAL_assertion(N2.is_valid());

  cgal_nef3_timer_on = true;

#if defined CGAL_NEF3_UNION
  N1=N1.union(N2);
#elif defined CGAL_NEF3_INTERSECTION
  N1=N1.intersection(N2);
#elif defined CGAL_NEF3_DIFFERENCE
  N1=N1.difference(N2);
#else
  N1=N1.symmetric_difference(N2);
#endif
}
