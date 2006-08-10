#ifdef CGAL_NEF3_USE_LEDA_INTEGER
#include <CGAL/leda_integer.h>
typedef leda_integer NT;
#endif

#ifdef CGAL_NEF3_USE_GMPZ
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz NT;
#endif

#ifdef CGAL_NEF3_USE_LAZY_EXACT_NT
#include <CGAL/Lazy_nt>
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

#include <CGAL/basic.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include "worst_case_sphere_map_generator.h"

typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Vector_3 Vector_3;
typedef Nef_polyhedron::Aff_transformation_3 Aff_transformation_3;
typedef Nef_polyhedron::Vertex_const_iterator Vertex_const_iterator;

bool cgal_nef3_timer_on = false;

int main(int argc, char* argv[]) {

  CGAL_assertion(argc==2);

  Nef_polyhedron N1,N2;
  std::ostringstream out;
  worst_case_sphere_map_generator<Nef_polyhedron> wgen(out);
  wgen.print(std::atoi(argv[1]));
  std::istringstream in(out.str());
  in >> N1;
  CGAL_assertion(N1.is_valid());

  N2=N1;
  N2.transform(Aff_transformation_3(1,0,0, 0,0,-1, 0,1,0, 1));

  cgal_nef3_timer_on = true;
  N1=N2.join(N1);
};
