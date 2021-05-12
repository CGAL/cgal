#ifdef CGAL_NEF3_USE_LEDA_INTEGER
#include <CGAL/leda_integer.h>
typedef leda_integer NT;
#endif

#ifdef CGAL_NEF3_USE_LEDA_RATIONAL
#include <CGAL/leda_rational.h>
typedef leda_rational NT;
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
//#include <CGAL/double.h>
#include <CGAL/Filtered_exact.h>
//typedef CGAL::Filtered_exact<double,NT> RT;
#include <CGAL/Nef_3/Filtered_gcd.h>
#include <CGAL/Lazy_exact_nt.h>
typedef CGAL::Lazy_exact_nt<NT> RT;
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

#ifdef CGAL_NEF3_USE_SIMPLE_CARTESIAN
#include <CGAL/Simple_cartesian.h>
typedef CGAL::Simple_cartesian<RT> Kernel;
#endif

#ifdef CGAL_NEF3_USE_CARTESIAN
#include <CGAL/Cartesian.h>
typedef CGAL::Cartesian<RT> Kernel;
#endif

#ifdef CGAL_NEF3_USE_EXTENDED_CARTESIAN
#include <CGAL/Extended_cartesian.h>
typdef CGAL::Extended_cartesian<RT> Kernel;
#endif

#include <CGAL/basic.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/rational_rotation.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cmath>
#include <cstddef>

typedef CGAL::Nef_polyhedron_3<Kernel>      Nef_polyhedron;
typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;
typedef Kernel::Vector_3                    Vector_3;
typedef Kernel::Aff_transformation_3        Aff_transformation_3;

bool cgal_nef3_timer_on = false;

Aff_transformation_3 compute_transformation_matrix(double alpha) {

  double arc = CGAL_PI * alpha / 180.0;

  RT epsilon = 1;
  double sin_double = std::sin( arc);
  double cos_double = std::cos( arc);

  std::cout << "alpha: " << alpha << std::endl;
  std::cout << "arc: " << arc << std::endl;
  std::cout << "sin_double: " << sin_double << std::endl;
  std::cout << "cos_double: " << cos_double << std::endl;

  while(sin_double < 1000 || cos_double < 1000) {
    sin_double *= 10;
    cos_double *= 10;
    epsilon *= RT(10);
  }

  std::cout << "epsilon      : 1/" << epsilon << std::endl;

  RT sin_alpha(0);
  RT cos_alpha(0);
  RT w(0);

  CGAL::Timer t;
  t.start();
  CGAL::rational_rotation_approximation( arc,
                                         sin_alpha, cos_alpha, w,
                                         RT(1), RT(epsilon));
  t.stop();
  std::cout << "approx. time: " << t.time() << std::endl;

  Aff_transformation_3 aff( cos_alpha,-sin_alpha, RT(0),
                            sin_alpha, cos_alpha, RT(0),
                            RT(0), RT(0), w,
                            w);

  std::cout << "sin(alpha)*w: " << sin_alpha << std::endl;
  std::cout << "cos(alpha)*w: " << cos_alpha << std::endl;
  std::cout << "w: " << w << std::endl;

  return aff;
}

Aff_transformation_3 compute_transformation_matrix(RT sinus, RT cosinus, RT w) {

  double sin_double = CGAL::to_double(sinus) / CGAL::to_double(w);
  double arc = std::asin(sin_double);
  double alpha = arc * 180 / CGAL_PI;

  std::cout << "sin(alpha)*w: " << sinus << std::endl;
  std::cout << "cos(alpha)*w: " << cosinus << std::endl;
  std::cout << "w: " << w << std::endl;

  std::cout << "sin_double: " << sin_double << std::endl;
  std::cout << "arc: " << arc << std::endl;
  std::cout << "alpha: " << alpha << std::endl;

  Aff_transformation_3 aff( cosinus,-sinus, RT(0),
                            sinus, cosinus, RT(0),
                            RT(0), RT(0), w,
                            w);
  return aff;

}

int main(int argc, char* argv[]) {

  assert(argc == 3);

  std::ifstream in(argv[1]);
  CGAL_assertion_msg(in, "incorrect parameter 1");

  std::ifstream rotations(argv[2]);
  CGAL_assertion_msg(rotations, "incorrect parameter 2");

  std::string mode;
  rotations >> mode;
  CGAL_assertion(mode == "angle" || mode == "sinus");

  Polyhedron poly;
  in >> poly;
  Nef_polyhedron N1(poly);

  int runs;
  rotations >> runs;

  double alpha;
  RT sinus, cosinus, w;
  Aff_transformation_3 aff;
  for(int i=0; i<runs; i++) {

    if(mode=="angle") {
      rotations >> alpha;
      aff = compute_transformation_matrix(alpha);
    }
    else {
      rotations >> sinus;
      rotations >> cosinus;
      rotations >> w;
      aff = compute_transformation_matrix(sinus, cosinus, w);
    }

    Nef_polyhedron N2 = N1;
    N2.transform(aff);
    N2.transform(Aff_transformation_3(CGAL::TRANSLATION, Vector_3(0,0,1)));

    cgal_nef3_timer_on = true;
#if defined CGAL_NEF3_SYMDIFF
  N1.symmetric_difference(N2);
#elif defined CGAL_NEF3_INTERSECTION
  N1.intersection(N2);
#elif defined CGAL_NEF3_DIFFERENCE
  N1.difference(N2);
#else
  N1.join(N2);
#endif
  }
}
