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
#include <LEDA/leda_rational.h>
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
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/rational_rotation.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cmath>
#include <cstddef>

typedef CGAL::Nef_polyhedron_3<Kernel>      Nef_polyhedron;
typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;
typedef Kernel::Aff_transformation_3        Aff_transformation_3;


Aff_transformation_3 compute_transformation_matrix(double alpha) {
  
  double arc = M_PI * alpha / 180.0;

  NT epsilon = 1;
  double sin_double = std::sin( arc);
  double cos_double = std::cos( arc);

  std::cerr << std::endl << "alpha        = " << alpha;
  std::cerr << std::endl << "arc          = " << arc;
  std::cerr << std::endl << "sin_double " << sin_double;
  std::cerr << std::endl << "cos_double " << cos_double;

  while(sin_double < 1000 || cos_double < 1000) {
    sin_double *= 10;
    cos_double *= 10;
    epsilon *= 10;
  }

  std::cerr << std::endl << "epsilon      = 1/" << epsilon; 

  NT sin_alpha(0);
  NT cos_alpha(0);
  NT w(0);
  
  CGAL::Timer t;
  t.start(); 
  CGAL::rational_rotation_approximation( arc,
					 sin_alpha, cos_alpha, w,
					 NT(1), NT(epsilon));
  t.stop();
  std::cerr << std::endl << "approx. time = " << t.time();

  Aff_transformation_3 aff( cos_alpha,-sin_alpha, NT(0),
			    sin_alpha, cos_alpha, NT(0),
			    NT(0), NT(0), w,
			    w);

  std::cerr << std::endl << "sin(alpha)*w = " << sin_alpha; 
  std::cerr << std::endl << "cos(alpha)*w = " << cos_alpha; 
  std::cerr << std::endl << "w            = " << w;

  return aff;
}

int main() {

  std::ifstream in("off/ngon1000.off");
  Polyhedron poly;
  in >> poly;
  Nef_polyhedron N1(poly);

  double alpha = 1;
  for(int i=0; i<10; i++) {

    Aff_transformation_3 aff = compute_transformation_matrix(alpha);
    Nef_polyhedron N2 = N1;
    N2.transform(aff);
    CGAL::Timer t;
    t.start();
    N2 = N2.join(N1);
    t.stop();

    CGAL_assertion(N2.number_of_vertices() == 800);

    std::cerr << std::endl << "time         = " << t.time();
    alpha/=10;
    std::cerr << std::endl;
  }
}
