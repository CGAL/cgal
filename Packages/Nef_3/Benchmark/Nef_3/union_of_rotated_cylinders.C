#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/leda_integer.h>
#include <CGAL/gmpxx.h>
#include <CGAL/Simple_homogeneous.h>
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
#include <string>

typedef leda_integer                           NT;
// typedef CGAL::Gmpz                             NT;
// typedef mpz_class                        NT;
typedef CGAL::Simple_homogeneous<NT>           Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel>         Nef_polyhedron;
typedef CGAL::Polyhedron_3<Kernel>             Polyhedron;
typedef Kernel::Aff_transformation_3           Aff_transformation_3;


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

Aff_transformation_3 compute_transformation_matrix(NT sinus, NT cosinus, NT w) {

  double sin_double = sinus.to_double() / w.to_double();
  double arc = std::asin(sin_double);
  double alpha = arc * 180 / M_PI;

  std::cerr << std::endl << "sin(alpha)*w = " << sinus; 
  std::cerr << std::endl << "cos(alpha)*w = " << cosinus; 
  std::cerr << std::endl << "w            = " << w;

  std::cerr << std::endl << "sin_double   = " << sin_double;
  std::cerr << std::endl << "arc          = " << arc;
  std::cerr << std::endl << "alpha        = " << alpha;
  
  Aff_transformation_3 aff( cosinus,-sinus, NT(0),
			    sinus, cosinus, NT(0),
			    NT(0), NT(0), w,
			    w);
  return aff;

}

int main(int argc, char* argv[]) {

  CGAL_assertion(argc == 3);
  
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
  NT sinus, cosinus, w;
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
    CGAL::Timer t;
    t.start();
    N2 = N2.join(N1);
    t.stop();
    
    std::cerr << std::endl << "time         = " << t.time();
    std::cerr << std::endl;
  }
}
