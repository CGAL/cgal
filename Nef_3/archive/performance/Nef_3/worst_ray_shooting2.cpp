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

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Nef_3/Polygon_constructor.h>
#include "tetrahedron_generator.h"

#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>

typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Point_3 Point_3;
typedef Nef_polyhedron::Vector_3 Vector_3;
typedef Nef_polyhedron::Aff_transformation_3 Aff_transformation_3;
typedef tetrahedron_generator<Kernel> tgen;

const int scale = 1000000000;
const double PI = 3.1415926; // 53589793238462643383280;

bool cgal_nef3_timer_on = false;

Nef_polyhedron create_complex_facet(int n) {
  
  typedef std::list<Point_3> pointlist;
  typedef pointlist::const_iterator pointiterator;
  typedef std::pair<pointiterator,pointiterator> pointrange;
  typedef std::list<pointrange> polygonlist;
  typedef polygonlist::const_iterator pi;

  typedef CGAL::Polygon_constructor<Nef_polyhedron, pi> Polygon_constructor;

  CGAL_assertion(n > 2);
  polygonlist poly;
  pointlist points;

  for(int i = 0; i<n; ++i)
    points.push_back(Point_3(sin(2*PI*i/n) * scale, cos(2*PI*i/n) * scale, 0));
  poly.push_back(pointrange(points.begin(),points.end()));

  Polygon_constructor pc(poly.begin(), poly.end());
  Nef_polyhedron N;
  N.delegate(pc,true);
  return N;
}

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

  Aff_transformation_3 aff( cos_alpha, NT(0), sin_alpha,
			    NT(0), w, NT(0),
			    -sin_alpha, NT(0), cos_alpha,
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
  
  Aff_transformation_3 aff( cosinus, NT(0), sinus,
			    NT(0), w, NT(0),
			    -sinus, NT(0), cosinus,
			    w);
  return aff;

}

int main(int argc, char* argv[]) {
  
  CGAL_assertion(argc>2 && argc<6);

  std::ifstream rotations(argv[1]);
  CGAL_assertion_msg(rotations, "incorrect parameter 2");

  std::string mode;
  rotations >> mode;
  CGAL_assertion(mode == "angle" || mode == "sinus");
  int dummy;
  rotations >> dummy;

  int n = argc > 2 ? std::atoi(argv[2]) : 50;
  std::cerr << "runs: " << n << std::endl;
  int step = argc > 3 ? std::atoi(argv[3]) : 10;
  std::cerr << "step: " << step << std::endl;
  int s = argc > 4 ? std::atoi(argv[4]) : 100;
  std::cerr << "size of cube: " << s << std::endl;
  
  double alpha;
  RT sinus, cosinus, w;
  Aff_transformation_3 aff;

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
  
  for(int i=step; i<=n*step; i+=step) {
    Nef_polyhedron C = create_complex_facet(i*8);
    C.transform(aff);

    std::ostringstream out;
    tgen t(out,s);
    t.create_tetrahedra(1,2*i,1);
    std::istringstream in(out.str());
    Nef_polyhedron NT;
    in >> NT;
    CGAL_assertion(NT.is_valid());
    
    NT.transform(Aff_transformation_3(CGAL::TRANSLATION,Vector_3(scale*2,(-s*i)/2,-s)));
  
    cgal_nef3_timer_on = true;
    C+=NT;
    cgal_nef3_timer_on = false;
  }
}
