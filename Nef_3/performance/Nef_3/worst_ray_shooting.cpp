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

int main(int argc, char* argv[]) {
  
  CGAL_assertion(argc>1 && argc<6);

  int n = argc > 1 ? std::atoi(argv[1]) : 50;
  std::cerr << "runs: " << n << std::endl;
  int step = argc > 2 ? std::atoi(argv[2]) : 10;
  std::cerr << "step: " << step << std::endl;
  int s = argc > 3 ? std::atoi(argv[3]) : 100;
  std::cerr << "size of cube: " << s << std::endl;
  
  for(int i=step; i<=n*step; i+=step) {
    Nef_polyhedron C = create_complex_facet(i*8);

    std::ostringstream out;
    tgen t(out,s);
    t.create_tetrahedra(1,1,2*i);
    std::istringstream in(out.str());
    Nef_polyhedron NT;
    in >> NT;
    CGAL_assertion(NT.is_valid());
    
    NT.transform(Aff_transformation_3(CGAL::TRANSLATION,Vector_3(scale*2,0,-s*i)));
  
    cgal_nef3_timer_on = true;
    C+=NT;
    cgal_nef3_timer_on = false;
  }
}
