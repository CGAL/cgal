#ifndef CGAL_USE_QT
#include <iostream>
int main(int, char*){
  std::cout << "Sorry, this demo needs QT..." << std::endl; return 0;}
#else
#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Gmpz.h>
#include <CGAL/random_selection.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/Nef_S2/SM_items.h>
#include <CGAL/IO/Qt_widget_Nef_S2.h>

typedef CGAL::Gmpz NT;
typedef CGAL::Homogeneous<NT> Kernel;
typedef Kernel::Point_3       Point_3;
typedef Kernel::Plane_3       Plane_3;

typedef CGAL::Nef_polyhedron_S2<Kernel> Nef_polyhedron_S2;
typedef Nef_polyhedron_S2::Sphere_circle  Sphere_circle;

typedef CGAL::Creator_uniform_3<NT,Point_3>  Creator;
typedef CGAL::Random_points_in_cube_3<Point_3,Creator> Point_source;


int main(int argc, char **argv) {
  int n(5), r(0);
  if ( argc > 1 ) n = atoi( argv[1] );  // reads number of Sphere_segments
  if ( argc > 2 ) r = atoi( argv[2] );  // reads seed for random points
  srand(r);

  std::list<Sphere_circle> L;
  Point_source S(5);
  Point_3 ph;
  Point_3 o(0,0,0);
  while ( n-- > 0 ) {
    do { ph = *S++; } while ( ph == o );
    Plane_3 h(o,(ph-CGAL::ORIGIN).direction());
    L.push_back( Sphere_circle(h) );
  }
  
  // partition input into two lists
  Nef_polyhedron_S2 Ni, N;
  bool first(true);
  std::list<Sphere_circle>::const_iterator it;
  CGAL_forall_iterators(it,L) {
    if ( first ) {
      N = Nef_polyhedron_S2(*it);
      first = false;
    } else {
      Ni = Nef_polyhedron_S2(*it);
      N = N.symmetric_difference(Ni);
    }
  }

  QApplication a(argc, argv);
  CGAL::Qt_widget_Nef_S2<Nef_polyhedron_S2>* w = 
    new CGAL::Qt_widget_Nef_S2<Nef_polyhedron_S2>(N);
  a.setMainWidget(w);
  w->show();
  return a.exec();
}
#endif

