// file: examples/Nef_3/visualization_SM.C

#ifndef CGAL_USE_QT
#include <iostream>
int main(int, char*){
  std::cout << "Sorry, this demo needs QT..." << std::endl; return 0;}
#else
#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/IO/Qt_widget_Nef_S2.h>
#include <qapplication.h>

typedef CGAL::Homogeneous<CGAL::Gmpz> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron_3;
typedef Nef_polyhedron_3::Vertex_const_iterator Vertex_const_iterator;
typedef Nef_polyhedron_3::Nef_polyhedron_S2 Nef_polyhedron_S2;

int main(int argc, char* argv[]) {
  Nef_polyhedron_3 N;
  std::cin >> N;
  Vertex_const_iterator v = N.vertices_begin();
  Nef_polyhedron_S2 S(N.get_sphere_map(v));

  QApplication a(argc, argv);
  CGAL::Qt_widget_Nef_S2<Nef_polyhedron_S2>* w = 
    new CGAL::Qt_widget_Nef_S2<Nef_polyhedron_S2>(S);
  a.setMainWidget(w);
  w->show();
  return a.exec();
}
#endif
