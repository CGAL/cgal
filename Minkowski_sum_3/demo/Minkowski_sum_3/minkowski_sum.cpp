#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Lazy_kernel.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/minkowski_sum_3.h>
#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>
#include <fstream>
#include <iostream>

typedef CGAL::Lazy_kernel<CGAL::Simple_cartesian<CGAL::Gmpq> > Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel,CGAL::SNC_indexed_items>     Nef_polyhedron;
typedef CGAL::Polyhedron_3<Kernel>     Polyhedron;

bool loadFile(char* filename, Nef_polyhedron& N) {
  std::ifstream in(filename);
  std::ifstream test(filename);
  char c;
  test >> c;
  if(c!='S' && c!='O') return false;
  if(c == 'S')
    in >> N;
  else {
    Polyhedron P;
    in >> P;
    N = Nef_polyhedron(P);
  }
  return true;
}

int main(int argc, char* argv[]) {

  CGAL_assertion(argc==3);

  Nef_polyhedron N0, N1;
  if(!loadFile(argv[1], N0)) {
    std::cerr << "parameter 1 is not a valid input file" << std::endl;
    return 0;
  }
  if(!loadFile(argv[2], N1)) {
    std::cerr << "parameter 2 is not a valid input file" << std::endl;
    return 0;
  }

  Nef_polyhedron result = CGAL::minkowski_sum_3(N0, N1);

  QApplication a(argc, argv);
  CGAL::Qt_widget_Nef_3<Nef_polyhedron>* w = 
    new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(result);
  a.setMainWidget(w);
  w->show();
  a.exec();
}
