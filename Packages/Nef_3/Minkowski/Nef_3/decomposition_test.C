#include <algorithm>
#include <map>
#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/leda_integer.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>
#include <CGAL/Nef_3/Single_wall_creator.h>
#include <CGAL/Nef_3/External_structure_builder.h>
#include <CGAL/Nef_3/SNC_io_parser.h>
#include <fstream>

typedef leda_integer NT;
typedef CGAL::Homogeneous<NT> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel>     Nef_polyhedron;
typedef Nef_polyhedron::SNC_structure  SNC_structure;
typedef CGAL::SNC_decorator<SNC_structure>  SNC_decorator;
typedef SNC_structure::Halfedge_handle  Halfedge_handle;
typedef Nef_polyhedron::Halfedge_iterator  Halfedge_iterator;
typedef Nef_polyhedron::Vector_3           Vector_3;

typedef CGAL::Single_wall_creator<Nef_polyhedron> Single_wall;
typedef CGAL::External_structure_builder<Nef_polyhedron> External_structure_builder;

int main(int argc, char* argv[]) {

  CGAL_assertion(argc==2);
  std::ifstream in(argv[1]);
  Nef_polyhedron N;
  in >> N;
  CGAL_NEF_SETDTHREAD(43);
  std::cerr << N;
  SNC_decorator D(*const_cast<SNC_structure*>(N.sncp()));

  int i=0;
  Halfedge_iterator e = D.halfedges_begin();
  for(;e!=D.halfedges_end();++e) {
    if(e->is_twin() && e->twin() != Halfedge_handle()) {
      std::cerr << "edge: " << e->source()->point();
      std::cerr << "->" 
		<< e->twin()->source()->point() << std::endl;
    
    Single_wall W(e,Vector_3(-1,0,0));
    N.delegate(W);
    }
  }

  /* 
  for(D.halfedges_begin();e!=D.halfedges_end();++e) {
    if(e->is_twin()) {

    std::cerr << "edge: " << e->source()->point() << "->" 
	      << e->twin()->source()->point() << std::endl;
    
    Single_wall W(e,Vector_3(1,0,0));
    N.delegate(W);
    }
  }
  */

  External_structure_builder esb;
  N.delegate(esb);

  QApplication a(argc, argv);
  CGAL::Qt_widget_Nef_3<Nef_polyhedron>* w = 
    new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(N);
  a.setMainWidget(w);
  w->show();
  a.exec();
}
