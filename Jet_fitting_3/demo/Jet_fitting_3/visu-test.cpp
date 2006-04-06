//visu, lib
#include <qapplication.h>
#include <qmainwindow.h>
#include "Observer.h"
#include "widget.h"

//visu, local
#include "SketchSample.h"


//geom, local
#include "visu_poly.h"

#include <iostream>
#include <fstream>




void parse_facet(Facet_handle f, const Vertex_index& vi)
{
  int idx, n=0;

  Halfedge_around_facet_const_circulator hc = f->facet_begin();
  Halfedge_around_facet_const_circulator hc_end = hc;
  
  //std::cerr << "here\n";
  do
  {
    idx = vi[Vertex_const_iterator(hc->vertex())];
    ++hc; n++;
  }
  while (hc != hc_end);
  //std::cerr << n << '\n';
}


void parse_triangles(Mesh& m){
  Facet_handle f;

  Vertex_index vi(m.vertices_begin(), m.vertices_end());
  
  Facet_iterator fib = m.facets_begin(), fie = m.facets_end();
  for (; fib != fie; ++fib){
    f = fib;
    parse_facet(f, vi);
  }
}

void load_geom(int argc, char* argv[])
{	
  assert(argc==2);
  const char* file_off = argv[1];
  //const char* file_res=argv[2];
 
  Mesh m;


  // initialisation du maillage
  std::ifstream f(file_off, std::ifstream::in);
  if(!f){
      exit(0);
    }
  f >> m;
  parse_triangles(m);
}





int main(int argc, char** argv) {

  load_geom(argc, argv);
  


  QApplication app(argc, argv);

  if ( !QGLFormat::hasOpenGL() ) {
    qWarning( "This system has no OpenGL support. Exiting." );
    return EXIT_FAILURE;
  }

  QMainWindow main;

  SketchSample* sketch = new SketchSample();
  Observer* obs = new Observer(sketch);

  GLWidget w1(&main, "OpenGL", obs);
  main.setCentralWidget( &w1 );
  main.resize( 500, 500 );
  app.setMainWidget(&main);

  main.show();

  return app.exec();
}
