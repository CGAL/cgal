//visu, lib
#include <qapplication.h>
#include <qmainwindow.h>
#include "Observer.h"
#include "widget.h"

//visu, local
#include "SketchSample.h"


//geom, local
#include "visu_poly.h"
//marc
#include "enriched_polyhedron.h" 

#include <iostream>
#include <fstream>

Mesh m_mesh;

DS ppal_data;

void clean_DS(DS& L)
{
  DS_iterator current=L.begin();
  for (DS_iterator it=++L.begin();it!=L.end();it++)
    {
      delete *current;
      L.pop_front();
      current=it;
    }
  delete *current;
  L.pop_front();
}

void 
read_line(FILE* file, Point& P1,Point& P2,
	  Vector&  D1,Vector& D2, double& k1,double& k2)
{
  double x,y,z;
  fscanf(file, "%lf%lf%lf",&x,&y,&z);
  P1=Point(x,y,z);
  fscanf(file, "%lf%lf%lf",&x,&y,&z);
  P2=Point(x,y,z);
  fscanf(file, "%lf%lf%lf",&x,&y,&z);
  D1=Vector(x,y,z);
  fscanf(file, "%lf%lf%lf",&x,&y,&z);
  D2=Vector(x,y,z);	
  fscanf(file, "%lf%lf",&k1,&k2);
}

void load_data_from_file(DS& l, char* file_res)
{
  FILE *file;
  if((file = fopen(file_res, "r")) != NULL) 
    {
      while (!feof(file))
	{
	  Point P1,P2;	
	  Vector D1,D2;
	  double k1,k2;
	  read_line(file,P1,P2,D1,D2,k1,k2);
	  if (feof(file)) break;
	  l.push_front(new data_line(P1,P2,D1,D2,k1,k2));
	}
      fclose(file);
    }
  else
    std::cout << "Cannot open file" << std::endl;	
}	


void load_geom(int argc, char* argv[])
{	
  assert(argc==3);
  const char* file_off = argv[1];
  const char* file_res = argv[2];
 
  // initialisation du maillage
  std::ifstream f(file_off, std::ifstream::in);
  if(!f){
      exit(0);
    }
  f >> m_mesh;
  m_mesh.compute_normals();
  
  // init of the ppal dir
  ppal_data.clear();
  load_data_from_file(ppal_data, file_res);

}





int main(int argc, char** argv) {

  load_geom(argc, argv);
  
  QApplication app(argc, argv);

  if ( !QGLFormat::hasOpenGL() ) {
    qWarning( "This system has no OpenGL support. Exiting." );
    return EXIT_FAILURE;
  }

  QMainWindow main;

  SketchSample* sketch = new SketchSample(&m_mesh, &ppal_data);
  Observer* obs = new Observer(sketch);

  GLWidget w1(&main, "OpenGL", obs);
  main.setCentralWidget( &w1 );
  main.resize( 500, 500 );
  app.setMainWidget(&main);

  main.show();

  return app.exec();
}



// void parse_facet(Facet_handle f, const Vertex_index& vi)
// {
//   int idx, n=0;

//   Halfedge_around_facet_const_circulator hc = f->facet_begin();
//   Halfedge_around_facet_const_circulator hc_end = hc;
  
//   //std::cerr << "here\n";
//   do
//   {
//     idx = vi[Vertex_const_iterator(hc->vertex())];
//     ++hc; n++;
//   }
//   while (hc != hc_end);
//   //std::cerr << n << '\n';
// }


// void parse_triangles(Mesh& m){
//   Facet_handle f;

//   Vertex_index vi(m.vertices_begin(), m.vertices_end());
  
//   Facet_iterator fib = m.facets_begin(), fie = m.facets_end();
//   for (; fib != fie; ++fib){
//     f = fib;
//     parse_facet(f, vi);
//   }
// }
