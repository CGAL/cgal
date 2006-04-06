
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>

#include <iostream>
#include <fstream>


typedef  CGAL::Cartesian<double>               Kernel;
typedef  Kernel::Point_3                       Point;
typedef  CGAL::Polyhedron_3<Kernel>            Mesh;


typedef Kernel::FT	FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef Mesh::Vertex_handle Vertex_handle;
typedef Mesh::Facet_handle Facet_handle;
typedef Mesh::Vertex Vertex;
typedef Mesh::Facet Facet;
typedef Mesh::Face_handle Face_handle;
typedef Mesh::Halfedge_handle Halfedge_handle;
typedef Mesh::Facet_iterator Facet_iterator;
typedef Mesh::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
typedef Mesh::Halfedge_iterator Halfedge_iterator;
typedef Mesh::Point_iterator Point_iterator;
typedef Mesh::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
typedef Mesh::Vertex_iterator Vertex_iterator;
typedef Mesh::Vertex_const_iterator   Vertex_const_iterator;

typedef Mesh::Edge_iterator Edge_iterator ;
typedef Mesh::Halfedge_around_facet_const_circulator Halfedge_around_facet_const_circulator;

typedef CGAL::Inverse_index < Vertex_const_iterator > Vertex_index;


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

int main(int argc, char* argv[])
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

