#include <CGAL/basic.h>

#include <cstring>
#include <iostream>
#include <fstream>
#include <strstream.h>

#include <list>

#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>

#include <CGAL/Triangulation_iterators_3.h>
#include <CGAL/Triangulation_circulators_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_geom_traits_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/IO/Geomview_stream.h>

typedef CGAL::Cartesian<double>  Rep;

typedef CGAL::Triangulation_geom_traits_3<Rep> Gt;
typedef CGAL::Triangulation_vertex_base_3<Gt> Vb;
typedef CGAL::Triangulation_cell_base_3<Gt>  Cb;

typedef CGAL::Triangulation_data_structure_3<Vb,Cb> TDS;
typedef CGAL::Triangulation_3<Gt,TDS> Triangulation;
typedef CGAL::Delaunay_triangulation_3<Gt,TDS> Delaunay;

typedef CGAL::Triangulation_vertex_iterator_3<Gt,TDS> Vertex_iterator;
typedef CGAL::Triangulation_edge_iterator_3<Gt,TDS> Edge_iterator;
typedef CGAL::Triangulation_cell_iterator_3<Gt,TDS> Cell_iterator;
typedef CGAL::Triangulation_facet_iterator_3<Gt,TDS> Facet_iterator;
typedef CGAL::Triangulation_cell_circulator_3<Gt,TDS> Cell_circulator;

typedef typename Triangulation::Cell Cell;
typedef typename Triangulation::Vertex Vertex;
typedef typename Triangulation::Cell_handle Cell_handle;
typedef typename Triangulation::Vertex_handle Vertex_handle;
typedef typename Triangulation::Locate_type Locate_type;

typedef Gt::Point Point;
//typedef CGAL::Point_3<Rep>  Point;

////////////////////// 
// VISU GEOMVIEW
////////////////////// 
template<class TRIANGULATION>
void visu_cells(CGAL::Geomview_stream & os, const TRIANGULATION & T)
{
  Cell_iterator cit = T.finite_cells_begin();
  Cell_iterator cdone = T.cells_end();
  
  if ( cit == cdone ) { cout << "debut=fin" << endl ;}
  else {
    while(cit != cdone) {
      os << T.tetrahedron(&(*cit));
    ++cit;
    }
  }
}
void visu_cell(CGAL::Geomview_stream & os, Cell_handle c)
{
  os << Gt::Tetrahedron(c->vertex(0)->point(),
			   c->vertex(1)->point(),
			   c->vertex(2)->point(),
			   c->vertex(3)->point());
}
template<class TRIANGULATION>
void visu_facets(CGAL::Geomview_stream & os, const TRIANGULATION & T)
{
  Facet_iterator fit = T.finite_facets_begin();
  Facet_iterator fdone = T.facets_end();
  
  if ( fit == fdone ) { cout << "debut=fin" << endl ;}
  else {
    while(fit != fdone) {
      os << T.triangle(*fit);
      ++fit;
    }
  }
}
template<class TRIANGULATION>
void visu_edges(CGAL::Geomview_stream & os, const TRIANGULATION & T)
{
  Edge_iterator eit = T.finite_edges_begin();
  Edge_iterator edone = T.edges_end();
  
  if ( eit == edone ) { cout << "debut=fin" << endl ;}
  else {
    while(eit != edone) {
      os << T.segment(*eit);
      ++eit;
    }
  }
}
template<class TRIANGULATION>
void visu_vertices(CGAL::Geomview_stream & os, const TRIANGULATION & T)
{
  Vertex_iterator vit = T.finite_vertices_begin();
  Vertex_iterator vdone = T.vertices_end();
  
  if ( vit == vdone ) { cout << "debut=fin" << endl ;}
  else {
    while(vit != vdone) {
      os << vit->point();
      ++vit;
    }
  }
}

////////////////////// 
// INSERTION
////////////////////// 

template<class TRIANGULATION>
void insere(CGAL::Geomview_stream & os, TRIANGULATION & T, Point p)
{
  cout << p << endl;
  os << p;
  cout << "localisation" << endl;
  Triangulation::Locate_type lt;
  int li, lj;
  Cell_handle c = T.locate( p, lt, li, lj ) ;
  switch ( T.dimension() ) {
  case 0:
    {
      pp_vertex(c->vertex(0));
      break;
    }
  case 1:
    {
      pp_edge(CGAL::make_triple(c,0,1));
      break;
    }
  case 2:
    {
      pp_facet(make_pair(c,3));
      break;
    }
  case 3:
    {
      pp_cell(c);
      break;
    }
  }
  cout << (int) lt << " " << li << " " << lj << endl;
  cout << "insertion " << endl;
  T.insert( p );  
  affiche_sommets(T);
  cout << "validite " << T.is_valid(true) << endl;
}


CGAL::Geomview_stream gv(CGAL::Bbox_3(0,0,0, 2, 2, 2));
Delaunay T;

int main(int argc, char* argv[])
{

  CGAL::Geomview_stream gv(CGAL::Bbox_3(0,0,0, 2, 2, 2));

  gv.set_line_width(4);
  gv.set_trace(false);
  gv.set_bg_color(CGAL::Color(200, 200, 200));
  gv.set_face_color(CGAL::RED);
  gv.set_edge_color(CGAL::GREEN);
  gv.set_vertex_color(CGAL::BLUE);

  Point p0(0,0,0);
  Point px(1,0,0);
  Point py(0,1,0);
  Point pz(0,0,1);
  
  ifstream iFile("data",ios::in);
  if (iFile) cout <<"                              reading file "
		  << "data" << endl ;
  Point nouv;
  if (iFile) {
    while ( iFile >> nouv ) {
      T.insert(nouv);
    }
  }

  visu_cells(gv,T);
  visu_vertices(gv,T);
  visu_edges(gv,T);

  cout << T.is_valid(true);

  ofstream oFileT("output",ios::out);
  cout <<"                              writing file "
       << "output" << endl << flush;
  oFileT << T;

  char ch;
  cout << "donner caractere de fin" << endl;
  cin >> ch;

  return 1;
}
