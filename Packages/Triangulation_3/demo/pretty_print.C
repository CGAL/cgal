
#include <CGAL/basic.h>
#include <stdio.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>
#include <strstream.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>

#include <CGAL/Triangulation_geom_traits_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_ds_iterators_3.h>

#include <CGAL/Triangulation_iterators_3.h>
#include <CGAL/Triangulation_circulators_3.h>
#include <CGAL/Triangulation_3.h>

#include <CGAL/IO/Geomview_stream.h>

typedef double coord_type;
typedef CGAL_Cartesian<coord_type>  Rep;

typedef CGAL_Point_3<Rep>  Point;

typedef CGAL_Triangulation_geom_traits_3<Rep> Gt;
typedef CGAL_Triangulation_vertex_base_3<Gt> Vb;
typedef CGAL_Triangulation_cell_base_3<Gt>  Cb;

typedef CGAL_Triangulation_data_structure_3<Vb,Cb> TDS;

typedef CGAL_Triangulation_ds_vertex_iterator_3<TDS> Vertex_iterator;
typedef CGAL_Triangulation_ds_edge_iterator_3<TDS> Edge_iterator;
typedef CGAL_Triangulation_ds_facet_iterator_3<TDS> Facet_iterator;

typedef typename TDS::Cell TDSCell;
typedef typename TDS::Facet TDSFacet;
typedef typename TDS::Edge TDSEdge;
typedef typename TDS::Vertex TDSVertex;

typedef CGAL_Triangulation_3<Gt,TDS> Triangulation_3;

typedef typename Triangulation_3::Cell Cell;
typedef typename Triangulation_3::Facet Facet;
typedef typename Triangulation_3::Edge Edge;
typedef typename Triangulation_3::Vertex Vertex;

typedef typename Triangulation_3::Vertex_handle Vertex_handle;
typedef typename Triangulation_3::Cell_handle Cell_handle;

void
pp_point(const Point & p)
{
  cerr << p << endl;
}

void
pp_tds_vertex(const TDSVertex* v)
{
  cerr  << v->point() << endl;
}

void
pp_tds_edge(const TDSEdge e)
{
  cerr  << e.first->vertex(e.second)->point() << ", "
	<< e.first->vertex(e.third)->point() << endl;
}

void
pp_tds_3_facet(const TDSFacet f)
{
  int i=f.second;
  if (i != 0) {
  cerr  << f.first->vertex(0)->point() << ", " ;
  }
  if (i != 1) {
  cerr  << f.first->vertex(1)->point() << ", " ;
  }
  if (i != 2) {
    cerr  << f.first->vertex(2)->point() << ", " ;
  }
  if (i != 3) {
    cerr  << f.first->vertex(3)->point() << ", " ;
  }
  cerr << endl;
}

void
pp_tds_cell(const TDSCell* f)
{
  cerr  << f->vertex(0)->point() << ", "
	<< f->vertex(1)->point() << ", "
	<< f->vertex(2)->point() << ", "
	<< f->vertex(3)->point() << endl;
}

void
pp_vertex(const Vertex_handle v)
{
  cerr  << v->point() << endl;
}

void
pp_edge(const Edge e)
{
  cerr  << e.first->vertex(e.second)->point() << ", "
	<< e.first->vertex(e.third)->point() << endl;
}
void
pp_edge(const Cell_handle c, int i, int j)
{
  pp_edge(CGAL_make_triple(c,i,j));
}

void
pp_facet(const Facet f)
{
  int i=f.second;
  if (i != 0) {
  cerr  << f.first->vertex(0)->point() << ", " ;
  }
  if (i != 1) {
  cerr  << f.first->vertex(1)->point() << ", " ;
  }
  if (i != 2) {
    cerr  << f.first->vertex(2)->point() << ", " ;
  }
  if (i != 3) {
    cerr  << f.first->vertex(3)->point() << ", " ;
  }
  cerr << endl;
}
void
pp_facet(const Cell_handle c, int i)
{
  pp_facet(make_pair(c,i));
}
void
pp_cell(const Cell_handle c)
{
  cerr  << c->vertex(0)->point() << ", "
	<< c->vertex(1)->point() << ", "
	<< c->vertex(2)->point() << ", "
	<< c->vertex(3)->point() << endl;
}

void visu_face(CGAL_Geomview_stream & os, Triangulation_3 & T, 
	       Cell_handle c, int i)
{
  if ( ! T.is_infinite(c,i) ) {
    os.set_face_color(CGAL_GREEN);
    os << T.triangle(make_pair(c,i));
  }
  else {
    os.set_vertex_color(CGAL_GREEN);
    cout << "face infinie" << endl;
    for (int j=0; j<4; j++) {
      if (j!=i) 
	{ os << c->vertex(j)->point();}
    }
  }
}
