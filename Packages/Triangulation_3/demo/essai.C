#include <CGAL/basic.h>

#include <stdio.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>
#include <strstream.h>

#include <list.h>

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

#include <CGAL/Triangulation_short_names_3.h>

#include <CGAL/IO/Geomview_stream.h>
#include "pretty_print.h"

typedef double coord_type;
typedef CGAL_Cartesian<coord_type>  Rep;
typedef CGAL_Triangulation_geom_traits_3<Rep> Gt;

typedef CGAL_Triangulation_vertex_base_3<Gt> Vb;
typedef CGAL_Triangulation_cell_base_3<Gt>  Cb;

typedef CGAL_Triangulation_data_structure_3<Vb,Cb> TDS;

typedef typename TDS::Cell TDSCell;
typedef typename TDS::Vertex TDSVertex;
typedef typename TDS::Edge TDSEdge;

typedef CGAL_Triangulation_ds_vertex_iterator_3<TDS> TDSVertex_iterator;
typedef CGAL_Triangulation_ds_edge_iterator_3<TDS> TDSEdge_iterator;
typedef CGAL_Triangulation_ds_cell_iterator_3<TDS> TDSCell_iterator;
typedef CGAL_Triangulation_ds_facet_iterator_3<TDS> TDSFacet_iterator;
typedef CGAL_Triangulation_ds_cell_circulator_3<TDS> TDSCell_circulator;

typedef CGAL_Triangulation_3<Gt,TDS> Triangulation_3;
typedef CGAL_Delaunay_triangulation_3<Gt,TDS> Delaunay;

typedef CGAL_Triangulation_vertex_iterator_3<Gt,TDS> Vertex_iterator;
typedef CGAL_Triangulation_edge_iterator_3<Gt,TDS> Edge_iterator;
typedef CGAL_Triangulation_cell_iterator_3<Gt,TDS> Cell_iterator;
typedef CGAL_Triangulation_facet_iterator_3<Gt,TDS> Facet_iterator;
typedef CGAL_Triangulation_cell_circulator_3<Gt,TDS> Cell_circulator;

typedef typename Triangulation_3::Cell Cell;
typedef typename Triangulation_3::Vertex Vertex;
typedef typename Triangulation_3::Cell_handle Cell_handle;
typedef typename Triangulation_3::Vertex_handle Vertex_handle;
typedef typename Triangulation_3::Locate_type Locate_type;

typedef Gt::Point Point;

////////////////////// 
// TDS
////////////////////// 

void affiche_sommets(TDS& T)
{
  cout << "sommets TDS : " << endl ;
  TDSVertex_iterator vit = T.vertices_begin(), vdone = T.vertices_end();

  if ( vit == vdone ) { cout << "debut=fin" << endl ;}
  else {
    while(vit != vdone) {
      cout << vit->point() << endl;
      ++vit;
    }
  }
  cout << endl;
}
void affiche_cellules(TDS & T)
{
  cout << "cellules TDS : " << endl ;
  TDSCell_iterator it = T.cells_begin();
  TDSCell_iterator done = T.cells_end();

  if ( it == done ) { cout << "debut=fin" << endl ;}
  else {
    while(it != done) {
      cout << it->vertex(0)->point() << ", "
	   << it->vertex(1)->point() << ", "
	   << it->vertex(2)->point() << ", "
	   << it->vertex(3)->point() << endl;
      ++it;
    }
  }
  cout << endl;
}
void affiche_validite(TDS& T)
{
  bool b = T.is_valid(true,0);
  cout << "validite " << b << endl << endl;
}

////////////////////// 
// Triangulation
////////////////////// 
template<class TRIANGULATION>
void affiche_validite(TRIANGULATION& T)
{
  bool b = T.is_valid(true,0);
  cout << "validite " << b << endl << endl;
}
template<class TRIANGULATION>
void affiche_sommets(TRIANGULATION& T)
{
  cout << "sommets : " << endl ;
  Vertex_iterator vit = T.all_vertices_begin(), vdone = T.vertices_end();

  if ( vit == vdone ) { cout << "debut=fin" << endl ;}
  else {
    while(vit != vdone) {
      cout << vit->point() << endl;
      ++vit;
    }
  }
  cout << endl;
}
template<class TRIANGULATION>
void affiche_aretes_plus(TRIANGULATION& T)
{
  cout << "aretes plus : " << endl ;
  Edge_iterator eit = T.all_edges_begin(), edone = T.edges_end();
  if ( eit == edone ) { cout << "debut=fin" << endl ;}
  else {
    while(eit != edone) {
      cout << (*eit).first->vertex((*eit).second)->point() << ", "
	   << (*eit).first->vertex((*eit).third)->point() << endl;
      ++eit;
    }
  }
  cout << endl;
}
template<class TRIANGULATION>
void affiche_aretes_finies(TRIANGULATION& T)
{
  cout << "aretes plus : " << endl ;
  Edge_iterator eit = T.finite_edges_begin(), edone = T.edges_end();
  if ( eit == edone ) { cout << "debut=fin" << endl ;}
  else {
    while(eit != edone) {
      cout << (*eit).first->vertex((*eit).second)->point() << ", "
	   << (*eit).first->vertex((*eit).third)->point() << endl;
      ++eit;
    }
  }
  cout << endl;
}
template<class TRIANGULATION>
void affiche_aretes_minus(TRIANGULATION& T)
{
  cout << "aretes minus : " << endl ;
  Edge_iterator edone = T.all_edges_begin(), eit = T.edges_end();
  if ( eit == edone ) { cout << "debut=fin" << endl ;}
  else {
    while(eit != edone) {
      --eit;
      cout << (*eit).first->vertex((*eit).second)->point() << ", "
	   << (*eit).first->vertex((*eit).third)->point() << endl;
    }
  }
  cout << endl;
}
template<class TRIANGULATION>
void affiche_faces_plus(TRIANGULATION& T)
{
  cout << "faces plus : " << endl ;
  Facet_iterator fit = T.all_facets_begin(), fdone = T.facets_end();
  int i;
  if ( fit == fdone ) { cout << "debut=fin" << endl ;}
  else {
    while(fit != fdone) {
      for ( i=0; i<4 ; i++ ) {
	if ( (*fit).second != i )
	  { cout << (*fit).first->vertex(i)->point() << ", "; }
      }
      cout << endl;
      ++fit;
    }
  }
  cout << endl;
}
template<class TRIANGULATION>
void affiche_faces_minus(TRIANGULATION& T)
{
  cout << "faces minus : " << endl ;
  Facet_iterator fdone = T.all_facets_begin(), fit = T.facets_end();
  int i;
  if ( fit == fdone ) { cout << "debut=fin" << endl ;}
  else {
    while(fit != fdone) {
      --fit;
      for ( i=0; i<4 ; i++ ) {
	if ( (*fit).second != i )
	  { cout << (*fit).first->vertex(i)->point() << ", "; }
      }
      cout << endl;
    }
  }
  cout << endl;
}
template<class TRIANGULATION>
void affiche_faces_voisins(TRIANGULATION& T)
{
  cout << "faces (dim 2 seulement) : " << endl ;
  Facet_iterator fit = T.all_facets_begin(), fdone = T.facets_end();
  int i;
  if ( fit == fdone ) { cout << "debut=fin" << endl ;}
  else {
    while(fit != fdone) {
      for ( i=0; i<3 ; i++ ) {
	cout << (*fit).first->vertex(i)->point() << ", "; 
      }
      cout << endl << "    voisins " << endl;
      for ( i=0; i<3 ; i++ ) {
	for (int j=0; j<3; j++ ) {
	  cout << (*fit).first->neighbor(i)->vertex(j)->point() << ", "; 
	}
	cout << endl;
      }
      cout << endl;
      ++fit;
    }
  }
  cout << endl;
}
template<class TRIANGULATION>
void affiche_cellules( TRIANGULATION& T)
{
  cout << "cellules : " << endl ;
  Cell_iterator it = T.all_cells_begin();
  Cell_iterator done = T.cells_end();

  if ( it == done ) { cout << "debut=fin" << endl ;}
  else {
    while(it != done) {
      cout << it->vertex(0)->point() << ", "
	   << it->vertex(1)->point() << ", "
	   << it->vertex(2)->point() << ", "
	   << it->vertex(3)->point() << endl;
      ++it;
    }
  }
  cout << endl;
}

void visu_cellules(CGAL_Geomview_stream & os, Triangulation_3 & T)
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
void visu_cellule(CGAL_Geomview_stream & os, Cell_handle c)
{
  os << Gt::Tetrahedron(c->vertex(0)->point(),
			   c->vertex(1)->point(),
			   c->vertex(2)->point(),
			   c->vertex(3)->point());
}
void visu_faces(CGAL_Geomview_stream & os, Triangulation_3 & T)
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
void visu_aretes(CGAL_Geomview_stream & os, Triangulation_3 & T)
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
void visu_sommets(CGAL_Geomview_stream & os, Triangulation_3 & T)
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

void affiche_cellules_voisins(Triangulation_3& T)
{
  cout << "cellules : " << endl ;
  Cell_iterator cit = T.all_cells_begin(), cdone = T.cells_end();

  if ( cit == cdone ) { cout << "debut=fin" << endl ;}
  else {
    while(cit != cdone) {
      cout << endl
	   << cit->vertex(0)->point() << ", "
	   << cit->vertex(1)->point() << ", "
	   << cit->vertex(2)->point() << ", "
	   << cit->vertex(3)->point() << endl;
      for (int i=0; i<4; i++) {
	cout << "voisin " << i << endl
	     << cit->neighbor(i)->vertex(0)->point() << ", "
	     << cit->neighbor(i)->vertex(1)->point() << ", "
	     << cit->neighbor(i)->vertex(2)->point() << ", "
	     << cit->neighbor(i)->vertex(3)->point() << endl;
      }
      ++cit;
    }
  }
}

void insere_dehors(Triangulation_3& T, Point p)
{
  cout << p << endl;
  T.insert_outside_affine_hull( p );
  cout << T.number_of_vertices() << " sommets, dimension " 
       << T.dimension() << endl;
}

template<class TRIANGULATION>
void insere(CGAL_Geomview_stream & os, TRIANGULATION & T, Point p)
{
  cout << p << endl;
  os << p;
  cout << "localisation" << endl;
  Triangulation_3::Locate_type lt;
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
      pp_edge(CGAL_make_triple(c,0,1));
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

CGAL_Geomview_stream gv(CGAL_Bbox_3(0,0,0, 2, 2, 2));
Delaunay* D;

int main(int argc, char* argv[])
{

  //  Cell essaicell;

  //  CGAL_Geomview_stream gv(CGAL_Bbox_3(0,0,0, 2, 2, 2));
  gv.set_line_width(4);
  gv.set_trace(false);
  gv.set_bg_color(CGAL_Color(200, 200, 200));
  gv.set_face_color(CGAL_RED);
  gv.set_edge_color(CGAL_GREEN);
  gv.set_vertex_color(CGAL_BLUE);

  Point p0(0,0,0);
  //  pp_point(p0);
  Point px(1,0,0);
  //  pp_point(p0);
  Point py(0,1,0);
  //  pp_point(p0);
  Point pz(0,0,1);
  //  pp_point(p0);
  
  //  Delaunay T;
  //  D = &T;

  Triangulation_3 T;
  T.insert(p0);
  T.insert(px);
  T.insert(py);
  T.insert(pz);

  affiche_validite(T);
//   list<Point> l;
//   T.insert(l.begin(),l.end());
//   affiche_validite(T);
//   affiche_aretes_finies(T);

//   TDS tds=T.tds();
//   TDSVertex_iterator vit;
//   for ( vit = tds.vertices_begin(); vit != tds.vertices_end(); ++vit ) {
//     cout << "vertex ";
//     pp_tds_vertex(&(*vit));

//     set<TDSCell*,less<TDSCell*> > C;
//     tds.incident_cells(&(*vit), C);
//     set<TDSCell*, less<TDSCell*> >::const_iterator it;
//     cout << "incident cells " << endl ;
//     for ( it = C.begin(); it != C.end(); ++it ) {
//       pp_tds_cell(*it);
//     }
    
//     set<TDSVertex*,less<TDSVertex*> > V;
//     tds.incident_vertices(&(*vit), V);
//     set<TDSVertex*, less<TDSVertex*> >::const_iterator itv;
//     cout << "incident vertices " << endl ;
//     for ( itv = V.begin(); itv != V.end(); ++itv ) {
//       pp_tds_vertex(*itv);
//     }
//   }

  Vertex_iterator tvit;
  for ( tvit = T.all_vertices_begin(); tvit != T.vertices_end(); ++tvit ) {
    cout << "vertex ";
    pp_vertex(&(*tvit));

    set<Cell*,less<Cell*> > Ch;
    T.incident_cells(&(*tvit), Ch);
    set<Cell*, less<Cell*> >::const_iterator ith;
    cout << "incident cells " << endl ;
    for ( ith = Ch.begin(); ith != Ch.end(); ++ith ) {
      pp_cell((*ith)->handle());
    }

    set<Vertex*,less<Vertex*> > Vh;
    T.incident_vertices(&(*tvit), Vh);
    set<Vertex*, less<Vertex*> >::const_iterator itvh;
    cout << "incident vertices " << endl ;
    for ( itvh = Vh.begin(); itvh != Vh.end(); ++itvh ) {
      pp_vertex((*itvh)->handle());
    }
  }
//   T.insert(Point(0,0,1));
//   T.insert(Point(0,0,2));
//   T.insert(Point(0,0,3));
//   affiche_validite(T);
//   affiche_sommets(T);

  // test dim 3
  //  Triangulation_3 T(p0,px,py,pz);
//   T.insert_outside_affine_hull(p0);
//   T.insert_outside_affine_hull(px);
//   T.insert_outside_affine_hull(py);
//   T.insert_outside_affine_hull(pz);

  //Point q(0.2,0.2,0.2); // ok repond CELL
  // Point q(.5,.5,0); // ok repond edge
  // Point q(0.2,0.2,0); // ok repond facet
  // Point q(0,0,0); // ok repond vertex
  // Point q(1,1,1); // ok repond OUTSIDE_CONVEX_HULL
  // Point q(1,1,0); // ok repond OUTSIDE_CONVEX_HULL

  // test dim 2
//   T.insert_outside_affine_hull(p0);
//   T.insert_outside_affine_hull(px);
//   T.insert_outside_affine_hull(py);

  // Point q(0.25,0.25,0); // ok facet
  // Point q(0,0,0); // ok vertex
  // Point q(0,0.25,0); // ok edge
  // Point q(1,1,0); // ok outside_convex_hull
  // Point q(1,1,1); // ok outside_affine_hull

  // test dim 1
//   Triangulation_3 T;
//   T.insert_outside_affine_hull(p0);
//   T.insert_outside_affine_hull(px);

  // Point q(.5,0,0); // ok edge
  // Point q(0,1,0);// ok outside_affine_hull
  // Point q(1,0,0);// ok vertex
  // Point q(2,0,0); // ok outside_convex_hull

  // test dim 0
  //Triangulation_3 T;
//    T.insert_outside_affine_hull(p0);
  
  // Point q(0,1,0);// ok outside_affine_hull
  // Point q(p0);// ok vertex

  //  TDS Ds = T.Triangulation_data_structure_3();
  //  affiche_sommets(T);
  //  affiche_sommets(Ds);

  //  affiche_aretes(T);
  // affiche_faces(T);
  //  affiche_cellules(T);
  //  affiche_cellules_voisins(T);

  //  visu_cellules( gv, T );
  // visu_faces( gv, T );
  // visu_aretes(gv,T);
  //  visu_sommets(gv,T);

  //  cout << "validite " << T.is_valid(true) << endl;

  // test locate
//   gv << q;
//   Locate_type lt;
//   int i,j;
//   Cell_handle c = T.locate(q,lt,i,j);
//   switch ( lt ) {
//   case 3: // CELL
//     {
//       cout << q << " CELL" << endl;
//       gv.set_face_color(CGAL_GREEN);
//       gv << T.tetrahedron(c);
//       break;
//     }
//   case 2: // FACET
//     {
//       cout << q << " FACET" << endl;
//       gv.set_face_color(CGAL_GREEN);
//       gv << T.triangle(c,i);
//       break;
//     }
//   case 1: // EDGE
//     {
//       cout << q << " EDGE" << endl;
//       gv << T.segment(c,i,j);
//       break;
//     }
//   case 0: // VERTEX
//     {
//       cout << q << " VERTEX" << endl;
//       gv.set_vertex_color(CGAL_GREEN);
//       gv << c->vertex(i)->point();
//        break;
//     }
//   case 4: // OUTSIDE_CONVEX_HULL
//     { 
//       cout << q << " OUTSIDE_CONVEX_HULL" << endl;
//       switch (T.dimension()) {
//       case 0:
// 	// impossible
// 	break;
//       case 1:
// 	{
// 	  gv.set_vertex_color(CGAL_GREEN);
// 	  gv << c->vertex(i)->point();
// 	  break;
// 	}
//       case 2:
// 	{
// 	  gv.set_edge_color(CGAL_GREEN);
// 	  gv << T.segment(c,i,j);
// 	  break;
// 	}
//       case 3:
// 	{
// 	  gv.set_face_color(CGAL_GREEN);
// 	  gv << T.triangle(c,i);
//       	  break;
// 	}
//       }
//       break;
//     }
//   case 5:
//     {
//       cout << q << " OUTSIDE_AFFINE_HULL" << endl;
//     }
//   }
  
//   char entree;
//   ifstream is(entree, os::in, filebuf::openprot);

//   if( ! is )
//     {
//       cerr << "unable to open " << opt.finname << " for input" << endl;
//       return false;
//     }

//   CGAL_set_ascii_mode(is);
//   is >> nouv;
  
//   Triangulation_3 T2;
//   int x,y,z=0;
//   for (z=0 ; z<2 ; z++)
//     for (y=0 ; y<5 ; y++)
//       for (x=0 ; x<5 ; x++) {
// 	T2.insert(Point(x,y,z));
// 	assert(T2.is_valid(true));
//       }
//   //   assert(T2.number_of_vertices()==125);
//   //   assert(T2.dimension()==3);
  
//   affiche_aretes_plus(T2);
//   affiche_aretes_minus(T2);
  
//   affiche_faces_plus(T2);
//   affiche_faces_minus(T2);
   
  Point nouv;

  cout << "point ? " << endl;
  cin >> nouv ;

  while ( nouv.x() != 3000 ) {
    insere(gv,T,nouv);
    gv.clear();
    switch (T.dimension()) {
    case 0:
      {
	visu_sommets(gv,T);
	break;
      }
    case 1:
      {
	visu_aretes(gv,T);
	break;
      }
    case 2:
      {
	visu_faces(gv,T);
	break;
      }
    case 3:
      {
	visu_cellules(gv,T);
	break;
      }
    }
    cout << "point ? " << endl;
    cin >> nouv ;
  }

//   cout << "copy constructeur" << endl;
//   cout << "validite T " << T.is_valid(true) << endl;
//   Delaunay T1(T);
//   cout << "validite copie " << T1.is_valid(true) << endl;
//   cout << "assert T valid" << endl;
//   assert(T.is_valid(true));
//   cout << "assert T1 valid" << endl;
//   assert(T1.is_valid(true));


//   Cell_handle c = T2.infinite_cell();
//   Edge e = CGAL_make_triple(c,0,1);
//   Cell_circulator ccir = T2.incident_cells(e,c);

// //   TDS tds1(T2.tds());
// //   TDS tds2;
// //   tds1.clear();
// //   tds2.clear();

//   TDSVertex* v1 = new TDSVertex(Point(1,0,0));
//   TDSVertex* v2 = new TDSVertex(Point(2,0,0));
//   TDSVertex* v3 = new TDSVertex(Point(3,0,0));
//   TDSVertex* v4 = new TDSVertex(Point(4,0,0));
//   TDSVertex* v5 = new TDSVertex(Point(5,0,0));

//   cout << "    insertions" << endl;
//   TDS tds0;
//   tds0.insert_outside_affine_hull(v1, NULL);
//   cout << "avec v1" << endl;
//   affiche_sommets(tds0);
//   affiche_validite(tds0);
//   tds0.insert_outside_affine_hull(v2,v1);
//   cout << "avec v2" << endl;
//   affiche_sommets(tds0);
//   affiche_validite(tds0);
//   tds0.insert_outside_affine_hull(v3,v2);
//   cout << "avec v3" << endl;
//   affiche_sommets(tds0);
//   affiche_validite(tds0);
//   tds0.insert_outside_affine_hull(v4,v3);
//   cout << "avec v4" << endl;
//   affiche_sommets(tds0);
//   affiche_validite(tds0);
//   tds0.insert_outside_affine_hull(v5,v4);
//   cout << "avec v5" << endl;
//   affiche_sommets(tds0);
//   affiche_validite(tds0);

//   v1 = new TDSVertex(Point(10,0,0));
//   v2 = new TDSVertex(Point(20,0,0));
//   v3 = new TDSVertex(Point(30,0,0));
//   v4 = new TDSVertex(Point(40,0,0));
//   v5 = new TDSVertex(Point(50,0,0));

//   cout << "    constructors" << endl;
//   TDS tds1;
//   TDS tds2;
// //   cout << "tds1" << endl;
// //   affiche_sommets(tds1);
//   cout << "tds1" << endl;
//   affiche_validite(tds1);

//   tds2.insert_outside_affine_hull(v1, NULL);
//   cout << "tds2" << endl;
//   affiche_validite(tds2);
//   affiche_sommets(tds2);

//   TDS tds3(tds2);
//   cout << "tds3 avant" << endl;
//   affiche_validite(tds3);
//   affiche_sommets(tds3);
//   tds3.insert_outside_affine_hull(v2,v1);
//   cout << "tds3 apres" << endl;
//   affiche_validite(tds3);
//   affiche_sommets(tds3);

//   tds3.insert_outside_affine_hull(v3,v2);
//   cout << "tds3 + v3"<< endl;
//   affiche_validite(tds3);
//   affiche_sommets(tds3);

//   TDS tds4 = tds3;
//   cout << "tds4 avant" << endl;
//   affiche_validite(tds4);
//   affiche_sommets(tds4);
//   tds4.insert_outside_affine_hull(v3,v2);
//   cout << "tds4 apres" << endl;
//   affiche_validite(tds4);
//   affiche_sommets(tds4);

//   TDS tds5;
//   tds5.swap(tds4);
//   cout << "tds5 avant" << endl;
//   affiche_validite(tds5);
//   affiche_sommets(tds5);
//   tds5.insert_outside_affine_hull(v4,v3);
//   cout << "tds5 apres" << endl;
//   affiche_validite(tds5);
//   affiche_sommets(tds5);
  
//   TDS tds6;
//   tds6.swap(tds5);
//   cout << "tds6 avant" << endl;
//   affiche_validite(tds6);
//   affiche_sommets(tds6);
//   tds6.insert_outside_affine_hull(v5,v4);
//   cout << "tds6 apres" << endl;
//   affiche_validite(tds6);
//   affiche_sommets(tds6);

// // Test constructors
//   cout << "    constructors" << endl;
//   TDS tds1;
//   TDS tds2;

//   TDSVertex_iterator vit;

//   tds2.insert_outside_affine_hull(v1, NULL);
//   cout << "tds2" << endl;
//   affiche_validite(tds2);
//   affiche_sommets(tds2);

//   TDS tds3(tds2);
//   cout << "tds3 avant" << endl;
//   affiche_validite(tds3);
//   affiche_sommets(tds3);

//   vit=tds3.vertices_begin();
//   tds3.insert_outside_affine_hull(v2,&*vit);
//   cout << "tds3 apres" << endl;
//   affiche_validite(tds3);
//   affiche_sommets(tds3);

//   TDS tds4 = tds3;
//   cout << "tds4 avant" << endl;
//   affiche_validite(tds4);
//   affiche_sommets(tds4);

//   vit=tds4.vertices_begin();
//   tds4.insert_outside_affine_hull(v3,&*vit);
//   cout << "tds4 apres" << endl;
//   affiche_validite(tds4);
//   affiche_sommets(tds4);
  
//   TDS tds5;
//   tds5.swap(tds4);
//   cout << "tds5 avant" << endl;
//   affiche_validite(tds5);
//   affiche_sommets(tds5);

//   vit=tds5.vertices_begin();
//   tds5.insert_outside_affine_hull(v4,&*vit);
//   cout << "tds5 apres" << endl;
//   affiche_validite(tds5);
//   affiche_sommets(tds5);

//   TDS tds6;
//   tds6.swap(tds5);
//   cout << "tds6 avant" << endl;
//   affiche_validite(tds6);
//   affiche_sommets(tds6);

//   vit=tds6.vertices_begin();
//   tds6.insert_outside_affine_hull(v5,&*vit);
//   cout << "tds6 apres" << endl;
//   affiche_validite(tds6);
//   affiche_sommets(tds6);

  char ch;
  cout << "donner caractere de fin" << endl;
  cin >> ch;

  return 1;
}
