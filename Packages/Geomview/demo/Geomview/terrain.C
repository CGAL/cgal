#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <unistd.h> // sleep

#include <CGAL/Gmpz.h>

#include <CGAL/Cartesian.h>

#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_euclidean_traits_xy_3.h>

#include <CGAL/Triangulation_data_structure_using_list_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_geom_traits_3.h>
#include <CGAL/Triangulation_data_structure_3.h>

#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_2.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>

typedef CGAL::Gmpz NT;
typedef CGAL::Cartesian<NT>  Rep;

typedef CGAL::Triangulation_euclidean_traits_2<Rep> Gt2;
typedef Gt2::Point Point2;
typedef CGAL::Triangulation_euclidean_traits_xy_3<Rep> Gt3;
typedef Gt3::Point Point3;

typedef CGAL::Triangulation_geom_traits_3<Rep> Gt3d;

typedef CGAL::Triangulation_vertex_base_2<Gt2> Vb2;
typedef CGAL::Triangulation_face_base_2<Gt2> Fb2;
typedef CGAL::Triangulation_vertex_base_2<Gt3> Vb3;
typedef CGAL::Triangulation_face_base_2<Gt3> Fb3;

typedef CGAL::Triangulation_vertex_base_3<Gt3d> Vb3d;
typedef CGAL::Triangulation_cell_base_3<Gt3d> Ce3d;

typedef CGAL::Triangulation_data_structure_using_list_2<Vb2,Fb2> Tds2;
typedef CGAL::Triangulation_data_structure_using_list_2<Vb3,Fb3> Tds3;

typedef CGAL::Triangulation_data_structure_3<Vb3d,Ce3d> Tds3d;

typedef CGAL::Delaunay_triangulation_2<Gt2, Tds2> Delaunay;
typedef CGAL::Delaunay_triangulation_2<Gt3, Tds3> Terrain;

typedef CGAL::Delaunay_triangulation_3<Gt3d, Tds3d> Delaunay3d;

void p_tetra()
{
  typedef CGAL::Cartesian<double>                               R;
  typedef CGAL::Halfedge_data_structure_polyhedron_default_3<R> HDS;
  typedef CGAL::Polyhedron_default_traits_3<R>                  Traits;
  typedef CGAL::Polyhedron_3<Traits,HDS>                        Polyhedron;
  typedef Polyhedron::Point                                     Point;
  typedef Polyhedron::Vertex_iterator                           Vertex_iterator;
  Point p( 1.0, 0.0, 0.0);
  Point q( 0.0, 1.0, 0.0);
  Point r( 0.0, 0.0, 1.0);
  Point s( 0.0, 0.0, 0.0);

  Polyhedron P;
  P.make_tetrahedron( p, q, r, s);
  CGAL::Geomview_stream gv(CGAL::Bbox_3(150,150,150, 350, 350, 350));
  gv.set_line_width(4);
  gv.set_bg_color(CGAL::Color(0, 200, 200));
  gv.set_trace(true);
  gv << P;
  sleep(10);
}

int main()
{
  // p_tetra();
  CGAL::Geomview_stream gv(CGAL::Bbox_3(150,150,150, 350, 350, 350));
  gv.set_line_width(4);
  gv.set_trace(false);
  gv.set_bg_color(CGAL::Color(0, 200, 200));

  Delaunay D;
  Delaunay3d D3d;
  Terrain T;
  Point3 p;
  std::ifstream iFile("points3",std::ios::in);
  while ( iFile >> p ) 
    { 
      D.insert( Point2(p.x(),p.y()) ); 
      D3d.insert( p );
      T.insert( p );
    }

#if 1
  gv << D3d;
  CGAL::show_triangulation_edges(gv,D3d);
#else

#if 1
  // std::cout << D << std::endl; return 0;
  // gv.set_trace(true);
  // gv << D;
  CGAL::show_triangulation_edges(gv,D);
  // gv.set_trace(false);
#else
  Delaunay::Finite_edges_iterator dit = D.finite_edges_begin();
  for ( ; dit != D.edges_end() ; ++dit ) {
    gv << D.segment( *dit ) ;
  }
#endif

#if 1
  gv << T;
#else
  Terrain::Finite_faces_iterator tit = T.finite_faces_begin();
  for ( ; tit != T.faces_end() ; ++tit ) {
    gv << CGAL::Color(200, 0, 200);
    gv << T.triangle( &*tit ) ;
  }
#endif
#endif

  gv.set_trace(true);
  std::cout << "entrez un point" << std::endl;
  // char tmpbuf[1024];
  // gv << CGAL::ascii << "(echo (geomview-version))";
  // gv >> tmpbuf;
  // std::cout << tmpbuf << std::endl;
  // gv << "(pickable pickplane yes) (ui-target pickplane yes)(interest (pick world pickplane * nil nil nil nil nil nil nil))" ;
  while(1) {
    gv >> p;
    std::cout << p << std::endl;
  }

  std::cout << "fini" << std::endl;

  char ch;
  std::cin >> ch;

  return 0;
}
