#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <unistd.h> // sleep

#include <CGAL/Gmpz.h>

#include <CGAL/Cartesian.h>

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

int main()
{
  CGAL::Geomview_stream gv(CGAL::Bbox_3(150, 150, 150, 350, 350, 350));
  gv.set_line_width(4);
  gv.set_trace(false);
  gv.set_bg_color(CGAL::Color(0, 200, 200));

  Delaunay D;
  Delaunay3d D3d;
  Terrain T;
  std::ifstream iFile("data/points3", std::ios::in);
  Point3 p;

  while ( iFile >> p ) 
  { 
      D.insert( Point2(p.x(), p.y()) ); 
      D3d.insert( p );
      T.insert( p );
  }

  // use different colors, and put a few sleeps/clear.

  gv << CGAL::BLUE;
  std::cout << "Drawing 2D Delaunay triangulation in wired mode.\n";
  gv.set_wired(true);
  gv << D;
  sleep(5);
  std::cout << "Drawing 2D Delaunay triangulation in non-wired mode.\n";
  gv.set_wired(false);
  gv << D;
  sleep(5);

  std::cout << "Drawing 3D Delaunay triangulation in wired mode.\n";
  gv.set_wired(true);
  gv << D3d;
  sleep(5);
  std::cout << "Drawing 3D Delaunay triangulation in non-wired mode.\n";
  gv.set_wired(false);
  gv << D3d;
  sleep(5);

  std::cout << "Drawing Terrain in wired mode.\n";
  gv.set_wired(true);
  gv << T;
  sleep(5);
  std::cout << "Drawing Terrain in non-wired mode.\n";
  gv.set_wired(false);
  gv << T;

  std::cout << "Enter a key to finish" << std::endl;
  char ch;
  std::cin >> ch;

  return 0;
}
