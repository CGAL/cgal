// file          : demo/Triangulation3/demo_color.C
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>

// Geomview doesn't work on M$ at the moment, so we don't compile this file.
#if defined(__BORLANDC__) || defined(_MSC_VER)
#include <iostream>
int main()
{
  std::cerr << "Geomview doesn't work on Windows, so this demo doesn't work"
            << std::endl;
  return 0;
}
#else

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>

struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};

typedef CGAL::Triangulation_vertex_base_with_info_3<CGAL::Color, K>  Vb;
typedef CGAL::Triangulation_data_structure_3<Vb>                     Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>                       Delaunay;

typedef Delaunay::Point Point;

int main()
{
  CGAL::Geomview_stream gv(CGAL::Bbox_3(0,0,0, 2, 2, 2));
  gv.set_bg_color(CGAL::Color(0, 200, 200));
  gv.clear();

  Delaunay T;

  T.insert(Point(0,0,0));
  T.insert(Point(1,0,0));  
  T.insert(Point(0,1,0));  
  T.insert(Point(0,0,1));  
  T.insert(Point(2,2,2));  
  T.insert(Point(-1,0,1));  

  // Set the color of finite vertices of degree 6 to red.
  Delaunay::Finite_vertices_iterator vit;
  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
    if (T.degree(vit) == 6)
      vit->info() = CGAL::RED;

  std::cout << "           Visualization of T" << std::endl;
  gv.set_wired(true);
  gv << T;

  std::cout << "           Vertices of T with their own color" << std::endl
	    << "           red for degree 6 (counting infinite vertex)" 
	    << std::endl 
	    << "           white otherwise" << std::endl;
  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
    gv << vit->info() << vit->point();

  std::cout << "Enter any character to quit" << std::endl;
  char ch;
  std::cin >> ch;

  return 0;
}

#endif // if defined(__BORLANDC__) || defined(_MSC_VER)
