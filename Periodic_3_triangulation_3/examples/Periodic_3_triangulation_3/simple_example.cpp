#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/periodic_3_triangulation_3_io.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K>       Gt;
typedef CGAL::Periodic_3_Delaunay_triangulation_3<Gt>             P3DT3;

typedef P3DT3::Point             Point;
typedef P3DT3::Iso_cuboid        Iso_cuboid;
typedef P3DT3::Vertex_handle     Vertex_handle;
typedef P3DT3::Cell_handle       Cell_handle;
typedef P3DT3::Locate_type       Locate_type;

int main(int, char**)
{
  Iso_cuboid domain(-1,-1,-1,2,2,2);  // the fundamental domain

  // construction from a list of points :
  std::list<Point> L;
  L.push_front(Point(0,0,0));
  L.push_front(Point(1,0,0));
  L.push_front(Point(0,1,0));

  P3DT3 T(L.begin(), L.end(), domain); // put the domain with the constructor

  P3DT3::size_type n = T.number_of_vertices();

  // insertion from a vector :
  std::vector<Point> V(3);
  V[0] = Point(0,0,1);
  V[1] = Point(1,1,1);
  V[2] = Point(-1,-1,-1);

  n = n + T.insert(V.begin(), V.end());

  assert( n == 6 );       // 6 points have been inserted
  assert( T.is_valid() ); // checking validity of T

  Locate_type lt;
  int li, lj;
  Point p(0,0,0);
  Cell_handle c = T.locate(p, lt, li, lj);
  // p is the vertex of c of index li :
  assert( lt == P3DT3::VERTEX );
  assert( c->vertex(li)->point() == p );

  Vertex_handle v = c->vertex( (li+1)&3 );
  // v is another vertex of c
  Cell_handle nc = c->neighbor(li);
  // nc = neighbor of c opposite to the vertex associated with p
  // nc must have vertex v :
  int nli;
  assert( nc->has_vertex( v, nli ) );
  // nli is the index of v in nc

  // writing file output
  std::ofstream oFileT("output.tri", std::ios::out); // as a .tri file
  oFileT << T;

  std::ofstream to_off("output_regular.off"); // as a .off file
  CGAL::write_triangulation_to_off(to_off, T);

  std::ofstream d_to_off("output_dual.off");
  draw_dual_to_off(d_to_off, T);

  // reading file output
  P3DT3 T1;
  std::ifstream iFileT("output.tri",std::ios::in);
  iFileT >> T1;
  assert( T1.is_valid() );
  assert( T1.number_of_vertices() == T.number_of_vertices() );
  assert( T1.number_of_cells() == T.number_of_cells() );

  return 0;
}
