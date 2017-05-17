#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Periodic_3_regular_triangulation_3.h>
#include <CGAL/periodic_3_triangulation_3_io.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Periodic_3_regular_triangulation_traits_3<K>        Gt;
typedef CGAL::Periodic_3_regular_triangulation_3<Gt>              P3RT3;

typedef P3RT3::Bare_point        Point;
typedef P3RT3::Weighted_point    Weighted_point;
typedef P3RT3::Iso_cuboid        Iso_cuboid;
typedef P3RT3::Vertex_handle     Vertex_handle;
typedef P3RT3::Cell_handle       Cell_handle;
typedef P3RT3::Locate_type       Locate_type;

int main(int, char**)
{
  Iso_cuboid domain(-1,-1,-1, 2,2,2);  // the cube for the periodic domain

  // construction from a list of weighted points :
  std::list<Weighted_point> L;
  L.push_front(Weighted_point(Point(0,0,0), 0.01));
  L.push_front(Weighted_point(Point(1,0,0), 0.02));
  L.push_front(Weighted_point(Point(0,1,0), 0.03));

  P3RT3 T(L.begin(), L.end(), domain); // put the domain with the constructor

  P3RT3::size_type n = T.number_of_vertices();

  // insertion from a vector :
  std::vector<Weighted_point> V(3);
  V[0] = Weighted_point(Point(0,0,1), 0.04);
  V[1] = Weighted_point(Point(1,1,1), 0.05);
  V[2] = Weighted_point(Point(-1,-1,-1), 0.06);

  n = n + T.insert(V.begin(), V.end());

  assert( n == 6 );       // 6 points have been inserted
  assert( T.is_valid() ); // checking validity of T

  Locate_type lt;
  int li, lj;
  Weighted_point p(Point(0,0,0), 1.);
  Cell_handle c = T.locate(p, lt, li, lj);
  // p is the vertex of c of index li :
  assert( lt == P3RT3::VERTEX );
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
  std::ofstream oFileT("output_regular.tri", std::ios::out); // as a .tri file
  oFileT << T;

  std::ofstream to_off("output_regular.off");
  write_triangulation_to_off(to_off, T);

  std::ofstream d_to_off("output_regular_dual.off"); // as a .off file
  draw_dual_to_off(d_to_off, T);

  // reading file output
  P3RT3 T1;
  std::ifstream iFileT("output_regular.tri",std::ios::in);
  iFileT >> T1;
  assert( T1.is_valid() );
  assert( T1.number_of_vertices() == T.number_of_vertices() );
  assert( T1.number_of_cells() == T.number_of_cells() );

  return 0;
}
