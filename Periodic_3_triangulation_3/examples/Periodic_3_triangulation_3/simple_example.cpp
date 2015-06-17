#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K> GT;

typedef CGAL::Periodic_3_Delaunay_triangulation_3<GT> PDT;

typedef PDT::Cell_handle       Cell_handle;
typedef PDT::Vertex_handle     Vertex_handle;
typedef PDT::Locate_type       Locate_type;
typedef PDT::Point             Point;
typedef PDT::Iso_cuboid        Iso_cuboid;

int main()
{
  Iso_cuboid domain(-1,-1,-1,2,2,2);  // The cube for the periodic domain

  // construction from a list of points :
  std::list<Point> L;
  L.push_front(Point(0,0,0));
  L.push_front(Point(1,0,0));
  L.push_front(Point(0,1,0));

  PDT T(L.begin(), L.end(), domain); // Put the domain with the constructor

  PDT::size_type n = T.number_of_vertices();

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
  assert( lt == PDT::VERTEX );
  assert( c->vertex(li)->point() == p );

  Vertex_handle v = c->vertex( (li+1)&3 );
  // v is another vertex of c
  Cell_handle nc = c->neighbor(li);
  // nc = neighbor of c opposite to the vertex associated with p
  // nc must have vertex v :
  int nli;
  assert( nc->has_vertex( v, nli ) );
  // nli is the index of v in nc

  std::ofstream oFileT("output.tri",std::ios::out);
  // writing file output; 
  oFileT << T; 

  PDT T1;
  std::ifstream iFileT("output.tri",std::ios::in);
  // reading file output; 
  iFileT >> T1; 
  assert( T1.is_valid() );
  assert( T1.number_of_vertices() == T.number_of_vertices() );
  assert( T1.number_of_cells() == T.number_of_cells() );

  return 0;
}
