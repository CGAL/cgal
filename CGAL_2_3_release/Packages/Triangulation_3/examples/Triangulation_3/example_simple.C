// Triangulation_3/example_simple.C
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

#include <CGAL/Triangulation_3.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>

typedef CGAL::Filtered_kernel<CGAL::Simple_cartesian<double> > my_K;

// This is just to shorten some symbol names for VC++
struct K : public my_K {};

typedef CGAL::Triangulation_3<K> Triangulation;

typedef Triangulation::Cell_handle Cell_handle;
typedef Triangulation::Vertex_handle Vertex_handle;
typedef Triangulation::Locate_type Locate_type;

typedef K::Point_3 Point;

int main()
{
  Triangulation T;

  // insertion from a list :
  std::list<Point> L;
  L.push_front(Point(0,0,0));
  L.push_front(Point(1,0,0));
  L.push_front(Point(0,1,0));

  int n = T.insert(L.begin(), L.end());

  // insertion from a vector :
  std::vector<Point> V(3);
  V[0] = Point(0,0,1);
  V[1] = Point(1,1,1);
  V[2] = Point(2,2,2);

  n = n + T.insert(V.begin(), V.end());

  // 6 points have been inserted :
  assert( n == 6 );

  // checking validity of T :
  assert( T.is_valid(false) );

  Locate_type lt;
  int li, lj;
  Point p(0,0,0);
  Cell_handle c = T.locate(p, lt, li, lj);
  // p is the vertex of c of index li :
  assert( lt == Triangulation::VERTEX );
  assert(  c->vertex(li)->point() == p );

  Vertex_handle v = c->vertex( (li+1)&3 );
  // v is another vertex of c
  Cell_handle nc = c->neighbor(li);
  // nc = neighbor of c opposite to the vertex associated with p
  // nc must have vertex v :
  int nli;
  assert(  nc->has_vertex( v, nli ) );
  // nli is the index of v in nc

  std::ofstream oFileT("output",std::ios::out);
  // writing file output; 
  oFileT << T; 

  Triangulation T1;
  std::ifstream iFileT("output",std::ios::in);
  // reading file output; 
  iFileT >> T1; 
  assert( T1.is_valid() );
  assert( T1.number_of_vertices() == T.number_of_vertices() );
  assert( T1.number_of_cells() == T.number_of_cells() );

  return 0;
}
