#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>

#include <cassert>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_3<K>             Vb;
typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb>  Vbh;
typedef CGAL::Triangulation_data_structure_3<Vbh>        Tds;
typedef CGAL::Delaunay_triangulation_3<K,Tds>            Dt;
typedef CGAL::Triangulation_hierarchy_3<Dt>              Dh;

typedef Dh::Finite_vertices_iterator Finite_vertices_iterator;
typedef Dh::Vertex_handle            Vertex_handle;
typedef Dh::Point                    Point;

int main()
{
  // generating points on a grid.
  std::vector<Point> P;

  for (int z=0 ; z<5 ; z++)
    for (int y=0 ; y<5 ; y++)
      for (int x=0 ; x<5 ; x++)
	  P.push_back(Point(x,y,z));

  Dh T;

  // using the range insert (it is faster than inserting points one by one)
  T.insert (P.begin(), P.end());

  assert( T.is_valid() );
  assert( T.number_of_vertices() == 125 );
  assert( T.dimension() == 3 );

  // get the vertices.
  std::vector<Vertex_handle> V;
  for (Finite_vertices_iterator v = T.finite_vertices_begin();
          v != T.finite_vertices_end(); ++v)
      V.push_back (v);

  // removal of the vertices in random order
  std::random_shuffle(V.begin(), V.end());

  for (int i=0; i<125; ++i)
    T.remove(V[i]);

  assert( T.is_valid() );
  assert( T.number_of_vertices() == 0 );

  return 0;
}
