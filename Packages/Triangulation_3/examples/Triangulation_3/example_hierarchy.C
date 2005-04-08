// file: examples/Triangulation_3/example_hierarchy.C

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

typedef Dh::Vertex_iterator Vertex_iterator;
typedef Dh::Vertex_handle   Vertex_handle;
typedef Dh::Point           Point;

int main()
{
  Dh T;

  // insertion of points on a 3D grid
  std::vector<Vertex_handle> V;
  
  for (int z=0 ; z<5 ; z++)
    for (int y=0 ; y<5 ; y++)
      for (int x=0 ; x<5 ; x++) 
	  V.push_back(T.insert(Point(x,y,z)));

  assert( T.is_valid() );
  assert( T.number_of_vertices() == 125 );
  assert( T.dimension() == 3 );

  // removal of the vertices in random order
  std::random_shuffle(V.begin(), V.end());

  for (int i=0; i<125; ++i)
    T.remove(V[i]);

  assert( T.is_valid() );
  assert( T.number_of_vertices() == 0 );

  return 0;
}
