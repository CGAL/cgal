// Triangulation_3/example_hierarchy.C
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>

#include <cassert>
#include <vector>

typedef CGAL::Filtered_kernel<CGAL::Simple_cartesian<double> > my_K;

// This is just to shorten some symbol names for VC++
struct K : public my_K {};

typedef CGAL::Triangulation_vertex_base_3<K>             Vb;
typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb>  Vbh;
typedef CGAL::Triangulation_cell_base_3<void>            Cb;
typedef CGAL::Triangulation_data_structure_3<Vbh,Cb>     Tds;
typedef CGAL::Delaunay_triangulation_3<K,Tds>            Dt;
typedef CGAL::Triangulation_hierarchy_3<Dt>              Dh;

typedef Dh::Vertex_iterator Vertex_iterator;
typedef Dh::Vertex_handle Vertex_handle;
typedef K::Point_3 Point;

int main()
{
  Dh T;

  // insertion of points on a 3D grid
  int x,y,z;
  std::vector<Vertex_handle> V(125);
  int i=0;
  
  for (z=0 ; z<5 ; z++)
    for (y=0 ; y<5 ; y++)
      for (x=0 ; x<5 ; x++) 
	  V[i++] = T.insert(Point(x,y,z));

  assert( T.is_valid() );
  assert( T.number_of_vertices() == 125 );
  assert( T.dimension() == 3 );

  // removal of the vertices in random order
  std::random_shuffle(V.begin(), V.end());

  for (i=0; i<125; ++i)
    assert( T.remove(V[i]) );

  assert( T.is_valid() );
  assert( T.number_of_vertices() == 0 );

  return 0;
}
