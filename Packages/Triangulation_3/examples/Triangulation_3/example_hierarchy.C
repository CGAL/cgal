// Triangulation_3/example_hierarchy.C
#include <CGAL/Cartesian.h>

// Workaround for VC++ necessary for Filtered_exact.
#ifdef CGAL_CFG_MATCHING_BUG_2
#  define CGAL_IA_CT double
#  define CGAL_IA_PROTECTED true
#  define CGAL_IA_CACHE No_Filter_Cache
#  define CGAL_IA_ET CGAL::MP_Float
#endif

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Filtered_exact.h>

#include <cassert>
#include <vector>

// using Filtered_exact number type is advised :

typedef CGAL::Filtered_exact<double, CGAL::MP_Float> NT;

// chosing the representation (Cartesian or homogeneous)

typedef CGAL::Cartesian<NT> K;

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
