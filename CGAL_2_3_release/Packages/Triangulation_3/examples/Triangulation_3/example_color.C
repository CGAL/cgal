// Triangulation_3/example_color.C
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <cassert>

template < class Traits >
class My_vertex_base
  : public CGAL::Triangulation_vertex_base_3<Traits>
{
public :
  CGAL::Color color;
  typedef typename Traits::Point_3 Point;

  My_vertex_base() 
    : CGAL::Triangulation_vertex_base_3<Traits>(), color(CGAL::WHITE)
    {}
 
  My_vertex_base(const Point & p) 
    : CGAL::Triangulation_vertex_base_3<Traits>(p), color(CGAL::WHITE)
    {}

  My_vertex_base(const Point & p, void* c) 
    : CGAL::Triangulation_vertex_base_3<Traits>(p,c), color(CGAL::WHITE)
    {}

  My_vertex_base(void* c)
    : CGAL::Triangulation_vertex_base_3<Traits>(c), color(CGAL::WHITE)
    {}
};

typedef CGAL::Filtered_kernel<CGAL::Simple_cartesian<double> > my_K;

// This is just to shorten some symbol names for VC++
struct K : public my_K {};

typedef K::Point_3 Point;

typedef CGAL::Triangulation_cell_base_3<K> Cb;
typedef My_vertex_base<K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb,Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;

typedef Delaunay::Vertex_iterator Vertex_iterator;
typedef Delaunay::Vertex_handle Vertex_handle;

int main()
{
  Delaunay T;

  T.insert(Point(0,0,0));
  T.insert(Point(1,0,0));  
  T.insert(Point(0,1,0));  
  T.insert(Point(0,0,1));  
  T.insert(Point(2,2,2));  
  T.insert(Point(-1,0,1));  

  Vertex_iterator vit;
  std::set<Vertex_handle> adjacent;
  for (vit = T.finite_vertices_begin(); vit != T.vertices_end(); ++vit) {
    T.incident_vertices( &*vit, adjacent);
    if (adjacent.size() == 6) 
      vit->color = CGAL::RED;
  }

  return 0;
}
