// Triangulation_3/example_color.C
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>

template < class Traits >
class My_vertex_base
  : public CGAL::Triangulation_vertex_base_3<Traits>
{
public :
  CGAL::Color color;

  My_vertex_base() 
    : CGAL::Triangulation_vertex_base_3<Traits>(), color(CGAL::WHITE)
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

int main()
{
  Delaunay T;

  T.insert(Point(0,0,0));
  T.insert(Point(1,0,0));  
  T.insert(Point(0,1,0));  
  T.insert(Point(0,0,1));  
  T.insert(Point(2,2,2));  
  T.insert(Point(-1,0,1));  

  Delaunay::Finite_vertices_iterator vit;
  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit) {
    std::vector<Delaunay::Vertex_handle> adjacent;
    T.incident_vertices( vit->handle(), std::back_inserter(adjacent));
    if (adjacent.size() == 6)
      vit->color = CGAL::RED;
  }

  return 0;
}
