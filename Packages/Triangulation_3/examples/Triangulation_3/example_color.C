#include <CGAL/Cartesian.h>

#include <CGAL/Delaunay_triangulation_3.h>

#include <cassert>

template < class Gt >
class My_vertex_base
  : public CGAL::Triangulation_vertex_base_3<Gt>
{
public :
  CGAL::Color color;

  My_vertex_base() 
    : CGAL::Triangulation_vertex_base_3<Gt>() 
    {color=CGAL::WHITE;}

  My_vertex_base(const Point & p) 
    : CGAL::Triangulation_vertex_base_3<Gt>(p) 
    {color=CGAL::WHITE;}

  My_vertex_base(const Point & p, void* c) 
    : CGAL::Triangulation_vertex_base_3<Gt>(p,c) 
    {color=CGAL::WHITE;}

  My_vertex_base(void* c)
    : CGAL::Triangulation_vertex_base_3<Gt>(c) 
    {color=CGAL::WHITE;}
};

typedef CGAL::Cartesian<double> Gt;
typedef Gt::Point_3 Point;

typedef CGAL::Triangulation_cell_base_3<Gt> Cb;
typedef My_vertex_base<Gt> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb,Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<Gt, Tds> Delaunay;

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
