// Triangulation_3/example_color.C
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>

template < class Traits, class Vb = CGAL::Triangulation_vertex_base_3<Traits> >
struct My_vertex_base
  : public Vb
{
  CGAL::Color color;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef My_vertex_base<Traits, Vb2>                    Other;
  };

  My_vertex_base() 
    : color(CGAL::WHITE) {}
};

typedef CGAL::Filtered_kernel<CGAL::Simple_cartesian<double> > my_K;

// This is just to shorten some symbol names for VC++
struct K : public my_K {};

typedef CGAL::Triangulation_data_structure_3<My_vertex_base<K> > Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>                   Delaunay;

typedef Delaunay::Point   Point;

int main()
{
  Delaunay T;

  T.insert(Point(0,0,0));
  T.insert(Point(1,0,0));  
  T.insert(Point(0,1,0));  
  T.insert(Point(0,0,1));  
  T.insert(Point(2,2,2));  
  T.insert(Point(-1,0,1));  

  // Set the color of finite vertices of degree 6 to red.
  Delaunay::Finite_vertices_iterator vit;
  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
    if (T.degree(vit) == 6)
      vit->color = CGAL::RED;

  return 0;
}
