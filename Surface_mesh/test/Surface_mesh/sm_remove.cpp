
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <iostream>
typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Sm;


int main()
{
  Sm m;
  Sm::vertex_index u;
  std::cout << m.num_vertices() << "  " << m.number_of_removed_vertices() << std::endl;
  for(int i=0; i < 10; i++){
    u = m.add_vertex(Point_3(0,0,0));
    m.remove_vertex(u);
  }

  std::cout << m.num_vertices() << "  " << m.number_of_removed_vertices() << std::endl;
  return 0;
}
