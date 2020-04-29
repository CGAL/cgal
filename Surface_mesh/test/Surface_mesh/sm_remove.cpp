
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

  assert(m.num_vertices() == 0);
  assert(m.number_of_removed_vertices() == 0);
  for(int i=0; i < 10; i++){
    u = m.add_vertex(Point_3(0,0,0));
    m.remove_vertex(u);
  }
  assert(m.num_vertices() == 1);
  assert(m.number_of_removed_vertices() == 1);


  assert(m.does_recycle_garbage());
  m.set_recycle_garbage(false);
  assert(! m.does_recycle_garbage());

  m.add_vertex(Point_3(0,0,0));
  assert(m.num_vertices() == 2);
  assert(m.number_of_removed_vertices() == 1);
  
  m.set_recycle_garbage(true);
  m.add_vertex(Point_3(0,0,0));
  assert(m.num_vertices() == 2);
  assert(m.number_of_removed_vertices() == 0);
  
  return 0;
}
