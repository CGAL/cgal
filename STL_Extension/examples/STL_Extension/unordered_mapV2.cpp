
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Small_unordered_mapV2.h>
#include <CGAL/Surface_mesh.h>
#include <fstream>

typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point_3;
typedef CGAL::Surface_mesh<Point_3>     Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

int main()
{
  typedef CGAL::Small_unordered_mapV2<vertex_descriptor, int, boost::hash<vertex_descriptor>,8,6> Small;
  Small small;
  Mesh mesh;
  
  vertex_descriptor v0  = mesh.add_vertex();
  vertex_descriptor v1  = mesh.add_vertex();
  vertex_descriptor v2  = mesh.add_vertex();
  vertex_descriptor v3  = mesh.add_vertex();
  vertex_descriptor v4  = mesh.add_vertex();
  vertex_descriptor v5  = mesh.add_vertex();
  vertex_descriptor v6  = mesh.add_vertex();
  vertex_descriptor v7  = mesh.add_vertex();
  vertex_descriptor v8  = mesh.add_vertex();
  vertex_descriptor v9  = mesh.add_vertex();
  vertex_descriptor v10 = mesh.add_vertex();
  vertex_descriptor v11 = mesh.add_vertex();

  small.set(v4, 4);
  small.set(v11, 11);
  small.set(v2, 2);
  small.set(v9, 9);
  small.set(v8, 8);

  std::cout << "After insertion of 5 elements" << std::endl;
  for(const std::pair<vertex_descriptor,int>& vi : small){
    std::cout << vi.first << " " << vi.second << std::endl;
  }

  std::cout << "After insertion of 7 more elements" << std::endl;
  small.set(v3, 3);
  small.set(v6, 6);
  small.set(v5, 5);
  small.set(v0, 0);
  small.set(v10, 10);
  small.set(v7, 7);
  small.set(v1, 1);

  for(const std::pair<vertex_descriptor,int>& vi : small){
    std::cout << vi.first << " " << vi.second << std::endl;
  }

  std::cout << "Iterate over all vertices" << std::endl;
  for(vertex_descriptor vd : vertices(mesh)){
    std::cout << vd << " " << small.get(vd) << std::endl;
  }
  

  return 0;
}
