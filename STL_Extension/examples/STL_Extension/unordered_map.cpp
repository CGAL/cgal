
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <boost/unordered_map.hpp>


typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point_3;
typedef CGAL::Surface_mesh<Point_3>     Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

int main()
{
  boost::unordered_map<vertex_descriptor, int> bum;
  Mesh mesh;
  vertex_descriptor vd  = mesh.add_vertex();
  bum[vd] = 7812;
  return 0;
}
