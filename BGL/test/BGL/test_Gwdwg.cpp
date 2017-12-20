#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Graph_with_descriptor_with_graph.h>
#include <fstream>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> SM;
typedef CGAL::Graph_with_descriptor_with_graph<SM> Mesh;

int main()
{
  SM sm;
  Mesh mesh(sm);
  std::ifstream in("data/cube.off");
  in >> sm;

  assert( num_vertices(mesh) == num_vertices(sm) );

  SM sm2;
  sm2 = sm;
  Mesh mesh2(sm2);
  try {
    if( target( *(halfedges(mesh).first), mesh2) == *(vertices(mesh).first)){
      CGAL_error_msg("The previous lie should have throw a exception");
    }
  } catch(...){
    std::cerr << "we caught it" << std::endl;
  }
    std::cout << "done" << std::endl;
  return 0;
}

