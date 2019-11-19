#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/Graph_with_descriptor_with_graph.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/Euler_operations.h>

#include <fstream>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> SM;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Graph_with_descriptor_with_graph<SM> Mesh;
typedef CGAL::Graph_with_descriptor_with_graph<Polyhedron> PMesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<PMesh>::vertex_descriptor pvertex_descriptor;
int main()
{
  SM sm;
  Mesh mesh(sm);

  CGAL::make_hexahedron(
        Point_3(-1,-1,-1),
        Point_3(1,-1,-1),
        Point_3(1,1,-1),
        Point_3(-1,1,-1),
        Point_3(-1,1,1),
        Point_3(-1,-1,1),
        Point_3(1,-1,1),
        Point_3(1,1,1),
        sm
        );


  vertex_descriptor vd = * vertices(mesh).first;
  std::cout << "Mesh is aware that it is based on a "<<typeid(*vd.graph).name()<<"."<< std::endl;


  Polyhedron poly;
  PMesh pmesh(poly);

  CGAL::make_hexahedron(
        Point_3(-0.5,-0.5,-0.5),
        Point_3(0.5,-0.5,-0.5),
        Point_3(0.5,0.5,-0.5),
        Point_3(-0.5,0.5,-0.5),
        Point_3(-0.5,0.5,0.5),
        Point_3(-0.5,-0.5,0.5),
        Point_3(0.5,-0.5,0.5),
        Point_3(0.5,0.5,0.5),
        poly
        );
  pvertex_descriptor pvd = *  vertices(pmesh).first;
  std::cout << "Pmesh is aware that it is based on a "<<typeid(*pvd.graph).name()<<"."<< std::endl;

  return 0;
}

