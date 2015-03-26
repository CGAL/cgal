#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;

typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef CGAL::Surface_mesh<Point> Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;

void test(const char* file_name)
{
  Surface_mesh mesh;
  std::ifstream input(file_name);
  if (!(input >> mesh))
  {
    std::cerr << "Error: cannot read Surface_mesh : " << file_name << "\n";
    CGAL_assertion(false);
  }

  Surface_mesh::Property_map<face_descriptor, Vector> fnormals;
  bool created;
  boost::tie(fnormals, created) = mesh.add_property_map<face_descriptor,Vector>("f:normals",Vector(0,0,0));
  CGAL::Polygon_mesh_processing::compute_face_normals(mesh, fnormals);
  CGAL::Polygon_mesh_processing::compute_face_normals(mesh, fnormals,
    CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh.points()));
  CGAL::Polygon_mesh_processing::compute_face_normals(mesh, fnormals,
    CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh.points()).kernel(K()));

 Surface_mesh::Property_map<vertex_descriptor,Vector> vnormals;

  boost::tie(vnormals, created) = mesh.add_property_map<vertex_descriptor,Vector>("v:normals",Vector(0,0,0));
  CGAL::Polygon_mesh_processing::compute_vertex_normals(mesh, vnormals);
  CGAL::Polygon_mesh_processing::compute_vertex_normals(mesh, vnormals,
    CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh.points())); 
  CGAL::Polygon_mesh_processing::compute_vertex_normals(mesh, vnormals,
    CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh.points()).kernel(K()));

  CGAL::Polygon_mesh_processing::compute_normals(mesh, vnormals, fnormals);
  CGAL::Polygon_mesh_processing::compute_normals(mesh, vnormals, fnormals,
    CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh.points()));
  CGAL::Polygon_mesh_processing::compute_normals(mesh, vnormals, fnormals,
    CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh.points()).kernel(K()));

  BOOST_FOREACH(face_descriptor fd , faces(mesh)){
    std::cout << fnormals[fd] << std::endl;
  }
}

int main()
{
  test("data/elephant.off");

  std::cerr << "All done." << std::endl;
}
