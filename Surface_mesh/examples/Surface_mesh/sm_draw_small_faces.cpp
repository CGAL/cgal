#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <fstream>
#include <string>
#include "draw_surface_mesh_small_faces.h"

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef Mesh::Vertex_index             vertex_descriptor;
typedef Mesh::Face_index               face_descriptor;
typedef K::FT                          FT;

int main(int argc, char* argv[])
{
  Mesh sm;
  std::ifstream input((argc>1)?argv[1]:"data/elephant.off");
  input>>sm;

  CGAL::Polygon_mesh_processing::triangulate_faces(sm);

  Mesh::Property_map<face_descriptor, FT> faces_size;
  bool created;
  boost::tie(faces_size, created)=sm.add_property_map<face_descriptor, FT>("f:size",0.);
  CGAL_assertion(created);

  for(face_descriptor fd : sm.faces())
  { faces_size[fd]=CGAL::Polygon_mesh_processing::face_area(fd, sm); }

  draw_surface_mesh_with_small_faces(sm);

  sm.remove_property_map(faces_size);

  return EXIT_SUCCESS;
}

