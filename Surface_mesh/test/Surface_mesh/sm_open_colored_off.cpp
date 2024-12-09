#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> SMesh;
typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;


void OpenOFF(int i)
{
  std::string path;
  switch(i)
  {
   case 1:
     path = CGAL::data_file_path("meshes/mesh_with_colors.off");
     break;
  case 2:
    path = "test2.off";
    break;
  case 3:
    path = "test3.off";
    break;
  }
  std::ifstream in(path.c_str());
  SMesh surface_mesh;
  in >> surface_mesh;
  assert(in && !surface_mesh.is_empty());


  SMesh::Property_map<face_descriptor, CGAL::IO::Color> fcolors =
      surface_mesh.property_map<face_descriptor, CGAL::IO::Color >("f:color").value();

  SMesh::Property_map<vertex_descriptor, CGAL::IO::Color> vcolors =
    surface_mesh.property_map<vertex_descriptor, CGAL::IO::Color >("v:color").value();

  CGAL::IO::Color c = fcolors[*(surface_mesh.faces().begin())];
  assert(c== CGAL::IO::Color(229,0,0));
  c = fcolors[*(--surface_mesh.faces().end())];
  assert(c== CGAL::IO::Color(0,0,229));

  c = vcolors[*(surface_mesh.vertices().begin())];
  assert((c== CGAL::IO::Color(229,0,0)));
  c = vcolors[*(--surface_mesh.vertices().end())];
  assert((c== CGAL::IO::Color(0,0,229)));
}


int main()
{
  OpenOFF(1);
  OpenOFF(2);
  OpenOFF(3);
  std::cerr << "done" << std::endl;
  return 0;
}
