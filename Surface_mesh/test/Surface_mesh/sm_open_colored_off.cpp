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


  auto [fcolors, created_fcolors] = surface_mesh.add_property_map<face_descriptor, CGAL::IO::Color >("f:color");
  auto [vcolors, created_vcolors] = surface_mesh.add_property_map<vertex_descriptor, CGAL::IO::Color >("v:color");

  // Both color maps should have already existed, because they were loaded from the file
  assert(!created_fcolors);
  assert(!created_vcolors);

  auto first_fcolor = fcolors[*(surface_mesh.faces().begin())];
  assert(first_fcolor == CGAL::IO::Color(229, 0, 0));
  auto last_fcolor = fcolors[*(--surface_mesh.faces().end())];
  assert(last_fcolor == CGAL::IO::Color(0, 0, 229));

  auto first_vcolor = vcolors[*(surface_mesh.vertices().begin())];
  assert((first_vcolor == CGAL::IO::Color(229, 0, 0)));
  auto last_vcolor = vcolors[*(--surface_mesh.vertices().end())];
  assert((last_vcolor == CGAL::IO::Color(0, 0, 229)));
}


int main()
{
  OpenOFF(1);
  OpenOFF(2);
  OpenOFF(3);
  std::cout << "done" << std::endl;
  return 0;
}
