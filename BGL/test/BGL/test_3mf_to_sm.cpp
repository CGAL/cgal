#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include "Model/COM/NMR_DLLInterfaces.h"

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// Use NMR namespace for the interfaces
using namespace NMR;
namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef std::vector<Point_3> PointRange;
typedef std::vector<std::size_t> Polygon;
typedef std::vector<Polygon> PolygonRange;
typedef std::vector<CGAL::IO::Color> ColorRange;

int main(int argc, char** argv)
{
#ifdef CGAL_LINKED_WITH_3MF
  const char* filename = (argc == 2) ? argv[1] : "data/test.3mf";

  std::vector<PointRange> all_points;
  std::vector<PolygonRange> all_polygons;
  std::vector<ColorRange> all_colors;
  std::vector<std::string> names;
  std::vector<Mesh> meshes;

  //testing reading functions.
  if(!CGAL::IO::read_3MF(filename, meshes))
    return 1;
  for(std::size_t i = 0; i< meshes.size(); ++i)
  {
    Mesh mesh = meshes[i];
    std::cout<<"mesh "<<i<<" is valid: "<<mesh.is_valid()<<std::endl;
    std::string outputName("output");
    outputName.append(std::to_string(i));
    outputName.append(".off");
    std::ofstream ofs(outputName);
    ofs << mesh;
    ofs.close();
  }

  // testing writing functions
  Mesh sphere, tube;
  CGAL::make_icosahedron<Mesh, Point_3>(sphere);
  CGAL::make_regular_prism(10, tube, Point_3(0,-10,0), 10);

  all_points.clear();
  all_polygons.clear();
  all_colors.clear();
  names.clear();

  PointRange points;
  PolygonRange triangles;
  ColorRange colors;

  typedef boost::property_map<Mesh, boost::vertex_point_t>::type VPMap;
  VPMap vpm = get(boost::vertex_point, sphere);

  std::unordered_map<boost::graph_traits<Mesh>::vertex_descriptor, std::size_t> vertex_id_map;
  std::size_t i = 0;
  for(auto v : sphere.vertices())
  {
    points.push_back(get(vpm, v));
    vertex_id_map[v] = i++;
  }
  all_points.push_back(points);

  for(auto f : sphere.faces())
  {
    Polygon triangle;
    for(auto vert : CGAL::vertices_around_face(halfedge(f, sphere), sphere))
    {
      triangle.push_back(vertex_id_map[vert]);
    }
    triangles.push_back(triangle);
    colors.push_back(CGAL::IO::Color(255,0,0,255));
  }

  all_polygons.push_back(triangles);
  all_colors.push_back(colors);
  points.clear();
  triangles.clear();
  colors.clear();
  vertex_id_map.clear();
  i = 0;

  vpm = get(boost::vertex_point, tube);
  for(auto v : tube.vertices())
  {
    points.push_back(get(vpm, v));
    vertex_id_map[v] = i++;
  }
  all_points.push_back(points);

  for(auto f : tube.faces())
  {
    Polygon triangle;
    for(auto vert : CGAL::vertices_around_face(halfedge(f, tube), tube))
    {
      triangle.push_back(vertex_id_map[vert]);
    }
    triangles.push_back(triangle);
    colors.push_back(CGAL::IO::Color(0,0,255,255));

  }
  all_polygons.push_back(triangles);
  all_colors.push_back(colors);
  names.push_back(std::string("sphere"));
  names.push_back(std::string("tube"));

  meshes.resize(2);
  meshes[0] = sphere;
  meshes[1] = tube;

  CGAL::IO::write_3MF("meshes.3mf", meshes, names);

  std::cout << "OK." << std::endl;
#endif //CGAL_LINKED_WITH_3MF

  return 0;
}

