//needed by functions
#include <iostream>
#include <vector>
#include <string>
#include "Model/COM/NMR_DLLInterfaces.h"
//needed by example
#include <CGAL/boost/graph/helpers.h>
#include <fstream>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO/3mf.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <unordered_map>
#include <CGAL/IO/read_3mf.h>
#include <CGAL/IO/write_3mf.h>
// Use NMR namespace for the interfaces
using namespace NMR;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef std::vector<Point_3> PointRange;
typedef std::vector<std::size_t> Polygon;
typedef std::vector<Polygon> PolygonRange;
typedef std::vector<CGAL::Color> ColorRange;

int main(int argc, char** argv)
{
  const char* file_name=(argc == 2) ? argv[1] : "data/test.3mf";

  std::vector<PointRange> all_points;
  std::vector<PolygonRange> all_polygons;
  std::vector<ColorRange> all_colors;
  std::vector<std::string> names;
  std::vector<Mesh> meshes;
  //testing reading functions.
  int nb_meshes =
      CGAL::read_3mf(file_name, meshes);
  if(nb_meshes <0)
    return 1;
  for(int i = 0; i< nb_meshes; ++i)
  {
    Mesh mesh = meshes[i];
    std::cout<<names[i]<<" is valid: "<<mesh.is_valid()<<std::endl;
    std::string outputName("output");
    outputName.append(std::to_string(i));
    outputName.append(".off");
    std::ofstream ofs(outputName);
    ofs << mesh;
    ofs.close();
  }
  int nb_polylines =
      CGAL::read_polylines_from_3mf(file_name, all_points, all_colors, names);

  if(nb_polylines == 0)
    std::cout<<"No polyline found."<<std::endl;
  else
  {
    std::cout<<nb_polylines<<" polylines found, of ";
    for(int i = 0; i< nb_polylines-1; ++i){
      std::cout<<all_points[i].size()<<", ";
    }
    std::cout<<all_points.back().size()<<" points."<<std::endl;
  }
  all_points.clear();
  all_colors.clear();
  int nb_point_sets =
      CGAL::read_point_clouds_from_3mf(file_name, all_points, all_colors, names);
  if(nb_point_sets == 0)
    std::cout<<"No point cloud found."<<std::endl;
  else
  {
    std::cout<<nb_point_sets<<" point clouds found, of ";
    for(int i = 0; i< nb_point_sets-1; ++i){
      std::cout<<all_points[i].size()<<", ";
    }
    std::cout<<all_points.back().size()<<" points."<<std::endl;
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
  std::unordered_map<boost::graph_traits<Mesh>::vertex_descriptor,
      std::size_t> vertex_id_map;
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
    colors.push_back(CGAL::Color(255,0,0,255));
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
    colors.push_back(CGAL::Color(0,0,255,255));

  }
  all_polygons.push_back(triangles);
  all_colors.push_back(colors);
  names.push_back(std::string("sphere"));
  names.push_back(std::string("tube"));

  meshes.resize(2);
  meshes[0] = sphere;
  meshes[1] = tube;
  CGAL::write_triangle_meshes_to_3mf("meshes.3mf", meshes, names);


  //testing of point clouds

  HRESULT hResult;
  NMR::PLib3MFModel * pModel;
  hResult = NMR::lib3mf_createmodel(&pModel);
  NMR::PLib3MFModelMeshObject* pMeshObject;
  if (hResult != LIB3MF_OK) {
    std::cout << "could not create model: " << std::hex << hResult << std::endl;
    return 1;
  }
  for(std::size_t i=0; i< names.size(); ++i)
  {
    CGAL::write_mesh_to_model(all_points[i], all_polygons[i],
                              all_colors[i], names[i], &pMeshObject, pModel);
  }
  CGAL::Color color(255,0,0);
  CGAL::write_point_cloud_to_model(all_points.front(),
                                   color, names.front(), &pMeshObject, pModel);
  CGAL::export_model_to_file("micro.3mf", pModel);
  //testing of polylines
  CGAL::write_polyline_to_model(all_points.back(),
                                color, names.back(), &pMeshObject, pModel);
  CGAL::export_model_to_file("micro.3mf", pModel);

  std::cout<<"OK."<<std::endl;
  return 0;
}

