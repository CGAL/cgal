#include <CGAL/Polygon_mesh_processing/refine_mesh_at_isolevel.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>

#include <fstream>
#include <iostream>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point_3;
typedef CGAL::Surface_mesh<Point_3>                          Triangle_mesh;

typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Triangle_mesh>::edge_descriptor edge_descriptor;

typedef Triangle_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;

int main()
{
  const std::string filename = CGAL::data_file_path("meshes/elephant.off");

  Triangle_mesh tm;
  if(!CGAL::IO::read_polygon_mesh(filename, tm) ||
     CGAL::is_empty(tm) || !CGAL::is_triangle_mesh(tm))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }
  //property map for the distance values to the source set
  Vertex_distance_map vertex_distance = tm.add_property_map<vertex_descriptor, double>("v:distance", 5).first;
  std::vector<int> zero_vids={
    /*a closed polyline*/ 2144,145,2690,1752,339,215,1395,338,77,2145,2052,2054,343,1936,22,1751,214,1499,142,358,2694,1750,301,65,59,2650,2060,205,2651,2061,2490,1939,898,13,298,
    /*two adjacent triangles*/ 532, 185, 534, 2735,
    /*another patch with missing crossed edges*/134,73,1883,2533,72,532,185,131,534
  };

  std::vector<int> minus_vids = {132, 364};

  for (int i : zero_vids)
    put(vertex_distance, Triangle_mesh::Vertex_index(i), 0);
  for (int i : minus_vids)
    put(vertex_distance, Triangle_mesh::Vertex_index(i), -5);

  // property map to flag new cut edge added in the mesh
  auto ecm = tm.add_property_map<edge_descriptor, bool>("e:is_constrained", 0).first;

  CGAL::Polygon_mesh_processing::refine_mesh_at_isolevel(tm, vertex_distance, 0, CGAL::parameters::edge_is_constrained_map(ecm));

  std::ofstream debug("edges.polylines.txt");
  for (Triangle_mesh::Edge_index e : edges(tm))
    if (get(ecm, e))
      debug << "2 " << tm.point(source(e, tm)) << " " <<  tm.point(target(e, tm)) << "\n";
  debug.close();

  // split the mesh in connected components bounded by the isocurves
  std::vector<Triangle_mesh> edges_split;
  CGAL::Polygon_mesh_processing::split_connected_components(tm, edges_split, CGAL::parameters::edge_is_constrained_map(ecm));

  // export each submesh in a file
  for(std::size_t i=0; i<edges_split.size(); ++i)
    std::ofstream("out_"+std::to_string(i)+".off") << std::setprecision(17) << edges_split[i];

  assert(edges_split.size() == 5);


  return 0;
}
