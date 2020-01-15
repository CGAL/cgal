#define CGAL_PMP_SNAP_DEBUG
#define CGAL_PMP_SNAP_DEBUG_PP
#define CGAL_PMP_SNAP_DEBUG_OUTPUT

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/IO/STL_reader.h>
#include <CGAL/IO/PLY_reader.h>

#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <boost/property_map/function_property_map.hpp>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                          Kernel;
typedef Kernel::Point_3                                         Point_3;
typedef CGAL::Surface_mesh<Point_3>                             Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor    vertex_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

template <typename K, typename Mesh>
void read_mesh(const char* filename, Mesh& sm)
{
  typedef typename K::Point_3                                   Point;

  std::ifstream in(filename, std::ios::binary);
  if(!in.good())
  {
    std::cerr << "Error: can't read file: " << filename << std::endl;
    std::exit(1);
  }

  std::string fn(filename);
  if(fn.substr(fn.find_last_of(".") + 1) == "stl")
  {
    std::vector<Point> points;
    std::vector<std::vector<int> > faces;
    if(!CGAL::read_STL(in, points, faces))
    {
      std::cerr << "Error: cannot open STL mesh\n";
      return;
    }

    std::cout << "Cleaning polygon soup..." << std::endl;
    PMP::repair_polygon_soup(points, faces);

    if(!CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces))
      std::cerr << "W: File does not describe a polygon mesh" << std::endl;

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces, sm);
  }
  else if(fn.substr(fn.find_last_of(".") + 1) == "off")
  {
    if(!in || !(in >> sm))
    {
      std::cerr << "Error: cannot OFF open mesh\n";
      return;
    }
  }
  else if(fn.substr(fn.find_last_of(".") + 1) == "ply")
  {
    std::vector<Point> points;
    std::vector<std::vector<std::size_t> > polygons;
    std::vector<CGAL::Color> fcolors;
    std::vector<CGAL::Color> vcolors;

    if(!(CGAL::read_PLY(in, points, polygons, fcolors, vcolors)))
    {
      std::cerr << "Error: cannot open PLY mesh\n";
      return;
    }

    std::cout << "Cleaning polygon soup..." << std::endl;
    PMP::repair_polygon_soup(points, polygons);

    if(!CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons))
      std::cerr << "W: File does not describe a polygon mesh" << std::endl;

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, sm);
  }
  else
  {
    std::cerr << "Unknown file type" << std::endl;
    return;
  }

  std::cout << "Input mesh: " << num_vertices(sm) << " nv "
                              << num_edges(sm) << " ne "
                              << num_faces(sm) << " nf" << std::endl;
}

int main(int /*argc*/, char** argv)
{
  Surface_mesh organic_mesh, fixed_mesh;

  const char* organic_filename = argv[1];
  read_mesh<Kernel>(organic_filename, organic_mesh);

  const char* fixed_filename = argv[2];
  read_mesh<Kernel>(fixed_filename, fixed_mesh);

  std::cout << "Organic mesh: " << num_vertices(organic_mesh) << std::endl;
  std::cout << "Fixed mesh: " << num_vertices(fixed_mesh) << std::endl;

  Surface_mesh::Property_map<vertex_descriptor, double> organic_tolerance_map;
  organic_tolerance_map = organic_mesh.add_property_map<vertex_descriptor, double>("v:t").first;
  for(vertex_descriptor v : vertices(organic_mesh))
    put(organic_tolerance_map, v, 0.3);

  Surface_mesh::Property_map<vertex_descriptor, double> fixed_tolerance_map;
  fixed_tolerance_map = fixed_mesh.add_property_map<vertex_descriptor, double>("v:t").first;
  for(vertex_descriptor v : vertices(fixed_mesh))
    put(fixed_tolerance_map, v, 0.3);

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  std::size_t nb_snapped = PMP::experimental::snap_borders(organic_mesh, organic_tolerance_map,
                                                           fixed_mesh, fixed_tolerance_map,
                                                           CGAL::parameters::do_simplify_border(true),
                                                           CGAL::parameters::do_lock_mesh(true));

  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

  std::cout << "Time elapsed: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
            << "ms" << std::endl;

  std::cout << "#Snapped: " << nb_snapped << std::endl;

  std::ofstream("results/snapped_organic.off") << std::setprecision(17) << organic_mesh;
  std::ofstream("results/snapped_fixed.off") << std::setprecision(17) << fixed_mesh;

  return EXIT_SUCCESS;
}
