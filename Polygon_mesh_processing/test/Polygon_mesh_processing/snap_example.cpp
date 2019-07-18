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
  Surface_mesh surface_mesh;

  const char* filename = argv[1];
  read_mesh<Kernel>(filename, surface_mesh);

  if(!CGAL::is_triangle_mesh(surface_mesh))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  auto tolerance_map =
      boost::make_function_property_map<Surface_mesh::Vertex_index>(
        [](Surface_mesh::Vertex_index){return 0.001;}
      );

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  std::size_t nb_snapped = PMP::experimental::snap_border_vertices_non_conforming(surface_mesh, tolerance_map);
  std::cout << "initial snapped: " << nb_snapped << std::endl;

  std::chrono::steady_clock::time_point snap_time = std::chrono::steady_clock::now();
  std::cout << "Time elapsed (snap): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(snap_time - start_time).count()
            << "ms" << std::endl;

  PMP::stitch_borders(surface_mesh);

  std::chrono::steady_clock::time_point stitch_time = std::chrono::steady_clock::now();
  std::cout << "Time elapsed (stitch): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(stitch_time - snap_time).count()
            << "ms" << std::endl;

  std::ofstream out("results/snapped.off");
  out.precision(17);
  out << surface_mesh;
  out.close();

  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

  std::cout << "#border: " << PMP::number_of_borders(surface_mesh) << std::endl;

  std::cout << "Time elapsed: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
            << "ms" << std::endl;

  return EXIT_SUCCESS;
}
