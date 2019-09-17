#include <iostream>
#include <fstream>
#include <chrono>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

#include <CGAL/IO/STL_reader.h>
#include <CGAL/IO/PLY_reader.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

namespace SMS = CGAL::Surface_mesh_simplification;

template <typename K, typename Mesh>
void read_mesh(const char* filename,
               Mesh& sm)
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
    std::vector<std::array<int, 3> > faces;
    CGAL::read_STL(in, points, faces);

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

    if(!CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons))
      std::cerr << "W: File does not describe a polygon mesh" << std::endl;

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, sm);
  }
  else
  {
    std::cerr << "Unknown file type" << std::endl;
    return;
  }
}

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;

  const char* filename = argv[1];
  read_mesh<Kernel>(filename, surface_mesh);

  std::chrono::steady_clock::time_point start_time
    = std::chrono::steady_clock::now();


  // In this example, the simplification stops when the number of undirected edges
  // drops below 10% of the initial count
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(0.05);

  int r = SMS::edge_collapse(surface_mesh, stop);

  std::chrono::steady_clock::time_point end_time
    = std::chrono::steady_clock::now();



  std::cout << "Time elapsed: "
   << std::chrono::duration_cast<std::chrono::milliseconds>(
         end_time - start_time
       ).count() << "ms" << std::endl;


  std::cout << "\nFinished...\n" << r << " edges removed.\n"
            << surface_mesh.number_of_edges() << " final edges.\n";

  std::ofstream os(argc > 2 ? argv[2] : "out.off");
  os.precision(17);
  os << surface_mesh;

  return EXIT_SUCCESS;
}
