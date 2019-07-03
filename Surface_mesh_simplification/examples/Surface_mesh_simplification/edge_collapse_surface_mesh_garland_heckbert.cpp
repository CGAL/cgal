#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

#include <CGAL/IO/STL_reader.h>

#include <CGAL/Surface_mesh_simplification/GarlandHeckbert_edge_collapse_visitor_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_cost_stop_predicate.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <chrono>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

namespace SMS = CGAL::Surface_mesh_simplification;

template <typename Kernel, typename Mesh>
void read_mesh(const char* filename,
               Mesh& sm)
{
  typedef typename Kernel::Point_3                                    Point;

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
  else // off reading
  {
    if(!in || !(in >> sm))
    {
      std::cerr << "Error: cannot open mesh\n";
      return;
    }
  }
}

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;

  const char* filename = argv[1];
  read_mesh<Kernel>(filename, surface_mesh);

  if(!CGAL::is_triangle_mesh(surface_mesh))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  double stop_ratio = 0.1;
  if(argc > 2) {
    stop_ratio = std::stod(argv[2]);
  }

  std::chrono::steady_clock::time_point start_time
    = std::chrono::steady_clock::now();

  // In this example, the simplification stops when the number of undirected edges
  // drops below 10% of the initial count
  double threshold = (argc > 2) ? std::atof(argv[2]) : 0.1;
//  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(threshold);
  SMS::GarlandHeckbert_cost_stop_predicate<double> stop(threshold);

  SMS::GarlandHeckbert_edge_collapse_visitor_base<Surface_mesh>::garland_heckbert_map_type map;
  SMS::GarlandHeckbert_edge_collapse_visitor_base<Surface_mesh> vis(map);

  int r = SMS::edge_collapse(surface_mesh, stop,
                             CGAL::parameters::get_cost(SMS::GarlandHeckbert_cost <Surface_mesh>(map))
                                              .get_placement(SMS::GarlandHeckbert_placement<Surface_mesh>(map))
                                              .visitor(vis));

  std::chrono::steady_clock::time_point end_time
    = std::chrono::steady_clock::now();



  std::cout << "Time elapsed: "
   << std::chrono::duration_cast<std::chrono::milliseconds>(
         end_time - start_time
       ).count() << "ms" << std::endl;


  std::cout << "\nFinished...\n" << r << " edges removed.\n"
            << surface_mesh.number_of_edges() << " final edges.\n";

  std::ofstream os(argc > 3 ? argv[3] : "out.off");
  os.precision(17);
  os << surface_mesh;

  return EXIT_SUCCESS;
}
