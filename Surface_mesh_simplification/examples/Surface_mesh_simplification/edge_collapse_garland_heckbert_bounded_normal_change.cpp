#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

#include <CGAL/IO/STL_reader.h>
#include <CGAL/IO/PLY_reader.h>

#include <CGAL/Surface_mesh_simplification/GarlandHeckbert_edge_collapse_visitor_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_cost_stop_predicate.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>

typedef CGAL::Simple_cartesian<double>                          Kernel;
typedef Kernel::Point_3                                         Point_3;
typedef CGAL::Surface_mesh<Point_3>                             Surface_mesh;

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

  if(!CGAL::is_triangle_mesh(surface_mesh))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input mesh: " << num_vertices(surface_mesh) << " nv "
                              << num_edges(surface_mesh) << " ne "
                              << num_faces(surface_mesh) << " nf" << std::endl;

  // In this example, the simplification stops when the number of undirected edges
  // drops below 10% of the initial count
  const double stop_ratio = (argc > 2) ? std::stod(argv[2]) : 0.1;

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  // Garland&Heckbert simplification requires a state to be shared between cost,
  // placement, and visitor function objects. This creates the state variable that is
  // required.
  SMS::GarlandHeckbert_edge_collapse_visitor_base<Surface_mesh>::garland_heckbert_state_type state;

  // stop
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(stop_ratio);

  // cost
  SMS::GarlandHeckbert_cost<Surface_mesh> cost(state);

  // placement
  SMS::GarlandHeckbert_placement<Surface_mesh> gh_placement(state);
  SMS::Bounded_normal_change_placement<SMS::GarlandHeckbert_placement<Surface_mesh> > placement(gh_placement);

  // visitor
  // note that the use of Garland&Heckbert visitor is required for it to run.
  SMS::GarlandHeckbert_edge_collapse_visitor_base<Surface_mesh> vis(state);

  int r = SMS::edge_collapse(surface_mesh, stop,
                             CGAL::parameters::get_cost(cost)
                                              .get_placement(placement)
                                              .visitor(vis));

  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

  std::cout << "Time elapsed: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
            << "ms" << std::endl;

  std::cout << "\nFinished...\n" << r << " edges removed.\n" << surface_mesh.number_of_edges() << " final edges.\n";

  std::ofstream os(argc > 3 ? argv[3] : "out.off");
  os.precision(17);
  os << surface_mesh;

  return EXIT_SUCCESS;
}
