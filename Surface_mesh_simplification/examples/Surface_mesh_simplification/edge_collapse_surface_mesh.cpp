#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

#include <chrono>
#include <fstream>
#include <iostream>

typedef CGAL::Simple_cartesian<double>               Kernel;
typedef Kernel::Point_3                              Point_3;
typedef CGAL::Surface_mesh<Point_3>                  Surface_mesh;

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

  if(!in || !(in >> sm))
  {
    std::cerr << "Error: cannot read OFF mesh\n";
    return;
  }
}

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;

  const char* filename = argv[1];
  read_mesh<Kernel>(filename, surface_mesh);

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  // In this example, the simplification stops when the number of undirected edges
  // drops below 10% of the initial count
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(0.05);

  int r = SMS::edge_collapse(surface_mesh, stop);

  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

  std::cout << "\nFinished...\n" << r << " edges removed.\n"
            << surface_mesh.number_of_edges() << " final edges.\n";

  std::cout << "Time elapsed: "
   << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << "ms" << std::endl;

  std::ofstream os(argc > 2 ? argv[2] : "out.off");
  os.precision(17);
  os << surface_mesh;

  return EXIT_SUCCESS;
}
