#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/cluster_point_set.h>
#include <CGAL/compute_average_spacing.h>

#include <CGAL/Random.h>
#include <CGAL/Real_timer.h>

#include <fstream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Point_set = CGAL::Point_set_3<Point_3>;

int main (int argc, char** argv)
{
  // Read input file
  std::ifstream ifile ((argc > 1) ? argv[1] : "data/hippo1.ply", std::ios_base::binary);
  Point_set points;
  ifile >> points;

  // Add a cluster map
  Point_set::Property_map<int> cluster_map = points.add_property_map<int> ("cluster", -1).first;

  // Compute average spacing
  double spacing = CGAL::compute_average_spacing<CGAL::Parallel_if_available_tag> (points, 12);
  std::cerr << "Spacing = " << spacing << std::endl;

  // Adjacencies stored in vector
  std::vector<std::pair<std::size_t, std::size_t> > adjacencies;

  // Compute clusters
  CGAL::Real_timer t;
  t.start();
  std::size_t nb_clusters
    = CGAL::cluster_point_set (points, cluster_map,
                               points.parameters().neighbor_radius(spacing).
                               adjacencies (std::back_inserter (adjacencies)));
  t.stop();
  std::cerr << "Found " << nb_clusters << " clusters with " << adjacencies.size()
            << " adjacencies in " << t.time() << " seconds" << std::endl;

  // Output a colored PLY file
  Point_set::Property_map<unsigned char> red = points.add_property_map<unsigned char>("red", 0).first;
  Point_set::Property_map<unsigned char> green = points.add_property_map<unsigned char>("green", 0).first;
  Point_set::Property_map<unsigned char> blue = points.add_property_map<unsigned char>("blue", 0).first;
  for (Point_set::Index idx : points)
  {
    // One color per cluster
    CGAL::Random rand (cluster_map[idx]);
    red[idx] = rand.get_int(64, 192);
    green[idx] = rand.get_int(64, 192);
    blue[idx] = rand.get_int(64, 192);
  }

  std::ofstream ofile ("out.ply", std::ios_base::binary);
  CGAL::set_binary_mode (ofile);
  ofile << points;

  return EXIT_SUCCESS;
}
