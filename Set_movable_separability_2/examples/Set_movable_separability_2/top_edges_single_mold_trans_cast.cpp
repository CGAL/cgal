#include <list>
#include <fstream>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Set_movable_separability_2/Single_mold_translational_casting/top_edges.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef Kernel::Direction_2                               Direction_2;
typedef Kernel::Vector_2                                        Vector_2;
typedef Kernel::Point_2                                        Point_2;

// A direction range is a closed range of directions on the unit circle.
typedef std::pair<Direction_2, Direction_2>               Direction_range;
typedef Polygon_2::Edge_const_iterator                    Edge_iter;

// A top edge is identified by the index to an edge of a polygon and the
// corresponding range of pullout directions.
typedef std::pair<Edge_iter, Direction_range>             Top_edge;

namespace SMS = CGAL::Set_movable_separability_2;
namespace casting = SMS::Single_mold_translational_casting;

// The main program:
int main(int  argc, char* argv[])
{
  Polygon_2 polygon;

  const char* filename = (argc > 1) ? argv[1] : "polygon.dat";
  std::ifstream input_file(filename);
  if (! input_file.is_open()) {
    std::cerr << "Failed to open the " << filename << std::endl;
    return -1;
  }
  input_file >> polygon;
  input_file.close();

  std::list<Top_edge> top_edges;

  // Example for top_edges_single_mold_translational_casting_2
  casting::top_edges(polygon, std::back_inserter(top_edges));
  if (top_edges.empty())
    std::cout << "The polygon is not castable!" << std::endl;
  else {
    std::cout << "There are " << top_edges.size() << " top edges:" << std::endl;
    for (const auto& top_edge : top_edges) {
      std::cout
        << "\tEdge: " << *top_edge.first<< std::endl
        << "\tPullout directions from: " << top_edge.second.first
        << " to " << top_edge.second.second
        << std::endl << std::endl;
    }
  }

  return 0;
}
