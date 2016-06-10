/*! \file find_single_mold_translational_casting_2.cpp
 * .
 */

#include <list>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/find_single_mold_translational_casting_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef Kernel::Direction_2                               Direction_2;
typedef Kernel::Point_2                                   Point_2;

typedef std::pair<Direction_2, Direction_2>               Direction_range;
typedef std::pair<size_t, Direction_range>                Top_edge;

// The main program:
int main()
{
  Polygon_2 pgn;
//  pgn.push_back(Point_2(0, 0));
//  pgn.push_back(Point_2(1, 0));
//  pgn.push_back(Point_2(1, 1));
//  pgn.push_back(Point_2(0, 1));


//  pgn.push_back(Point_2(0, 7));
//    pgn.push_back(Point_2(1, 0));
//    pgn.push_back(Point_2(1, 1));
//    pgn.push_back(Point_2(0, 1));
    pgn.push_back(Point_2(-1, 7));
        pgn.push_back(Point_2(0, 6));

    pgn.push_back(Point_2(3, 8));
    pgn.push_back(Point_2(1, 1));
    pgn.push_back(Point_2(1, 0));
    pgn.push_back(Point_2(0, 0));
  std::list<Top_edge> top_edges;
  find_single_mold_translational_casting_2(pgn, std::back_inserter(top_edges));

  if (top_edges.empty())
    std::cout << "The polygon is not castable!" << std::endl;
  else {
    std::cout << "There are " << top_edges.size() << " top edges:" << std::endl;
    for (const auto& top_edge : top_edges) {
      std::cout <<pgn.edge(top_edge.first) <<" X "<< top_edge.first << ", ("
                << top_edge.second.first << "," << top_edge.second.second
                << ")" << std::endl;
    }
  }

  return 0;
}
