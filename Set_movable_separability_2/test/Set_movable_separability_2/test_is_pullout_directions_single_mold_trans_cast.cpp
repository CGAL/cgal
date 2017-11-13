#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <utility>
#include <cctype>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Set_movable_separability_2/Single_mold_translational_casting/top_edges.h>
#include <CGAL/Set_movable_separability_2/Single_mold_translational_casting/is_pullout_direction.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef Kernel::Direction_2                               Direction_2;
typedef Kernel::Point_2                                   Point_2;

typedef std::pair<Direction_2, Direction_2>               Direction_range;
typedef Polygon_2::Edge_const_iterator                    Edge_iter;
typedef std::pair<Edge_iter, Direction_range>             Top_edge;

namespace SMS = CGAL::Set_movable_separability_2;
namespace casting = SMS::Single_mold_translational_casting;

bool test_one_file(std::ifstream& inp)
{
  Polygon_2 pgn;
  inp >> pgn;
  // std::cout << pgn << std::endl;

  std::vector<Top_edge> top_edges;
  casting::top_edges(pgn, std::back_inserter(top_edges));

  for (auto it = top_edges.begin(); it != top_edges.end(); ++it) {
    auto facet = it->first;
    const auto& d1 = it->second.first;
    const auto& d2 = it->second.second;

    if (! casting::is_pullout_direction(pgn, facet, d1)) return false;
    if (! casting::is_pullout_direction(pgn, facet, d2)) return false;

    auto od1 = CGAL::opposite(d1);
    if (casting::is_pullout_direction(pgn, facet, od1)) return false;
    auto od2 = CGAL::opposite(d2);
    if (casting::is_pullout_direction(pgn, facet, od2)) return false;
  }
  return true;
}

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Missing input file" << std::endl;
    return -1;
  }

  int success = 0;
  for (size_t i = 1; i < static_cast<size_t>(argc); ++i) {
    std::string str(argv[i]);
    if (str.empty()) continue;

    auto itr = str.end();
    --itr;
    while (itr != str.begin()) {
      auto tmp = itr;
      --tmp;
      if (!isspace(*itr)) break;
      str.erase(itr);
      itr = tmp;
    }
    if (str.size() <= 1) continue;
    std::ifstream inp(str.c_str());
    if (!inp.is_open()) {
      std::cerr << "Failed to open " << str << std::endl;
      return -1;
    }
    if (! test_one_file(inp)) {
      std::cout << str << ": ERROR" << std::endl;
      ++success;
    }
    else std::cout << str << ": succeeded" << std::endl;
    inp.close();
  }

  return success;
}
