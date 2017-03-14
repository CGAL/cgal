#include <list>
#include <fstream>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/top_edges_single_mold_translational_casting_2.h>
#include <CGAL/pullout_directions_single_mold_translational_casting_2.h>
#include <CGAL/is_pullout_direction_single_mold_translational_casting_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef Kernel::Direction_2                               Direction_2;
typedef Kernel::Vector_2                              	  Vector_2;
typedef Kernel::Point_2                              	  Point_2;
typedef std::pair<Direction_2, Direction_2>               Direction_range;
typedef Polygon_2::Edge_const_iterator Edge_iter;
typedef std::pair< Edge_iter, Direction_range>                Top_edge;

namespace SMS = CGAL::Set_movable_separability_2;

// The main program:
int main(int  argc, char* argv[])
{
  Polygon_2 pgn;

  const char* filename = (argc > 1) ? argv[1] : "polygon.dat";
  std::ifstream input_file(filename);
  if (! input_file.is_open()) {
      std::cerr << "Failed to open the " << filename << std::endl;
      return -1;
  }
  input_file >> pgn;
  input_file.close();

  auto poly_orientation = pgn.orientation();
  std::list<Top_edge> top_edges;


  //example for top_edges_single_mold_translational_casting_2
  SMS::top_edges_single_mold_translational_casting_2
  (pgn, std::back_inserter(top_edges));
  if (top_edges.empty())
    std::cout << "The polygon is not castable!" << std::endl;
  else {
      std::cout << "There are " << top_edges.size() << " top edges:" << std::endl;
      for (const auto& top_edge : top_edges) {
	  std::cout
	      << "\tEdge: "<< *top_edge.first<< std::endl
	      << "\tPullout directions from: "<< top_edge.second.first
	      << " to " << top_edge.second.second
	      << std::endl<< std::endl;
      }
  }
  std::cout << "-----------------------------------"<< std::endl;

  //example for pullout_directions_single_mold_translational_casting_2
  int index =0;
  for (auto e_it = pgn.edges_begin(); e_it != pgn.edges_end(); ++e_it, ++index)

    {

      std::pair<bool, std::pair< Kernel::Direction_2,
      Kernel::Direction_2> > res = SMS::pullout_directions_single_mold_translational_casting_2(pgn,e_it);
      if (res.first)
	{
	  std::cout << "The polygon is castable using edge "<<index<<" in range " << res.second.first
	      << " to " << res.second.second
	      << std::endl<< std::endl;
	}
      else {
	  std::cout << "The polygon is not castable using edge "<<index<<std::endl;
      }
    }
  std::cout << "-----------------------------------"<< std::endl;
  //example for is_pullout_directions_single_mold_translational_casting_2 that accepts the edge
   index =0;
  for (auto e_it = pgn.edges_begin(); e_it != pgn.edges_end(); ++e_it, ++index)
    {
      auto segment_outer_circle =
	  SMS::internal::get_segment_outer_circle<Kernel>(*e_it, poly_orientation);
      Direction_2  d = segment_outer_circle.first;
      d= d.perpendicular(CGAL::CLOCKWISE);
      bool res = SMS::is_pullout_direction_single_mold_translational_casting_2(pgn,e_it,d);
      std::cout << "The polygon is "<<(res?"":"not ") <<"castable using edge "<<index<<" in vartical translation ("<<d<<")"<< std::endl;

    }
  std::cout << "-----------------------------------"<< std::endl;

  //example for is_pullout_directions_single_mold_translational_casting_2 that do not accepts the edge
  {
    Vector_2 v (Point_2(0,0),Point_2(1,0));
    Direction_2  d(v);
    CGAL::Polygon_2<Kernel>::Edge_const_iterator
    res = SMS::is_pullout_direction_single_mold_translational_casting_2(pgn,d);
    if (res!= pgn.edges_end())
      {
	std::cout << "The polygon is castable in direction d ("<<d<<") using edge "<<*res<<std::endl;

      }
    else {
	std::cout << "The polygon is not castable in direction d ("<<d<<")"<<std::endl;
    }
  }

  return 0;
}
