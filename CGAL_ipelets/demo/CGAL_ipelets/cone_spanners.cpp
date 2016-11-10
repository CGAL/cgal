#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/CGAL_Ipelet_base.h>
#include <CGAL/Construct_theta_graph_2.h>
#include <CGAL/Construct_yao_graph_2.h>
#include <CGAL/Compute_cone_boundaries_2.h>
#include <CGAL/property_map.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace CGAL_cone_spanners {

typedef CGAL::Exact_predicates_inexact_constructions_kernel        Kernel;
typedef Kernel::Point_2                                            Point_2;
typedef Kernel::Direction_2                                        Direction_2;

/* Note: due to a bug in the boost library, using a directed graph
 * will cause a compilation error with g++ and clang++ when using c++11 standard.
 * See http://lists.boost.org/Archives/boost/2016/05/229458.php.
 */
typedef boost::adjacency_list<boost::listS,
                              boost::vecS,
                              #ifdef CGAL_CXX11
                              boost::undirectedS,
                              #else
                              boost::directedS,
                              #endif
                              Point_2
                             > Graph;

const std::string labels[] = {  "Theta-k-graph", "Yao-k-graph", "Half-theta-k-graph", "Half-Yao-k-graph", "k cones", "Help" };
const std::string hmsg[] = {
  "Compute a theta-graph with k cones.",
  "Compute a Yao-graph with k cones.",
  "Compute an half-theta-graph with k cones.",
  "Compute an half-Yao-graph with k cones.",
  "Draw k cones around the points.",
};

class Cone_spanners_ipelet
  : public CGAL::Ipelet_base<Kernel,6> {
public:
  Cone_spanners_ipelet()
    :CGAL::Ipelet_base<Kernel,6>("Cone Spanners",labels,hmsg){}
  void protected_run(int);
private:
};

void Cone_spanners_ipelet::protected_run(int fn)
{
  std::vector<Point_2> lst;
  int number_of_cones;
  switch (fn){
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    {
      std::vector<Point_2> points_read;
      read_active_objects(
        CGAL::dispatch_or_drop_output<Point_2>(std::back_inserter(points_read))
      );

      if (points_read.empty()) {
        print_error_message("No mark selected");
        return;
      }
      for(std::vector<Point_2>::iterator it = points_read.begin(); it != points_read.end(); it++) {
        if(std::find(points_read.begin(), it, *it) == it) {
          lst.push_back(*it);
        }
      }

      int ret_val;
      boost::tie(ret_val,number_of_cones)=request_value_from_user<int>("Enter the number of cones");
      if (ret_val < 0) {
        print_error_message("Incorrect value");
        return;
      }
      if(number_of_cones < 2) {
        print_error_message("The number of cones must be larger than 1!");
        return;
      }
      break;
    }
    case 5:
      show_help();
      return;
  }

  if(fn >= 0 && fn <= 3) {
    bool half = (fn == 2 || fn == 3);

    Graph g;
    switch (fn){
      case 0:
      case 2:
      {
        CGAL::Construct_theta_graph_2<Kernel, Graph> theta(number_of_cones, Direction_2(1,0), half);
        theta(lst.begin(), lst.end(), g);
        break;
      }
      case 1:
      case 3:
      {
        CGAL::Construct_yao_graph_2<Kernel, Graph> yao(number_of_cones, Direction_2(1,0), half);
        yao(lst.begin(), lst.end(), g);
        break;
      }
    }
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
      boost::graph_traits<Graph>::edge_descriptor e = *ei;
      boost::graph_traits<Graph>::vertex_descriptor  u = source(e, g);
      boost::graph_traits<Graph>::vertex_descriptor  v = target(e, g);
      draw_in_ipe(Segment_2(g[u], g[v]));
    }
    group_selected_objects_();
  }
  else if(fn == 4) {
    CGAL::Compute_cone_boundaries_2<Kernel> cones;
    std::vector<Direction_2> directions(number_of_cones);
    cones(number_of_cones, Direction_2(1,0), directions.begin());
    for(std::vector<Point_2>::iterator it = lst.begin(); it != lst.end(); it++) {
      for(std::vector<Direction_2>::iterator dir = directions.begin(); dir != directions.end(); dir++) {
        draw_in_ipe(Segment_2(*it,*it + 100*dir->to_vector()));
      }
      group_selected_objects_();
      get_IpePage()->deselectAll();
    }
  }
}


}

CGAL_IPELET(CGAL_cone_spanners::Cone_spanners_ipelet)
