#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>

#include <boost/graph/graph_test.hpp>
#include <boost/graph/random.hpp>
#include <boost/assign.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Triangulation_2<Kernel> Triangulation;
typedef typename Triangulation::Point Point;
typedef typename Triangulation::Vertex Vertex;
typedef typename Triangulation::Vertex Vertex_handle;

int test_main(int, char*[])
{
  using namespace boost::assign;

  Triangulation t;
  {
    std::vector<Point> v;
    // taken from the triangulation test cases, should be randomly generated
    v += Point(5,6,1), Point(1,9,1), Point(6,14,1), Point(4,12,1), Point(3,29,1),  Point(6,7,1),
      Point(6,39,1), Point(8,9,1), Point(10,18,1), Point(75625,155625,10000),
      Point(10,50,2), Point(6,15,2), Point(6,16,2), Point(10,11,1),
      Point(10,40,1), Point(60,-10,1);

    t.insert(v.begin(), v.end());
  }
  
  typedef typename boost::graph_traits<Triangulation>::vertex_descriptor vertex_t;
  std::vector<vertex_t> vv;
  for(typename Triangulation::Finite_vertices_iterator it = t.finite_vertices_begin(); 
      it != t.finite_vertices_end(); ++it) { 
    vv.push_back(it);
  }


  std::vector< std::pair<vertex_t, vertex_t> > e;
  for(Triangulation::All_edges_iterator it = t.all_edges_begin(); 
      it != t.all_edges_end(); ++it) { 
    e.push_back(
      std::make_pair(
        it->first->vertex((it->second + 2) % 3)
        , it->first->vertex((it->second + 1) % 3)));
  }

  boost::graph_test<Triangulation> gt;
  gt.test_bidirectional_graph(vv, e, t);
  return 0;
}


