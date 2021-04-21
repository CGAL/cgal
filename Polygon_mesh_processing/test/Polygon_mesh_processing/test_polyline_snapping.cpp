#include <cstdlib>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/polyline_snapping.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>

#include <boost/graph/adjacency_list.hpp>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;

struct Vertex_property { Point_3 point; };

using Graph = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, Vertex_property>;
using vertex_descriptor = boost::graph_traits<Graph>::vertex_descriptor;
using edge_descriptor = boost::graph_traits<Graph>::edge_descriptor;

class Visitor
{
  const Graph& graph;
  std::ofstream& ofile;
  std::vector<vertex_descriptor> vertices;

public:

  Visitor (const Graph& graph, std::ofstream& ofile) : graph (graph), ofile(ofile) { }

  void start_new_polyline()
  {
    vertices.clear();
  }

  void add_node (vertex_descriptor v)
  {
    vertices.push_back(v);
  }

  void end_polyline()
  {
    ofile << vertices.size();
    for (vertex_descriptor vd : vertices)
      ofile << " " << graph[vd].point;
    ofile << std::endl;
  }
};

int main (int argc, char** argv)
{
  std::vector<std::vector<Point_3> > polylines;

  if (argc == 1)
  {
    std::cerr << "Usage: " << argv[0] << " input1.polylines.txt input2.polylines.txt (...)" << std::endl;
    return EXIT_SUCCESS;
  }

  Graph graph;

  std::size_t nb_polylines = 0;
  std::map<Point_3, vertex_descriptor> map_p2v;

  for (int i = 1; i < argc; ++ i)
  {
    std::ifstream stream(argv[std::size_t(i)]);
    if (!stream)
      std::cerr << "Error: can't read " << argv[std::size_t(i)] << std::endl;

    std::string line;
    std::istringstream iss;
    while(getline(stream,line))
    {
      iss.clear();
      iss.str(line);
      int size = 0;

      if (iss >> size && (size > 0))
      {
        vertex_descriptor previous;
        for (int n = 0; n < size; ++ n)
        {
          double x, y, z;
          if (!(iss >> CGAL::iformat(x) >> CGAL::iformat(y) >> CGAL::iformat(z)))
          {
            std::cerr << "Reading error!" << std::endl;
            return EXIT_FAILURE;
          }

          Point_3 p (x, y, z);
          auto iter = map_p2v.insert (std::make_pair (p, vertex_descriptor()));
          if (iter.second)
          {
            iter.first->second = add_vertex (graph);
            graph[iter.first->second].point = Point_3 (x, y, z);
          }
          vertex_descriptor current = iter.first->second;
          if (n != 0)
            add_edge (previous, current, graph);
          previous = current;

          ++ nb_polylines;
        }
      }
    }
  }

  std::cerr << nb_polylines << " polyline(s) read" << std::endl
            << "  -> stored in graph with " << num_vertices(graph)
            << " vertices and " << num_edges(graph) << " edges" << std::endl;

  std::cerr << "Snapping..." << std::endl;

  CGAL::Polygon_mesh_processing::polyline_snapping (graph, get(&Vertex_property::point, graph), 0.1);

  std::cerr << "  -> resulting graph has " << num_vertices(graph)
            << " vertices and " << num_edges(graph) << " edges" << std::endl;

  std::ofstream ofile ("snapped.polylines.txt");
  Visitor visitor(graph, ofile);
  CGAL::split_graph_into_polylines (graph, visitor);

  return EXIT_SUCCESS;
}
