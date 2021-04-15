#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>

#include <CGAL/boost/graph/copy_face_graph.h>

#include <iostream>
#include <fstream>
#include <iterator>
#include <boost/unordered_map.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef CGAL::Polyhedron_3<Kernel>                       Source;
typedef boost::graph_traits<Source>::vertex_descriptor   sm_vertex_descriptor;
typedef boost::graph_traits<Source>::halfedge_descriptor sm_halfedge_descriptor;
typedef boost::graph_traits<Source>::face_descriptor     sm_face_descriptor;

typedef CGAL::Exact_predicates_exact_constructions_kernel Other_kernel;
typedef Other_kernel::Point_3                             Point;
typedef CGAL::Linear_cell_complex_traits<3, Other_kernel> LCC_traits;
typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
         <2, 3, LCC_traits>::type LCC;

int main(int argc, char* argv[])
{
  Source S;

  std::ifstream in((argc>1)?argv[1]:"cube.off");
  in >> S;

  // Note that the vertex_point property of the Source and Target1
  // come from different kernels.
  typedef LCC Target1;
  Target1 T1;
  {
    CGAL::copy_face_graph(S, T1);
    CGAL::write_off("lcc.off", T1);
  }

  S.clear();
  {
    typedef boost::graph_traits<Target1>::vertex_descriptor   source_vertex_descriptor;
    typedef boost::graph_traits<Target1>::halfedge_descriptor source_halfedge_descriptor;

    typedef boost::graph_traits<Source>::vertex_descriptor   tm_vertex_descriptor;
    typedef boost::graph_traits<Source>::halfedge_descriptor tm_halfedge_descriptor;

    boost::unordered_map<source_vertex_descriptor, tm_vertex_descriptor> v2v;
    boost::unordered_map<source_halfedge_descriptor, tm_halfedge_descriptor> h2h;

    CGAL::copy_face_graph(T1, S, CGAL::parameters::vertex_to_vertex_output_iterator(std::inserter(v2v, v2v.end()))
                                                  .halfedge_to_halfedge_output_iterator(std::inserter(h2h, h2h.end())));
    std::ofstream out("reverse.off");
    out.precision(17);
    out << S;
  }
  return 0;
}
