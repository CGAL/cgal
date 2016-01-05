#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <fstream>

#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                        Kernel;
typedef Kernel::Point_3                                       Point;
typedef CGAL::Polyhedron_3<Kernel>                            Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;

typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron> Skeletonization;
typedef Skeletonization::Skeleton                             Skeleton;

typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
typedef Skeleton::edge_descriptor                             Skeleton_edge;

//only needed for the display of the skeleton as maximal polylines
struct Display_polylines{
  const Skeleton& skeleton;
  std::ofstream& out;
  int polyline_size;
  std::stringstream sstr;

  Display_polylines(const Skeleton& skeleton, std::ofstream& out)
    : skeleton(skeleton), out(out)
  {}

  void start_new_polyline(){
    polyline_size=0;
    sstr.str("");
    sstr.clear();
  }
  void add_node(Skeleton_vertex v){
    ++polyline_size;
    sstr << " " << skeleton[v].point;
  }
  void end_polyline()
  {
    out << polyline_size << sstr.str() << "\n";
  }
};

// This example extracts a medially centered skeleton from a given mesh.
int main(int argc, char* argv[])
{
  std::ifstream input((argc>1)?argv[1]:"data/elephant.off");
  Polyhedron tmesh;
  input >> tmesh;

  Skeleton skeleton;

  CGAL::extract_mean_curvature_flow_skeleton(tmesh, skeleton);

  std::cout << "Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
  std::cout << "Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";

  // Output all the edges of the skeleton.
  std::ofstream output("skel.cgal");
  Display_polylines display(skeleton,output);
  CGAL::split_graph_into_polylines(skeleton, display);
  output.close();

  // Output skeleton points and the corresponding surface points
  output.open("correspondance.cgal");
  BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton))
    BOOST_FOREACH(vertex_descriptor vd, skeleton[v].vertices)
      output << "2 " << skeleton[v].point << " "
                     << get(CGAL::vertex_point, tmesh, vd)  << "\n";

  return 0;
}

