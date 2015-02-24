#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Mean_curvature_skeleton_functions.h>

#include <fstream>

typedef CGAL::Simple_cartesian<double>                               Kernel;
typedef Kernel::Point_3                                              Point;
typedef CGAL::Polyhedron_3<Kernel>                                   Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor           vertex_descriptor;

typedef CGAL::Mean_curvature_flow_skeletonization<Kernel,Polyhedron> Mean_curvature_skeleton;
typedef Mean_curvature_skeleton::Skeleton                            Skeleton;

typedef boost::graph_traits<Skeleton>::vertex_descriptor             vertex_desc;
typedef boost::graph_traits<Skeleton>::vertex_iterator               vertex_iter;
typedef boost::graph_traits<Skeleton>::edge_iterator                 edge_iter;


// This example extracts a medially centered skeleton from a given mesh.
int main()
{
  Polyhedron mesh;
  std::ifstream input("data/sindorelax.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Cannot open data/sindorelax.off" << std::endl;
    return 1;
  }

  Skeleton skeleton;

  CGAL::extract_mean_curvature_flow_skeleton(mesh, skeleton);

  std::cout << "vertices: " << boost::num_vertices(skeleton) << "\n";
  std::cout << "edges: " << boost::num_edges(skeleton) << "\n";

  // Output all the edges.
  edge_iter ei, ei_end;
  for (boost::tie(ei, ei_end) = boost::edges(skeleton); ei != ei_end; ++ei)
  {
    Point s = skeleton[source(*ei, skeleton)].point;
    Point t = skeleton[target(*ei, skeleton)].point;
    std::cout << s << " " << t << "\n";
  }


  // Output skeletal points and the corresponding surface points
  vertex_iter gvb, gve;
  for (boost::tie(gvb, gve) = boost::vertices(skeleton); gvb != gve; ++gvb)
  {
    vertex_desc gv = *gvb;
    Point skel = skeleton[gv].point;
    std::cout << skel << ": ";

    for (size_t i = 0; i < skeleton[gv].vertices.size(); ++i)
    {
      vertex_descriptor v = skeleton[gv].vertices[i];
      Point surf = v ->point();
      std::cout << surf << " ";
    }
    std::cout << "\n";
  }
  return 0;
}

