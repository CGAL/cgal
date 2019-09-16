#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Mean_curvature_flow_skeletonization.h>

#include <fstream>

typedef CGAL::Simple_cartesian<double>                        Kernel;
typedef Kernel::Point_3                                       Point;
typedef CGAL::Linear_cell_complex_traits<3, Kernel>           MyTraits;
typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
          <2, 3, MyTraits>::type LCC;
typedef boost::graph_traits<LCC>::vertex_descriptor    vertex_descriptor;

typedef CGAL::Mean_curvature_flow_skeletonization<LCC> Skeletonization;
typedef Skeletonization::Skeleton                             Skeleton;

typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
typedef Skeleton::edge_descriptor                             Skeleton_edge;


int main()
{
  LCC lcc;
  CGAL::read_off("data/elephant.off", lcc);

  Skeleton skeleton;
  Skeletonization mcs(lcc);

  // 1. Contract the mesh by mean curvature flow.
  mcs.contract_geometry();

  // 2. Collapse short edges and split bad triangles.
  mcs.collapse_edges();
  mcs.split_faces();

  // 3. Fix degenerate vertices.
  mcs.detect_degeneracies();

  // Perform the above three steps in one iteration.
  mcs.contract();

  // Iteratively apply step 1 to 3 until convergence.
  mcs.contract_until_convergence();

  // Convert the contracted mesh into a curve skeleton and
  // get the correspondent surface points
  mcs.convert_to_skeleton(skeleton);

  std::cout << "Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
  std::cout << "Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";

  // Output all the edges of the skeleton.
  std::ofstream output("skel-lcc.cgal");
  for(Skeleton_edge e : CGAL::make_range(edges(skeleton)))
  {
    const Point& s = skeleton[source(e, skeleton)].point;
    const Point& t = skeleton[target(e, skeleton)].point;
    output << "2 "<< s << " " << t << "\n";
  }
  output.close();

  // Output skeleton points and the corresponding surface points
  output.open("correspondance-lcc.cgal");
  for(Skeleton_vertex v : CGAL::make_range(vertices(skeleton)))
    for(vertex_descriptor vd : skeleton[v].vertices)
      output << "2 " << skeleton[v].point << "  " << get(CGAL::vertex_point, lcc, vd)  << "\n";

  return 0;
}

