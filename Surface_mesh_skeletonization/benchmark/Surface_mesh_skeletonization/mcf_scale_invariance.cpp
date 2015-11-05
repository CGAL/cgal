#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Mean_curvature_flow_skeletonization.h>

#include <fstream>
#include <sstream>

typedef CGAL::Simple_cartesian<double>                        Kernel;
typedef Kernel::Point_3                                       Point;
typedef CGAL::Polyhedron_3<Kernel>                            Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;

typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron> Skeletonization;
typedef Skeletonization::Skeleton                             Skeleton;

typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
typedef Skeleton::edge_descriptor                             Skeleton_edge;


int main(int argc, char* argv[])
{
  std::ifstream input((argc>1)?argv[1]:"data/elephant.off");
  Polyhedron tmesh;
  input >> tmesh;

  //scale off
  double factor = atof(argv[2]);

  CGAL::Bbox_3 bbox1=CGAL::bbox_3(tmesh.points_begin(), tmesh.points_end());
  Kernel::Vector_3 to_origin(-bbox1.xmin(), -bbox1.ymin(), -bbox1.zmin());
  
  for(Polyhedron::Point_iterator pit=tmesh.points_begin(),
                                 pit_end=tmesh.points_end();pit!=pit_end;++pit)
  {
    Point new_point = *pit+to_origin;
    double x=new_point.x();
    double y=new_point.y();
    double z=new_point.z();
    
    x = x * factor;
    y = y * factor;
    z = z * factor;
    
    *pit=Point(x, y, z) - to_origin ;
  }

  CGAL::Bbox_3 bbox=CGAL::bbox_3(tmesh.points_begin(), tmesh.points_end());
  double l=std::sqrt(
    CGAL::square(bbox.xmin()-bbox.xmax()) +
    CGAL::square(bbox.ymin()-bbox.ymax()) +
    CGAL::square(bbox.zmin()-bbox.zmax())
  );
  std::cout << "l is " << l << "\n";
  Skeleton skeleton;
  Skeletonization mcs(tmesh);

  mcs.set_is_medially_centered(false); //just to be sure
  //~ mcs.set_quality_speed_tradeoff(0.1 * std::pow(l,-0.01));
  mcs.set_quality_speed_tradeoff(20);
  //~ mcs.set_medially_centered_speed_tradeoff(0.2 * std::pow(l,-0.75));

  // Iteratively apply step 1 to 3 until convergence.
  mcs.contract_until_convergence();

  // Convert the contracted mesh into a curve skeleton and
  // get the correspondent surface points
  mcs.convert_to_skeleton(skeleton);

  std::cout << "Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
  std::cout << "Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";


//scale skelton  
  BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton))
  {
    Point new_point = skeleton[v].point+to_origin;
    double x=new_point.x();
    double y=new_point.y();
    double z=new_point.z();
    
    x = x / factor;
    y = y / factor;
    z = z / factor;
    
    skeleton[v].point=Point(x, y, z)-to_origin;
  }

  // Output all the edges of the skeleton.
  std::stringstream ss;
  ss << "skel." << factor << ".cgal"; 
  std::ofstream output(ss.str().c_str());
  BOOST_FOREACH(Skeleton_edge e, edges(skeleton))
  {
    const Point& s = skeleton[source(e, skeleton)].point;
    const Point& t = skeleton[target(e, skeleton)].point;
    output << "2 "<< s << " " << t << "\n";
  }
  output.close();

  // Output skeleton points and the corresponding surface points
  output.open("correspondance.cgal");
  BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton))
    BOOST_FOREACH(vertex_descriptor vd, skeleton[v].vertices)
      output << "2 " << skeleton[v].point << "  " << get(CGAL::vertex_point, tmesh, vd)  << "\n";

  return 0;
}

