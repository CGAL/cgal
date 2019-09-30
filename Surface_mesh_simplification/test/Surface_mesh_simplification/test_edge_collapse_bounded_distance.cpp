#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>

//AABB_tree
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

//Timer
#include <CGAL/Timer.h>

//filter{

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface;

typedef CGAL::AABB_face_graph_triangle_primitive<Surface> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

namespace SMS = CGAL::Surface_mesh_simplification ;

namespace CGAL {

namespace Surface_mesh_simplification
{

template<class Placement, class Tree>
class Bounded_distance_placement
{
public:

  typedef typename Placement::TM TM ;

public:

  Bounded_distance_placement(const double sq_dist,
                             const Tree& tree,
                             const Placement& placement = Placement() )
    : mPlacement(placement), tree(tree), threshold_sq_dist(sq_dist)
  {
  }

  template <typename Profile>
  boost::optional<typename Profile::Point>
  operator()( Profile const& aProfile) const
  {
    boost::optional<typename Profile::Point> op = mPlacement(aProfile);
    typedef typename Profile::Point Point;
    if(op){
      const Point* p = boost::get<Point>(&op);
      if(tree.squared_distance(*p) > threshold_sq_dist){
        return boost::none;
      }
    }
    return op;
  }

private:

  Placement  mPlacement ;
  const Tree& tree;
  double threshold_sq_dist;

};
}
}
//} end filter



int main( int argc, char** argv )
{
  Surface surface_mesh;
  std::ifstream is(argc > 1 ? argv[1] : "input.off");
  is >> surface_mesh;
  // This is a stop predicate (defines when the algorithm terminates).
  // In this example, the simplification stops when the number of undirected edges
  // left in the surface mesh drops below the specified number (1000)
  SMS::Count_stop_predicate<Surface> stop(num_halfedges(surface_mesh)/10);

  typedef SMS::Bounded_distance_placement<SMS::LindstromTurk_placement<Surface>, Tree > Filtered_placement;
  typedef SMS::LindstromTurk_placement<Surface> Placement;
  double copy_time, tree_time, placement1_time, placement2_time;

  std::cout<<"input has "<<num_vertices(surface_mesh)<<" vertices."<<std::endl;
  double input_sq_dist=0.000000001;
  std::cout<<"threshold squared distance = "<<input_sq_dist<<std::endl;
  CGAL::Timer timer;
  // This the actual call to the simplification algorithm.
  // The surface mesh and stop conditions are mandatory arguments.
  // The index maps are needed because the vertices and edges
  // of this surface mesh lack an "id()" field.
  timer.start();
  Surface copy_mesh = surface_mesh;
  timer.stop();
  copy_time=timer.time();
  timer.reset();
  timer.start();
  Tree tree( faces(copy_mesh).first, faces(copy_mesh).second, copy_mesh);
  tree.accelerate_distance_queries();
  timer.stop();
  tree_time=timer.time();
  timer.reset();
  Placement placement1;
  Filtered_placement placement2(input_sq_dist, tree, placement1);
  timer.start();

  SMS::edge_collapse( surface_mesh,
                      stop,
                      CGAL::parameters::get_cost (SMS::LindstromTurk_cost<Surface>())
                      .get_placement(placement2)
                      );
  timer.stop();
  placement2_time=timer.time();
  timer.reset();
  timer.start();
  SMS::edge_collapse( copy_mesh,
                      stop,
                      CGAL::parameters::get_cost (SMS::LindstromTurk_cost<Surface>())
                      .get_placement(placement1)
                      );
  timer.stop();
  placement1_time=timer.time();
  surface_mesh.collect_garbage();
  copy_mesh.collect_garbage();
  std::cout<<"Filtered placement total time = "<<copy_time+tree_time+placement2_time<<"s. The placement itself took "<<placement2_time<<"s."<<std::endl;
  std::cout<<" Not filtered placement took "<<placement1_time<<"s."<<std::endl;
  std::cout<<"There are "<<num_vertices(surface_mesh)<<" vertices left when filtered and "<<num_vertices(copy_mesh)<<" when not filtered."<<std::endl;
  std::ofstream os( "out_filtered.off" );
  os.precision(17);
  os << surface_mesh;
  os.close();
  os.open("out_clear.off");
  os << copy_mesh;
  os.close();
  return EXIT_SUCCESS;
}

// EOF //
