#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Timer.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/boost/graph/property_maps.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <boost/container/flat_map.hpp>

#include <fstream>
#include <ostream>

//typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3                     P_3;
//typedef CGAL::Surface_mesh<K::Point_3> Mesh;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef CGAL::Polyhedron_3<K> Mesh;


int main(int, char**)
{

  Mesh m1,m2;
  std::ifstream in1("/home/mgimeno/Bureau/Data/00_elephant.off");
  in1 >> m1;
  std::ifstream in2("/home/mgimeno/Bureau/Data/02_elephant.off");
  in2 >> m2;

  /*CGAL::make_tetrahedron(P_3(0,1,0),
                         P_3(3,3,0),
                         P_3(2,4,0),
                         P_3(2.5,3.5,1),
                         m1);

  CGAL::make_tetrahedron(P_3(-1,0,0),
                         P_3(-2,-3,0),
                         P_3(-2.5,-3.5,-1),
                         P_3(-3,-4,0),
                         m2);*/

  std::vector<P_3> test_points;
  /*test_points.push_back(P_3(0,0,0));
  test_points.push_back(P_3(1,0,3));
  test_points.push_back(P_3(-1,2,5));
  test_points.push_back(P_3(4.2,1.3,0.6));
  test_points.push_back(P_3(7.1,4.6,.8));
  test_points.push_back(P_3(5.6,6.1,1.7));
  test_points.push_back(P_3(4.6,7.3,2.1));
  test_points.push_back(P_3(1.14,5.20,4.12));
*/
  int size_RU(20000), size_MC(18000);

  CGAL::Timer timer;
  timer.start();
  double min(800), max(0);
  for(int i = 0; i<15; ++i)
  {
    double dist = CGAL::Polygon_mesh_processing::approximated_Hausdorff_distance<CGAL::Sequential_tag>(m1,m2,size_RU);
    if(dist > max) max = dist;
    if(dist < min) min = dist;
  }
  std::cerr<< " UNIFORM 15825: ["<<min<<", "<<max<<"]"<<std::endl;
  std::cerr << " in "
            << timer.time() << " seconds." << std::endl;
  timer.reset();
  min= 800; max = 0;
  for(int i = 0; i<15; ++i)
  {
    double dist = CGAL::Polygon_mesh_processing::approximated_Hausdorff_distance<CGAL::Sequential_tag, K>(m1,m2,0.01, get(CGAL::vertex_point, m1), get(CGAL::vertex_point, m2), CGAL::Polygon_mesh_processing::GRID);
    if(dist > max) max = dist;
    if(dist < min) min = dist;
  }
  std::cerr<< " GRID 0.01: ["<<min<<", "<<max<<"]"<<std::endl;
  std::cerr << " in "
            << timer.time() << " seconds." << std::endl;

  /*timer.reset();
  min= 800; max = 0;
  for(int i = 0; i<15; ++i)
  {
    double dist = CGAL::Polygon_mesh_processing::approximated_Hausdorff_distance<CGAL::Sequential_tag, Mesh>(m1,m2,size_MC, CGAL::Polygon_mesh_processing::MONTE_CARLO);
    if(dist > max) max = dist;
    if(dist < min) min = dist;
  }
  std::cerr<< " MONTE_CARLO 15825: ["<<min<<", "<<max<<"]"<<std::endl;
  std::cerr << " in "
            << timer.time() << " seconds." << std::endl;
  timer.reset();
*/

 /* std::cerr<< " Approximated Hausdorff symetric distance : "<<CGAL::Polygon_mesh_processing::approximated_symmetric_Hausdorff_distance<CGAL::Sequential_tag, Mesh>(m1,m2,40000)<<std::endl;

  std::cerr<< " Approximated Hausdorff distance with named parameters : "<<CGAL::Polygon_mesh_processing::approximated_Hausdorff_distance<CGAL::Sequential_tag>(m1,m2,40000.0)<<std::endl;

  std::cerr<<" Max distance to triangle mesh : "<<CGAL::Polygon_mesh_processing::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(test_points, m1,40000, CGAL::Polygon_mesh_processing::parameters::all_default())<<std::endl;
  */

/*

  PMP::sample_triangle_mesh<K>(m1,size_RU,test_points, PMP::RANDOM_UNIFORM);
  std::cerr<<"RU nb_pts : "<<test_points.size()<<std::endl;
  test_points.clear();
 PMP::sample_triangle_mesh<K>(m1,0.01,test_points, PMP::GRID);
  std::cerr<<"grid nb_pts : "<<test_points.size()<<std::endl;
  test_points.clear();
  PMP::sample_triangle_mesh<K>(m1,size_MC,test_points, PMP::MONTE_CARLO);
  std::cerr<<"MC nb_pts : "<<test_points.size()<<std::endl;
  test_points.clear();
  /**/
}


