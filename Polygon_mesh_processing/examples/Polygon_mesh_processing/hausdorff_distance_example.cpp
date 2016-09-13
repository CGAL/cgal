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
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
namespace PMP = CGAL::Polygon_mesh_processing;
//typedef CGAL::Polyhedron_3<K> Mesh;


int main(int, char**)
{

  Mesh m1,m2;
  /*std::ifstream in1("/home/mgimeno/Bureau/Data/00_elephant.off");
  in1 >> m1;
  std::ifstream in2("/home/mgimeno/Bureau/Data/02_elephant.off");
  in2 >> m2;*/

  CGAL::make_tetrahedron(P_3(.000,.000,.000),
                         P_3(.002,.000,.000),
                         P_3(.001,.001,.001),
                         P_3(.001,.000,.002),
                         m1);

  CGAL::make_tetrahedron(P_3(-1,0,0),
                         P_3(-2,-3,0),
                         P_3(-2.5,-3.5,-1),
                         P_3(-3,-4,0),
                         m2);

  std::vector<P_3> test_points;
  test_points.push_back(P_3(-1,0,0));
  test_points.push_back(P_3(1,0,3));
  test_points.push_back(P_3(-1,2,5));
  test_points.push_back(P_3(4.2,1.3,0.6));
  test_points.push_back(P_3(7.1,4.6,.8));
  test_points.push_back(P_3(5.6,6.1,1.7));
  test_points.push_back(P_3(4.6,7.3,2.1));
  test_points.push_back(P_3(1.14,5.20,4.12));


  std::cerr<<"distance = "<<PMP::max_distance_to_point_set(m1,test_points,0.01, CGAL::Polygon_mesh_processing::parameters::all_default())<<std::endl;
  std::cerr<<"distance = "<<PMP::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(test_points,m1,CGAL::Polygon_mesh_processing::parameters::all_default())<<std::endl;
  std::cerr<< " Approximated Hausdorff R_U distance : "<<CGAL::Polygon_mesh_processing::approximated_Hausdorff_distance<CGAL::Sequential_tag>(m1,m2,40)<<std::endl;
  std::cerr<< " Approximated Hausdorff GRID distance : "<<CGAL::Polygon_mesh_processing::approximated_Hausdorff_distance<CGAL::Sequential_tag, K>(m1,m2,0.01, get(CGAL::vertex_point, m1), get(CGAL::vertex_point, m2), CGAL::Polygon_mesh_processing::GRID)<<std::endl;
  std::cerr<< " Approximated Hausdorff MONTE_CARLO distance : "<< CGAL::Polygon_mesh_processing::approximated_Hausdorff_distance<CGAL::Sequential_tag, Mesh>(m1,m2,40, CGAL::Polygon_mesh_processing::MONTE_CARLO)<<std::endl;
  std::cerr<< " Approximated Hausdorff symetric distance : "<<CGAL::Polygon_mesh_processing::approximated_symmetric_Hausdorff_distance<CGAL::Sequential_tag>(m1,m2,40)<<std::endl;
  std::cerr<< " Approximated Hausdorff distance with named parameters : "<<CGAL::Polygon_mesh_processing::approximated_Hausdorff_distance<CGAL::Sequential_tag>(m1,m2,40.0)<<std::endl;

}


