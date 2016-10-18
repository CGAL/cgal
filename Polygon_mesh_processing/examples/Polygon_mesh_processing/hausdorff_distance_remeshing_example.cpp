#include <CGAL/Surface_mesh.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>


#if CGAL_LINKED_WITH_TBB
#define TAG CGAL::Parallel_tag
#else
#define TAG CGAL::Sequential_tag
#endif
typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3                     P_3;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
namespace PMP = CGAL::Polygon_mesh_processing;


int main(int, char**)
{

  Mesh m1, m2;
  CGAL::make_tetrahedron(P_3(.0,.0,.0),
                         P_3(2,.0,.0),
                         P_3(1,1,1),
                         P_3(1,.0,2),
                         m1);
  m2=m1;
  PMP::isotropic_remeshing(m2.faces(),.05, m2);
  std::cerr<< " Approximated Hausdorff distance : "<< CGAL::Polygon_mesh_processing::approximated_Hausdorff_distance<TAG>(m1,m2,4000)<<std::endl;

}


