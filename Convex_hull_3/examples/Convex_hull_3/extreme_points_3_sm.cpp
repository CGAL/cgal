#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Extreme_points_traits_adapter_3.h>
#include <CGAL/convex_hull_3.h>
#include <vector>
#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
typedef K::Point_3                                               Point_3;
typedef CGAL::Surface_mesh<Point_3>                              SM;
typedef boost::property_map<SM, CGAL::vertex_point_t>::type      PMAP;
typedef CGAL::Convex_hull_traits_3<K, SM>                        Base_traits;
typedef CGAL::Extreme_points_traits_adapter_3<PMAP, Base_traits> Adapter_traits;

int main(int argc, char* argv[])
{
  std::ifstream in( (argc>1)? argv[1] : "data/cross.off");
  SM sm;
  if(! in || !(in >> sm) || !sm.is_valid()){
    std::cerr<<"Not a valid off file."<<std::endl;
  }
  
  //Define the base ConvexHullTraits
   Base_traits b_traits;
  //get the vertex property map of the Surface_mesh
  PMAP pmap = get(CGAL::vertex_point, sm);
  //This will contain the extreme vertices
  std::vector<boost::graph_traits<SM>::vertex_descriptor> verts;
  //make the adapter traits. You can also directly pass the helper function to extreme_points_3.
  Adapter_traits traits = CGAL::make_extreme_points_traits_adapter(pmap, b_traits);
  
  //call the function with the traits adapter for vertices
  CGAL::extreme_points_3(vertices(sm), std::back_inserter(verts) ,
                   traits);
  //print the number of extreme vertices
  std::cout << "There are  " << verts.size() << " extreme vertices in this mesh." << std::endl;

  return 0;
}
