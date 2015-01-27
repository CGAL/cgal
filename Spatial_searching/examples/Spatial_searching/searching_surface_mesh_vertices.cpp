#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Surface_mesh.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;

typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor Point;
typedef boost::graph_traits<Mesh>::vertices_size_type size_type;

typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type Vertex_point_pmap;

typedef CGAL::Search_traits_3<Kernel>                                    Traits_base;
typedef CGAL::Search_traits_adapter<Point,Vertex_point_pmap,Traits_base> Traits;


typedef CGAL::Orthogonal_k_neighbor_search<Traits>                      K_neighbor_search;
typedef K_neighbor_search::Tree                                         Tree;
typedef Tree::Splitter                                                  Splitter;
typedef K_neighbor_search::Distance                                     Distance;

int main(int argc, char* argv[]) {
  Mesh mesh;
  std::cerr << argc <<std::endl;
  std::ifstream in((argc>1)?argv[1]:"data/tripod.off");
  in  >> mesh;
  const unsigned int K = 5;

  Vertex_point_pmap vppmap = get(CGAL::vertex_point,mesh);

  // Insert number_of_data_points in the tree
  Tree tree(
            vertices(mesh).begin(),
            vertices(mesh).end(),
            Splitter(),
            Traits(vppmap)
  );
  Point_3 query(0.0, 0.0, 0.0);
  Distance tr_dist(vppmap);

  // search K nearest neighbours
  K_neighbor_search search(tree, query, K,0,true,tr_dist);
  std::cout <<"The "<< K << " nearest vertices to the query point at (0,0,0) are:" << std::endl;
  for(K_neighbor_search::iterator it = search.begin(); it != search.end(); it++){
    std::cout << "vertex " << it->first << " : " << vppmap[it->first] << " at distance " 
	      << tr_dist.inverse_of_transformed_distance(it->second) << std::endl;
  }
  return 0;
}
