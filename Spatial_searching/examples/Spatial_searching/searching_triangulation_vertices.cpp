#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/function_objects.h>

#include <boost/property_map/function_property_map.hpp>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Delaunay_triangulation_3<Kernel> Delaunay;
typedef Delaunay::Vertex_handle Vertex_handle;


struct Project {
  typedef Vertex_handle   argument_type;
  typedef const Point_3&  Point;
  typedef Point           result_type;
  Point operator()( Vertex_handle  v) const { return v->point(); }
};


typedef boost::function_property_map<Project, Vertex_handle> Vertex_point_pmap;
typedef CGAL::Search_traits_3<Kernel>                        Traits_base;
typedef CGAL::Search_traits_adapter<Vertex_handle,Vertex_point_pmap,Traits_base> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits>                               K_neighbor_search;
typedef K_neighbor_search::Tree                                                  Tree;
typedef Tree::Splitter                                                           Splitter;
typedef K_neighbor_search::Distance                                              Distance;

int main(int argc, char* argv[])
{
  std::ifstream in((argc>1)?argv[1]:CGAL::data_file_path("points_3/kitten.xyz"));
  std::vector<Point_3> points;
  Point_3 p;
  while(in >> p){
    points.push_back(p);
  }

  Delaunay dt(points.begin(), points.end());

  const unsigned int K = 5;

  Project project;
  Vertex_point_pmap vppmap(project);

  // Insert number_of_data_points in the tree
  Tree tree(dt.finite_vertex_handles().begin(), dt.finite_vertex_handles().end(), Splitter(), Traits(vppmap));

  // search K nearest neighbors
  Point_3 query(0.0, 0.0, 0.0);
  Distance tr_dist(vppmap);

  K_neighbor_search search(tree, query, K,0,true,tr_dist);
  std::cout <<"The "<< K << " nearest vertices to the query point at (0,0,0) are:" << std::endl;
  for(K_neighbor_search::iterator it = search.begin(); it != search.end(); it++){
    std::cout << "vertex " << &*(it->first) << " : " << vppmap[it->first] << " at distance "
              << tr_dist.inverse_of_transformed_distance(it->second) << std::endl;
  }

  return 0;
}
