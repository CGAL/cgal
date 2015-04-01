#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef CGAL::Random_points_in_square_2<Point_2> Random_points_iterator;
typedef Kernel::Circle_2 Sphere_2;
typedef CGAL::Search_traits_2<Kernel> TreeTraits;
typedef CGAL::Euclidean_distance<TreeTraits> Distance;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits, Distance> Neighbor_search;
typedef TreeTraits::Construct_center_d Construct_center_d;
typedef TreeTraits::Compute_squared_radius_d Compute_squared_radius_d;
typedef Neighbor_search::Tree Tree;
int  main() {
  const int N = 1000;
  const unsigned int K = 10;
  Tree tree;
  Random_points_iterator rpg;
  for(int i = 0; i < N; i++){
    tree.insert(*rpg++);
  }
  Sphere_2 query(*rpg,0.5);
  Distance tr_dist;

  Point_2 center = Construct_center_d()(query);
  Neighbor_search N1(tree, center, K, 0.0, false); // eps=0.0, nearest=false
  std::cout << "For the query circle  " << std::endl
        <<  "The " << K << " approximate furthest neighbors are: " << std::endl;
  double radius = Compute_squared_radius_d()(query);
  for (Neighbor_search::iterator it = N1.begin();it != N1.end();it++) {
    std::cout << " Point " << it->first << " at distance = " << tr_dist.inverse_of_transformed_distance(it->second - radius) << std::endl;
  }
  return 0;
}
