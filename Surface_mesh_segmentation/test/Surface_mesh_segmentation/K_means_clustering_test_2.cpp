#include <iostream>
#include <string>
#include <CGAL/internal/Surface_mesh_segmentation/K_means_clustering.h>
typedef CGAL::internal::K_means_clustering K_means;
/**
 * Test where number of points equal to number of centers,
 * Specifically useful for testing random and plus initializations (i.e. Selector class)
 */

int main(void)
{
  for(std::size_t init_type = 0; init_type < 2; ++init_type)
  {
    K_means::Initialization_types init_type_enum = init_type == 0 ? K_means::RANDOM_INITIALIZATION : 
                                                                    K_means::PLUS_INITIALIZATION;
    // Test case: number of points equals to number of clusters
    //            and all points are unique
    for(std::size_t center_size = 1; center_size < 100; ++center_size)
    {
      // unique point generate
      std::vector<double> points;
      for(std::size_t i = 0; i < center_size; ++i) {
        points.push_back(i);
      }

      // test kmeans, expected result: each point has its own cluster
      K_means kmeans(center_size, points, init_type_enum);
      std::vector<std::size_t> center_ids;
      kmeans.fill_with_center_ids(center_ids);

      for(std::size_t i = 0; i < center_size -1; ++i) {
        if(center_ids[i] >= center_ids[i +1]) // center ids are ordered according to mean                                          
        {                                     // and since points are generated ascendingly ... 
          std::string init_type_s = init_type_enum == K_means::RANDOM_INITIALIZATION ? "Random " :
                                                                                       "Plus plus ";
          std::cerr << "Init type: " << init_type_s << "center size: "  << center_size << std::endl;
          CGAL_assertion(false);
        }
      }
    }
  }
}
