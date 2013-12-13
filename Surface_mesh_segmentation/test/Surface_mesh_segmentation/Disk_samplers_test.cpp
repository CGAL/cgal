#include <CGAL/internal/Surface_mesh_segmentation/Disk_samplers.h>

#include <boost/tuple/tuple.hpp>
#include <vector>

typedef boost::tuple<double, double, double> boost_tuple;

void print(const std::vector<boost_tuple>& samples)
{
  const std::size_t map_size = 31;
  const std::size_t map_size_2 = 45;
  std::vector<std::vector<bool> > sample_map(map_size, std::vector<bool>(map_size_2, false));

  for(std::vector<boost_tuple>::const_iterator sample_it = samples.begin(); 
    sample_it != samples.end(); ++sample_it)
  {
    double x = (sample_it->get<0>() +1)/2;
    double y = (sample_it->get<1>() +1)/2;
    x *= (map_size-1);
    y *= (map_size_2-1);
    std::size_t x_c  = static_cast<std::size_t>(x + 0.49);
    std::size_t y_c  = static_cast<std::size_t>(y + 0.49);
    sample_map[x_c][y_c] = true;
  }
  for(std::size_t i = 0; i < map_size; ++i)
  {
    for(std::size_t j = 0; j < map_size_2; ++j)
    {
      if(sample_map[i][j]){ std::cout << "*"; }
      else                { std::cout << " "; }
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}
/**
 * Uses disk sampling functors to sample points from unit-disk.
 * It also prints sampled points for visual debugging.
 *
 * Note that it always return EXIT_SUCCESS
 */
int main(void)
{
  CGAL::internal::Vogel_disk_sampling<boost_tuple> sampling_1;
  CGAL::internal::Vogel_disk_sampling<boost_tuple, true> sampling_2;
  CGAL::internal::Polar_disk_sampling<boost_tuple> sampling_3;
  CGAL::internal::Concentric_disk_sampling<boost_tuple> sampling_4;

  std::vector<boost_tuple> samples_1;
  std::vector<boost_tuple> samples_2;
  std::vector<boost_tuple> samples_3;
  std::vector<boost_tuple> samples_4;

  sampling_1(64, std::back_inserter(samples_1));   
  sampling_2(64, std::back_inserter(samples_2));  
  sampling_3(64, std::back_inserter(samples_3));
  sampling_4(64, std::back_inserter(samples_4));

  print(samples_1);
  print(samples_2);
  print(samples_3);
  print(samples_4);
}
