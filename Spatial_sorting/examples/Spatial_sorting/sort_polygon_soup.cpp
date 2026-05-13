#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/hilbert_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/property_map.h>
#include <vector>
#include <array>
#include <utility>
#include <string>
#include <filesystem>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Polygon = std::array<std::size_t,3>;
using Polygon_soup = std::pair<std::vector<Point_3>, std::vector<Polygon>>;

typedef CGAL::Spatial_sort_traits_adapter_3<Kernel,
          CGAL::Pointer_property_map<Point_3>::type > Point_traits_3;


template <typename Array, typename Points>
struct Point_of_polygon_property_map
{
  typedef Point_of_polygon_property_map<Array,Points> Self;

  typedef Array key_type; ///< typedef to `Pair`
  typedef typename Points::value_type value_type;
  typedef const value_type& reference;
  typedef boost::lvalue_property_map_tag category;

  Point_of_polygon_property_map(const Points& points)
  : points(points)
  {}

  const value_type& operator[](const key_type& array) const { return points[array[0]]; }

  friend reference get(const Self& pm, const key_type& array) { return pm.points[array[0]]; }

  private:
  const Points& points;
};


typedef CGAL::Spatial_sort_traits_adapter_3<Kernel,
  Point_of_polygon_property_map<Polygon, std::vector<Point_3>>
> Polygon_traits_3;


template <class T>
void permute_in_place(std::vector<T>& a, std::vector<std::size_t>& permutation) {
  const std::size_t n = a.size();

  for (std::size_t i = 0; i < n; ++i) {
    while (permutation[i] != i) {
      std::size_t j = permutation[i];
      std::swap(a[i], a[j]);
      std::swap(permutation[i], permutation[j]);
    }
  }
}


int main(int argc, char* argv[])
{

  const std::string input = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");
  const std::string output = std::filesystem::path(input).stem().string() + "_sorted.off";
  Polygon_soup soup;
  auto& points = soup.first;
  auto& polygons = soup.second;
  CGAL::IO::read_polygon_soup(input, soup.first, soup.second);

  std::vector<std::size_t> indices;
  indices.reserve(points.size());

  std::copy(boost::counting_iterator<std::size_t>(0),
            boost::counting_iterator<std::size_t>(points.size()),
            std::back_inserter(indices));

  CGAL::hilbert_sort( indices.begin(),
                      indices.end(),
                      Point_traits_3(CGAL::make_property_map(points)));

  const std::vector<std::size_t> permutation(indices);

  permute_in_place(points, indices);

  for (Polygon& polygon : polygons) {
    for(int i = 0; i < 3; ++i) {
      polygon[i] = permutation[polygon[i]];
    }
  }

  Point_of_polygon_property_map<Polygon, std::vector<Point_3>> polygon_map(points);
  CGAL::hilbert_sort( polygons.begin(),
                      polygons.end(),
                      Polygon_traits_3(polygon_map));

  CGAL::IO::write_polygon_soup(output, soup.first, soup.second);

  return EXIT_SUCCESS;
}

