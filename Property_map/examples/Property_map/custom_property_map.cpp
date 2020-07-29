#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;
using Generator = CGAL::Random_points_on_sphere_3<Point_3>;

// Example of readable property map to get CGAL::Point_3 objects from
// 3 coordinate arrays
struct Custom_point_map
{
  using key_type = std::size_t; // The iterator's value type is an index
  using value_type = Point_3;   // The object manipulated by the algorithm is a Point_3
  using reference = Point_3;    // The object does not exist in memory, so there's no reference
  using category = boost::readable_property_map_tag; // The property map is only used for reading

  double *x, *y, *z;

  Custom_point_map (double* x = nullptr, double* y = nullptr, double* z = nullptr)
    : x(x), y(y), z(z) { }

  // The get() function returns the object expected by the algorithm (here, Point_3)
  friend Point_3 get (const Custom_point_map& map, std::size_t idx)
  {
    return Point_3 (map.x[idx], map.y[idx], map.z[idx]);
  }
};

// Example of read-write property map to get CGAL::Vector_3 objects from
// a buffer array and put CGAL::Vector_3 values in this buffer
struct Custom_normal_map
{
  using key_type = std::size_t; // The iterator's value type is an index
  using value_type = Vector_3;  // The object manipulated by the algorithm is a Vector_3
  using reference = Vector_3;   // The object does not exist in memory, so there's no reference
  using category = boost::read_write_property_map_tag; // The property map is used both
                                                       // for reading and writing data
  double *buffer;

  Custom_normal_map (double* buffer = nullptr)
    : buffer (buffer) { }

  // The get() function returns the object expected by the algorithm (here, Vector_3)
  friend Vector_3 get (const Custom_normal_map& map, std::size_t idx)
  {
    return Vector_3 (map.buffer[idx * 3    ],
                     map.buffer[idx * 3 + 1],
                     map.buffer[idx * 3 + 2]);
  }

  // The put() function updated the user's data structure from the
  // object handled by the algorithm (here Vector_3)
  friend void put (const Custom_normal_map& map, std::size_t idx, const Vector_3& vector_3)
  {
    map.buffer[idx * 3    ] = vector_3.x();
    map.buffer[idx * 3 + 1] = vector_3.y();
    map.buffer[idx * 3 + 2] = vector_3.z();
  }
};


int main()
{
  constexpr std::size_t nb_points = 1000;

  // in this example, points are stored as separate coordinate arrays
  double x[nb_points];
  double y[nb_points];
  double z[nb_points];

  // generate random points
  Generator generator;
  for (std::size_t i = 0; i < nb_points; ++ i)
  {
    Point_3 p = *(generator ++ );
    x[i] = p.x();
    y[i] = p.y();
    z[i] = p.z();
  }

  // normals are stored as a contiguous double array
  double normals[3 *nb_points];

  // we use a vector of indices to access arrays
  std::vector<std::size_t> indices;
  indices.reserve (nb_points);
  for (std::size_t i = 0; i < nb_points; ++ i)
    indices.push_back(i);

  // estimate and orient normals using directly user's data structure
  // instead of creating deep copies using Point_3 and Vector_3
  CGAL::jet_estimate_normals<CGAL::Sequential_tag>
    (indices, 12,
     CGAL::parameters::point_map (Custom_point_map(x,y,z)).
     normal_map (Custom_normal_map(normals)));

  CGAL::mst_orient_normals
    (indices, 12,
     CGAL::parameters::point_map (Custom_point_map(x,y,z)).
     normal_map (Custom_normal_map(normals)));

  // Display first 10 points+normals
  for (std::size_t i = 0; i < 10; ++ i)
    std::cerr << "Point(" << i << ") = " << x[i] << " " << y[i] << " " << z[i]
              << "\tNormal(" << i << ") = "
              << normals[3*i] << " " << normals[3*i+1] << " " << normals[3*i+2] << std::endl;

  return EXIT_SUCCESS;
}
