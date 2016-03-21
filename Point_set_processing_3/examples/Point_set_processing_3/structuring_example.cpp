#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>

#include <CGAL/Shape_detection_3.h>

#include <CGAL/structure_point_set.h>

#include <iostream>
#include <fstream>

// Type declarations
typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_3                                      Point;
typedef std::pair<Kernel::Point_3, Kernel::Vector_3>         Point_with_normal;
typedef std::vector<Point_with_normal>                       Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

// Efficient RANSAC types
typedef CGAL::Shape_detection_3::Efficient_RANSAC_traits
  <Kernel, Pwn_vector, Point_map, Normal_map>                Traits;
typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits>    Efficient_ransac;
typedef CGAL::Shape_detection_3::Plane<Traits>               Plane;

typedef std::vector<boost::shared_ptr<Plane> >               Adjacency;

int main (int argc, char** argv)
{
  // Points with normals.
  Pwn_vector points;

  // Loading point set from a file. 
  std::ifstream stream(argc>1 ? argv[1] : "data/cube.pwn");

  if (!stream || 
    !CGAL::read_xyz_points_and_normals(stream,
      std::back_inserter(points),
      Point_map(),
      Normal_map()))
  {
      std::cerr << "Error: cannot read file cube.pwn" << std::endl;
      return EXIT_FAILURE;
  }

  std::cerr << points.size() << " point(s) read." << std::endl;
  
  // Shape detection
  Efficient_ransac ransac;
  ransac.set_input(points);
  ransac.add_shape_factory<Plane>();
  ransac.detect();

  Pwn_vector structured_pts;
  std::vector<Adjacency> adjacency;

  CGAL::structure_point_set (points.begin (), points.end (), // input points
                             std::back_inserter (structured_pts),
                             std::back_inserter (adjacency),
                             ransac, // shape detection engine
                             0.1); // epsilon for structuring points

  std::cerr << structured_pts.size ()
            << " structured point(s) generated." << std::endl;

  CGAL::cpp11::array<std::size_t, 4> nb_pts = {{ 0, 0, 0, 0 }};

  // Check degree of points and display statistics
  for (std::size_t i = 0; i < adjacency.size(); ++ i)
    nb_pts[adjacency[i].size()] ++;

  std::cerr << " * " << nb_pts[0] << " crease point(s) (unstructured)" << std::endl
            << " * " << nb_pts[1] << " planar point(s)" << std::endl
            << " * " << nb_pts[2] << " edge point(s)" << std::endl
            << " * " << nb_pts[3] << " corner point(s)" << std::endl;
    
  return EXIT_SUCCESS;
}
