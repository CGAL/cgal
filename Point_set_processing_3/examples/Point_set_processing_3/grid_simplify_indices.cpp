#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/property_map.h>

#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

int main(int argc, char*argv[])
{
  // Reads a .xyz point set file in points[].
  std::vector<Point> points;
  std::vector<Vector> normals;
  const char* fname = (argc>1)?argv[1]:"data/fin90_with_PCA_normals.xyz";
  std::ifstream stream(fname);
  Point p;
  Vector v;
  while(stream >> p >> v){
    points.push_back(p);
    normals.push_back(v);
  }

  std::cout << points.size() << " input points" << std::endl;
  std::vector<std::size_t> indices(points.size());
  for(std::size_t i = 0; i < points.size(); ++i){
    indices[i] = i;
  }
  // simplification by clustering using erase-remove idiom
  double cell_size = 0.05;
  std::vector<std::size_t>::iterator end;
  end = CGAL::grid_simplify_point_set(indices,
                                      cell_size,
                                      CGAL::parameters::point_map (CGAL::make_property_map(points)));

  std::size_t k = end - indices.begin();

  std::cerr << "Keep " << k << " of " << indices.size() <<  " indices" << std::endl;

  {
    std::vector<Point> tmp_points(k);
    std::vector<Vector> tmp_normals(k);
    for(std::size_t i=0; i<k; ++i){
      tmp_points[i] = points[indices[i]];
      tmp_normals[i] = normals[indices[i]];
    }
    points.swap(tmp_points);
    normals.swap(tmp_normals);
  }

  std::cout << points.size() << " points after the simplification" << std::endl;

  return EXIT_SUCCESS;
}

