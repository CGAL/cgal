#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/alpha_wrap_3.h>

#include <CGAL/IO/polygon_soup_io.h>

#include <array>
#include <iostream>
#include <string>
#include <vector>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;

using Mesh = CGAL::Surface_mesh<Point_3>;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  const int argc_check = argc - 1;
  const char* entry_name_ptr = nullptr;
  double relative_alpha_ratio = 20., relative_offset_ratio = 600.;

  for(int i=1; i<argc; ++i)
  {
    if(!strcmp("-i", argv[i]) && i < argc_check)
      entry_name_ptr = argv[++i];
    else if(!strcmp("-a", argv[i]) && i < argc_check)
      relative_alpha_ratio = std::stod(argv[++i]);
    else if(!strcmp("-d", argv[i]) && i < argc_check)
      relative_offset_ratio = std::stod(argv[++i]);
  }

  if(argc < 3 || relative_alpha_ratio <= 0.)
  {
    std::cerr << "Error: bad input parameters." << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<Point_3> points;
  std::vector<std::array<std::size_t, 3> > faces;
  if(!CGAL::IO::read_polygon_soup(entry_name_ptr, points, faces) || faces.empty())
  {
    std::cerr << "Error: Invalid input data." << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::Bbox_3 bbox;
  for(const Point_3& p : points)
    bbox += p.bbox();

  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));
  const double alpha = diag_length / relative_alpha_ratio;
  const double offset = diag_length / relative_offset_ratio;

  Mesh wrap;
  CGAL::alpha_wrap_3(points, faces, alpha, offset, wrap);

  return EXIT_SUCCESS;
}
