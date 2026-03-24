#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Simple_cartesian.h>

#include <vector>
#include <fstream>

int main(int argc, char** argv)
{
  std::vector<CGAL::Simple_cartesian<double>::Point_3> points;
  std::vector<std::vector<std::size_t>> polygons;

  if (argc!=3)
  {
    std::cerr << "Usage: " << argv[0] << " input.XXX output.YYY\n";
    std::cerr << "See https://doc.cgal.org/latest/Stream_support/index.html#title11 for the list of supported file formats.\n";
    return 0;
  }

  if (!CGAL::IO::read_polygon_soup(argv[1], points, polygons))
  {
    std::cerr << "Error reading " << argv[1] << "\n";
    return EXIT_FAILURE;
  }
  if (!CGAL::IO::write_polygon_soup(argv[2], points, polygons, CGAL::parameters::stream_precision(17)))
  {
    std::cerr << "Error writing " << argv[2] << "\n";
    return EXIT_FAILURE;
  }

  return 0;
}
