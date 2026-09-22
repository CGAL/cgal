#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/poisson_eliminate.h>
#include <CGAL/IO/write_points.h>
#include <CGAL/IO/read_points.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;

void sampling(const std::string& filename, double size_factor = 0.2) {
  if (size_factor >= 1.0) {
    std::cout << "usage poisson_eliminate_example filename size_factor" << std::endl
              << "0 < size_factor < 1" << std::endl;
    return;
  }
  std::vector<Point_3> points;

  if (!CGAL::IO::read_points(
    filename,
    std::back_inserter(points))) {

    std::cerr << "Error: cannot read file!" << std::endl;
    return;
  }

  std::size_t target_size = std::size_t(points.size() * size_factor);
  std::vector<Point_3> output;
  output.reserve(target_size);

  CGAL::poisson_eliminate(points, target_size, std::back_inserter(output));

  CGAL::IO::write_points("out.xyz", output, CGAL::parameters::stream_precision(17));
}


int main(int argc, char* argv[])
{
  if (argc < 2)
    sampling(CGAL::data_file_path("points_3/radar.xyz"));
  else if (argc < 3)
    sampling(argv[1]);
  else if (argc < 4)
    sampling(argv[1], std::atof(argv[2]));

  return 0;
}
