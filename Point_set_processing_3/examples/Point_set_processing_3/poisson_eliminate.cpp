#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/poisson_eliminate.h>
#include <CGAL/IO/write_points.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/Real_timer.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;

void sampling(const std::string& filename) {
  std::vector<Point_3> points;

  if (!CGAL::IO::read_points(
    filename,
    std::back_inserter(points))) {

    std::cerr << "Error: cannot read file!" << std::endl;
    return;
  }

  std::size_t target_size = points.size() / 5;

  bool progressive = true;
  bool weight_limiting = false;
  bool tiling = false;
  std::vector<Point_3> output;

#ifdef CGAL_USE_CY
  std::vector<cy::Vec3d> inputPoints, outputPoints;
  std::vector<Point_3> out;

  for (int i = 0; i < points.size(); ++i)
    inputPoints.push_back(cy::Vec3d(CGAL::to_double(points[i].x()), CGAL::to_double(points[i].y()), CGAL::to_double(points[i].z())));

  outputPoints.resize(target_size);

  CGAL::Bbox_3 bb = CGAL::bbox_3(points.begin(), points.end());
  cy::Vec3d bl(bb.xmin(), bb.ymin(), bb.zmin());
  cy::Vec3d tr(bb.xmax(), bb.ymax(), bb.zmax());

  cy::WeightedSampleElimination< cy::Vec3d, double, 3, int > wse;
  wse.SetBoundsMin(bl);
  wse.SetBoundsMax(tr);
  wse.SetTiling(tiling);
#endif

  output.reserve(target_size);
  CGAL::Real_timer timer;
  timer.start();
  CGAL::poisson_eliminate(points, target_size, std::back_inserter(output), CGAL::parameters::progressive(progressive).weight_limiting(weight_limiting).tiling(tiling));
  timer.stop();
  std::cout << std::endl << timer.time() << std::endl;
  CGAL::IO::write_points("out.xyz", output, CGAL::parameters::stream_precision(17));

#ifdef CGAL_USE_CY
  timer.reset();
  timer.start();
  wse.SetWeightLimiting(weight_limiting);
  wse.Eliminate(inputPoints.data(), int(inputPoints.size()),
    outputPoints.data(), int(outputPoints.size()), progressive);
  timer.stop();
  std::cout << timer.time() << std::endl;

  for (const cy::Vec3d& p : outputPoints) {
    out.push_back(Point_3(p.x, p.y, p.z));
  }
  CGAL::IO::write_points("out_ref.xyz", out, CGAL::parameters::stream_precision(17));
#endif
}


int main(int argc, char* argv[])
{
  if (argc < 2)
    sampling(CGAL::data_file_path("points_3/radar.xyz"));
  else
    sampling(argv[1]);

  return 0;
}
