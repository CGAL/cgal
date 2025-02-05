#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/poisson_eliminate.h>
#include <CGAL/IO/write_points.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/Real_timer.h>

#include <CGAL/Fuzzy_sphere.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef CGAL::Search_traits_3<K> Traits;
typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_box;

void compare(const char *filename) {
  std::vector<Point_3> points;

  if (!CGAL::IO::read_points(
    filename,
    std::back_inserter(points))) {

    std::cerr << "Error: cannot read file!" << std::endl;
    return;
  }

  //CGAL::Kd_tree<Traits> tree(points.begin(), points.end());

  std::vector<cy::Vec3d> inputPoints, outputPoints;
  std::vector<Point_3> output, out;

  for (int i = 0; i < points.size(); ++i) {
    inputPoints.push_back(cy::Vec3d(CGAL::to_double(points[i].x()), CGAL::to_double(points[i].y()), CGAL::to_double(points[i].z())));
  }

  outputPoints.resize(inputPoints.size() / 5);


  CGAL::Bbox_3 bb = CGAL::bbox_3(points.begin(), points.end());
  cy::Vec3d bl(bb.xmin(), bb.ymin(), bb.zmin());
  cy::Vec3d tr(bb.xmax(), bb.ymax(), bb.zmax());

  cy::WeightedSampleElimination< cy::Vec3d, double, 3, int > wse;
  wse.SetBoundsMin(bl);
  wse.SetBoundsMax(tr);
  bool progressive = false;
  bool weight_limiting = false;
  bool tiling = false;
  wse.SetTiling(tiling);

  output.reserve(outputPoints.size());
  CGAL::Real_timer timer;
  timer.start();
  CGAL::Polygon_mesh_processing::poisson_eliminate2(points, outputPoints.size(), std::back_inserter(output), CGAL::parameters::progressive(progressive).weight_limiting(weight_limiting).tiling(tiling));
  timer.stop();
  std::cout << std::endl << timer.time() << std::endl;
  CGAL::IO::write_points("my_cylinder.xyz", output, CGAL::parameters::stream_precision(17));

  std::cout << points.size() << std::endl;

/*
  timer.reset();
  std::cout << "before" << std::endl;
  timer.start();
  wse.SetWeightLimiting(weight_limiting);
  wse.Eliminate(inputPoints.data(), inputPoints.size(),
    outputPoints.data(), outputPoints.size(), progressive);
  timer.stop();
  std::cout << timer.time() << std::endl;

  for (const cy::Vec3d& p : outputPoints) {
    out.push_back(Point_3(p.x, p.y, p.z));
  }
  CGAL::IO::write_points("orig_cylinder_tiling.xyz", out, CGAL::parameters::stream_precision(17));*/
}


int main(int argc, char* argv[])
{
  std::string filename = std::filesystem::path(argv[1]).stem().string();
  compare(argv[1]);
  /*Mesh sm;
  CGAL::IO::read_polygon_mesh(argv[1], sm);

  std::vector<Point_3> points;

  CGAL::Polygon_mesh_processing::poisson_eliminate(sm, std::back_inserter(points));

  std::string poisson_points = filename+"-poisson.xyz";
  CGAL::IO::write_points(poisson_points, points, CGAL::parameters::stream_precision(17));*/

  return 0;
}
