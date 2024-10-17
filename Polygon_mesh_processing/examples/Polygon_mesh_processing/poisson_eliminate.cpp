#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>

#include <cyVector.h>
#include <cySampleElim.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;


std::vector< cy::Vec3d > inputPoints, outputPoints;

int main(int argc, char* argv[])
{
  std::string filename = std::filesystem::path(argv[1]).stem().string();
  Mesh sm;
  CGAL::IO::read_polygon_mesh(argv[1], sm);

  CGAL::Bbox_3 bb = CGAL::bbox_3(sm.points().begin(), sm.points().end());
  cy::Vec3d bl(bb.xmin(), bb.ymin(), bb.zmin());
  cy::Vec3d tr(bb.xmax(), bb.ymax(), bb.zmax());
  std::vector<Point_3> points;

  CGAL::Polygon_mesh_processing::sample_triangle_mesh(sm, std::back_inserter(points),
                                                    CGAL::parameters::number_of_points_on_faces(2* num_vertices(sm)).do_sample_vertices(false).do_sample_edges(false));
  double x, y, z;
  std::cout << "# samples = " << points.size() << std::endl;
  double area = CGAL::Polygon_mesh_processing::area(sm);

  std::cout << "area = " << area << std::endl;
  {
    std::string random_points = filename+"-sampled.xyz";
    std::ofstream out(random_points.c_str());
    out.precision(17);
    for (auto p : points) {
        out << p << std::endl;
    }
  }

  for(int i = 0; i < points.size(); ++i){
    inputPoints.push_back(cy::Vec3d(points[i].x(), points[i].y(), points[i].z()));
  }

  outputPoints.resize(num_vertices(sm)/2);

  cy::WeightedSampleElimination< cy::Vec3d, double, 3, int > wse;
  wse.SetBoundsMin(bl);
  wse.SetBoundsMax(tr);
  bool isProgressive = true;

  // independent from CGAL
  std::cout << wse.GetMaxPoissonDiskRadius(2, outputPoints.size(), area) << std::endl;

  double d_max = 2 * wse.GetMaxPoissonDiskRadius( 2, outputPoints.size(), area );

  wse.Eliminate( inputPoints.data(), inputPoints.size(),
                 outputPoints.data(), outputPoints.size(),
                 isProgressive,
                 d_max, 2 );



  std::string poisson_points = filename+"-poisson.xyz";
  std::ofstream out(poisson_points);
  out.precision(17);
  for(int i = 0; i < outputPoints.size(); ++i){
    out << outputPoints[i].x  << " " << outputPoints[i].y << outputPoints[i].z << std::endl;
  }
  return 0;
}
