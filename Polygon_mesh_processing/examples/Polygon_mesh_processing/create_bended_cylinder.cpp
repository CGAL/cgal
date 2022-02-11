#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

struct Face : public std::array<int,3>
{
  Face(int i, int j, int k)
  {
    (*this)[0] = i;
    (*this)[1] = j;
    (*this)[2] = k;
  }
};

int main()
{
  //create a mesh of a disk
  std::vector<Point_3> circle;
  for (int i=0; i<25; ++i)
    circle.push_back(Point_3(std::cos(2*CGAL_PI/25 * i), std::sin(2*CGAL_PI/25 * i), 0.));
  std::vector<Face> triangles;
  PMP::triangulate_hole_polyline(circle, std::back_inserter(triangles), CGAL::parameters::use_2d_constrained_delaunay_triangulation(true));
  Mesh mesh;
  PMP::polygon_soup_to_polygon_mesh(circle, triangles, mesh);

  std::vector<Point_3> guide;

  for (int i=0; i<12; ++i)
    guide.push_back(Point_3(10*std::cos(2*CGAL_PI/25 * i), 0, 10*std::sin(2*CGAL_PI/25 * i)));

  std::ofstream debug("guide.polylines.txt");
  debug << guide.size();
  for (auto p : guide)
    debug << " " << p;
  debug << std::endl;

  Mesh out;
  PMP::reverse_face_orientations(mesh);
  std::ofstream("in.off") << mesh;

  PMP::sweep_extrude(mesh, guide, out);

  std::ofstream("swept_volume.off") << out;

}