#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/triangulate_hole.h>

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <CGAL/utility.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;

struct My_triangle {
  int v0, v1, v2;
  My_triangle(int v0, int v1, int v2) : v0(v0), v1(v1), v2(v2) { }
};

int main() {
  std::vector<Point_3> polyline;
  polyline.push_back(Point_3( 1.,0.,0.));
  polyline.push_back(Point_3( 0.,1.,0.));
  polyline.push_back(Point_3(-1.,0.,0.));
  polyline.push_back(Point_3( 1.,1.,0.));
  // repeating first point (i.e. polyline.push_back(Point_3(1.,0.,0.)) ) is optional

  // any type, having Type(int, int, int) constructor available, can be used to hold output triangles
  std::vector<boost::tuple<int, int, int> > patch_1;
  std::vector<CGAL::Triple<int, int, int> > patch_2;
  std::vector<My_triangle>                  patch_3;

  patch_1.reserve(polyline.size() -2); // there will be exactly n-2 triangles in the patch
  CGAL::Polygon_mesh_processing::triangulate_hole_polyline(
    polyline.begin(), polyline.end(), back_inserter(patch_1));
  CGAL::Polygon_mesh_processing::triangulate_hole_polyline(
    polyline.begin(), polyline.end(), back_inserter(patch_2));
  CGAL::Polygon_mesh_processing::triangulate_hole_polyline(
    polyline.begin(), polyline.end(), back_inserter(patch_3));

  for(std::size_t i = 0; i < patch_1.size(); ++i) {
    std::cout << "Triangle " << i << ": " << patch_1[i].get<0>() << " " 
              << patch_1[i].get<1>() << " " << patch_1[i].get<2>() << std::endl;

    CGAL_assertion(patch_1[i].get<0>() == patch_2[i].first
                && patch_2[i].first  == patch_3[i].v0);
    CGAL_assertion(patch_1[i].get<1>() == patch_2[i].second
                && patch_2[i].second == patch_3[i].v1);
    CGAL_assertion(patch_1[i].get<2>() == patch_2[i].third
                && patch_2[i].third == patch_3[i].v2);
  }

  // note that no degenerate triangle is constructed in patch
  std::vector<Point_3> polyline_collinear;
  polyline_collinear.push_back(Point_3(1.,0.,0.));
  polyline_collinear.push_back(Point_3(2.,0.,0.));
  polyline_collinear.push_back(Point_3(3.,0.,0.));
  polyline_collinear.push_back(Point_3(4.,0.,0.));
  std::vector<My_triangle> patch_will_be_empty;
  CGAL::Polygon_mesh_processing::triangulate_hole_polyline(
    polyline_collinear.begin(), polyline_collinear.end(), 
    back_inserter(patch_will_be_empty));
  CGAL_assertion(patch_will_be_empty.empty());
}
