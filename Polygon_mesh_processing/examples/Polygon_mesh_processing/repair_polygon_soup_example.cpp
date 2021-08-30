#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef K::FT                                                   FT;
typedef K::Point_3                                              Point_3;

typedef CGAL::Surface_mesh<Point_3>                             Mesh;

typedef std::array<FT, 3>                                       Custom_point;
typedef std::vector<std::size_t>                                CGAL_Polygon;

namespace PMP = CGAL::Polygon_mesh_processing;

struct Array_traits
{
  struct Equal_3
  {
    bool operator()(const Custom_point& p, const Custom_point& q) const {
      return (p == q);
    }
  };

  struct Less_xyz_3
  {
    bool operator()(const Custom_point& p, const Custom_point& q) const {
      return std::lexicographical_compare(p.begin(), p.end(), q.begin(), q.end());
    }
  };

  Equal_3 equal_3_object() const { return Equal_3(); }
  Less_xyz_3 less_xyz_3_object() const { return Less_xyz_3(); }
};

int main(int, char**)
{
  // First, construct a polygon soup with some problems
  std::vector<std::array<FT, 3> > points;
  std::vector<CGAL_Polygon> polygons;

  points.push_back(CGAL::make_array<FT>(0,0,0));
  points.push_back(CGAL::make_array<FT>(1,0,0));
  points.push_back(CGAL::make_array<FT>(0,1,0));
  points.push_back(CGAL::make_array<FT>(-1,0,0));
  points.push_back(CGAL::make_array<FT>(0,-1,0));
  points.push_back(CGAL::make_array<FT>(0,1,0)); // duplicate point
  points.push_back(CGAL::make_array<FT>(0,-2,0)); // unused point

  CGAL_Polygon p;
  p.push_back(0); p.push_back(1); p.push_back(2);
  polygons.push_back(p);

  // degenerate face
  p.clear();
  p.push_back(0); p.push_back(0); p.push_back(0);
  polygons.push_back(p);

  p.clear();
  p.push_back(0); p.push_back(1); p.push_back(4);
  polygons.push_back(p);

  // duplicate face with different orientation
  p.clear();
  p.push_back(0); p.push_back(4); p.push_back(1);
  polygons.push_back(p);

  p.clear();
  p.push_back(0); p.push_back(3); p.push_back(5);
  polygons.push_back(p);

  // degenerate face
  p.clear();
  p.push_back(0); p.push_back(3); p.push_back(0);
  polygons.push_back(p);

  p.clear();
  p.push_back(0); p.push_back(3); p.push_back(4);
  polygons.push_back(p);

  // pinched and degenerate face
  p.clear();
  p.push_back(0); p.push_back(1); p.push_back(2); p.push_back(3);
  p.push_back(4); p.push_back(3); p.push_back(2); p.push_back(1);
  polygons.push_back(p);

  std::cout << "Before reparation, the soup has " << points.size() << " vertices and " << polygons.size() << " faces" << std::endl;
  PMP::repair_polygon_soup(points, polygons, CGAL::parameters::geom_traits(Array_traits()));
  std::cout << "After reparation, the soup has " << points.size() << " vertices and " << polygons.size() << " faces" << std::endl;

  Mesh mesh;
  PMP::orient_polygon_soup(points, polygons);
  PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh);

  std::cout << "Mesh has " << num_vertices(mesh) << " vertices and " << num_faces(mesh) << " faces" << std::endl;

  assert(num_vertices(mesh) == 5);
  assert(num_faces(mesh) == 4);

  return 0;
}
