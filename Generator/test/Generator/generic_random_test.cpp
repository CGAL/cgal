#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

#include <CGAL/Polygon_2.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/point_generators_3.h>

#include <iostream>
#include <iterator>
#include <vector>

using namespace CGAL::parameters;

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;

typedef K::FT                                                     FT;
typedef K::Point_2                                                Point_2;
typedef K::Point_3                                                Point_3;
typedef K::Weighted_point_3                                       Weighted_point_3;
typedef K::Triangle_2                                             Triangle_2;
typedef K::Triangle_3                                             Triangle_3;

typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K>                        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>              Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>        CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>                  Mesh_2_criteria;

typedef CGAL::Polygon_2<K>                                        Polygon_2;

typedef CGAL::Mesh_polyhedron_3<K>::type                          Polyhedron;

// Domain
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K>           Mesh_domain;
typedef Mesh_domain::Corner_index                                 Corner_index;
typedef Mesh_domain::Curve_index                                  Curve_index;

typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type             Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
          Tr, Corner_index, Curve_index>                          C3t3;
typedef CGAL::Mesh_criteria_3<Tr>                                 Mesh_criteria;
typedef C3t3::Point                                               Point;

int test_triangles_2(const FT eps)
{
  std::cout << "test_triangles_2 (tolerance: " << eps << ")" << std::endl;

  // Generated points are in that vector
  std::vector<Point_2> points;

  // Create input triangles
  std::vector<Triangle_2> triangles;
  triangles.push_back(Triangle_2(Point_2(0,0), Point_2(0.5,0), Point_2(0,0.5)));
  triangles.push_back(Triangle_2(Point_2(0,0.5), Point_2(0.5,0), Point_2(0.5,0.5)));

  // Create the generator, input is the vector of Triangle_2
  CGAL::Random_points_in_triangles_2<Point_2> g(triangles);

  // Get 100 random points in triangle range
  std::copy_n(g, 100, std::back_inserter(points));

  // Check that we have really created 100 points.
  assert( points.size() == 100);

  for(Point_2 p : points)
  {
    bool on_quad = p.x() > -eps && p.x() < 0.5 + eps &&
                   p.y() > -eps && p.y() < 0.5 + eps;
    if(!on_quad)
    {
      std::cerr << "ERROR : Generated point (" << p << ") is not on the square." << std::endl;
      return 0;
    }
  }

  return 1;
}

int test_triangles_3(const FT eps)
{
  std::cout << "test_triangles_3 (tolerance: " << eps << ")" << std::endl;

  // Generated points are in that vector
  std::vector<Point_3> points;

  // Create input triangles
  std::vector<Triangle_3> triangles;
  triangles.push_back(Triangle_3(Point_3(0,0,0), Point_3(0.5,0,0), Point_3(0,0.5,0)));
  triangles.push_back(Triangle_3(Point_3(0,0.5,0), Point_3(0.5,0,0), Point_3(0.5,0.5,0)));
  triangles.push_back(Triangle_3(Point_3(0.5,0,0), Point_3(0.5,0,0.5), Point_3(0.5,0.5,0)));
  triangles.push_back(Triangle_3(Point_3(0.5,0.5,0), Point_3(0.5,0.5,0.5), Point_3(0.5,0.,0.5)));

  // Create the generator, input is the vector of Triangle_3
  CGAL::Random_points_in_triangles_3<Point_3> g(triangles);

  // Get 100 random points in triangle range
  std::copy_n(g, 100, std::back_inserter(points));

  // Check that we have really created 100 points.
  assert(points.size() == 100);

  for(Point_3 p : points)
  {
    bool on_front = p.z() <  eps && p.z() > -eps &&
                    p.x() > -eps && p.x() < 0.5 + eps &&
                    p.y() > -eps && p.y() < 0.5 + eps;
    bool on_right = p.x() < 0.5 + eps && p.x() > 0.5 - eps &&
                    p.z() > -eps && p.z() < 0.5 + eps &&
                    p.y() > -eps && p.y() < 0.5 + eps;

    if(!on_front && !on_right)
    {
      std::cerr << "ERROR : Generated point (" << p << ") is not on a triangle of the range." << std::endl;
      return 0;
    }
  }

  return 1;
}

int test_T2(const FT eps)
{
  std::cout << "test_T2 (tolerance: " << eps << ")" << std::endl;

  std::vector<Point_2> points;

  Polygon_2 polygon1;

  polygon1.push_back(Point_2(0,0));
  polygon1.push_back(Point_2(2,0));
  polygon1.push_back(Point_2(2,2));
  polygon1.push_back(Point_2(0,2));

  // Insert the constraint...
  CDT cdt;
  cdt.insert_constraint(polygon1.vertices_begin(), polygon1.vertices_end(), true);

  // ... and four points outside the box to create faces outside the domain
  cdt.insert(Point_2(-1,-1));
  cdt.insert(Point_2(3,-1));
  cdt.insert(Point_2(3,3));
  cdt.insert(Point_2(-1,3));

  // Refine the triangulation
  CGAL::refine_Delaunay_mesh_2(cdt, Mesh_2_criteria(0.125, 0.5));

  CGAL::Random_points_in_triangle_mesh_2<Point_2, CDT> g(cdt);
  std::copy_n(g, 300, std::back_inserter(points));
  for(std::size_t i=0; i<points.size(); ++i)
  {
    Point_2 p = points[i];
    for(int j=0; j<2; ++j)
    {
      double coords[2] = {p.x(), p.y()};
      if(coords[j] > 2.0 + eps || coords[j] < -eps)
      {
        std::cerr << "ERROR : Generated point (" << p << ") is not on a face of the domain." << std::endl;
        return 0;
      }
    }
  }

  return 1;
}

bool on_face(int face, double coord[3], const FT eps)
{
  if(CGAL::abs(CGAL::abs(coord[face]) - 0.5) < eps &&
     CGAL::abs(coord[(face+1)%3]) - 0.5 < eps &&
     CGAL::abs(coord[(face+2)%3]) - 0.5 < eps)
    return true;

  return false;
}

int test_volume_mesh(Polyhedron& polyhedron, const FT eps)
{
  std::cout << "test_volume_mesh (tolerance: " << eps << ")" << std::endl;

  std::vector<Point_3> points;
  CGAL::Random_points_in_triangle_mesh_3<Polyhedron> g(polyhedron);
  std::copy_n(g, 300, std::back_inserter(points));
  for(std::size_t i=0; i<points.size(); ++i)
  {
    Point_3 p = points[i];
    double coords[3] = {p.x(), p.y(), p.z()};
    if(!(on_face(0, coords, eps) || on_face(1, coords, eps) || on_face(2, coords, eps)))
    {
      std::cerr << "ERROR : Generated point (" << p << ") is not on a face." << std::endl;
      return 0;
    }
  }

  return 1;
}

int test_on_c3t3(const Polyhedron& polyhedron, const FT eps)
{
  std::cout << "test_on_c3t3 (tolerance: " << eps << ")" << std::endl;

  std::vector<Point_3> points;

  // Create domain
  Mesh_domain domain(polyhedron);
  domain.detect_features();

  // Mesh criteria (no cell_size set)
  Mesh_criteria criteria(facet_angle=25,
                         facet_size=0.05,
                         facet_distance=0.008,
                         cell_radius_edge_ratio=3);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

  CGAL::Random_points_in_tetrahedral_mesh_boundary_3<C3t3> g(c3t3);
  std::copy_n(g, 300, std::back_inserter(points));
  for(std::size_t i=0; i<points.size(); ++i)
  {
    Point_3 p = points[i];
    double coords[3] = {p.x(), p.y(), p.z()};
    if(!(on_face(0, coords, eps) || on_face(1, coords, eps) || on_face(2, coords, eps)))
    {
      std::cerr << "ERROR : Generated point (" << p << ") is not on a face." << std::endl;
      return 0;
    }
  }

  return 1;
}

int test_in_c3t3(const Polyhedron& polyhedron, const FT eps)
{
  std::cout << "test_in_c3t3 (tolerance: " << eps << ")" << std::endl;

  std::vector<Point_3> points;

  // Create domain
  Mesh_domain domain(polyhedron);
  domain.detect_features();

  // Mesh criteria (no cell_size set)
  Mesh_criteria criteria(facet_angle=25,
                         facet_size=0.05,
                         facet_distance=0.008,
                         cell_radius_edge_ratio=3);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

  CGAL::Random_points_in_tetrahedral_mesh_3<C3t3> g(c3t3);
  std::copy_n(g, 300, std::back_inserter(points));
  for(std::size_t i=0; i<points.size(); ++i)
  {
    Point_3 p = points[i];
    double coords[3] = {p.x(), p.y(), p.z()};
    for(int j=0; j<3; ++j)
    {
      if(CGAL::abs(coords[j]) > 0.5 + eps)
      {
        std::cerr << "ERROR : Generated point (" << p << ") is not on the cube." << std::endl;
        return 0;
      }
    }
  }

  return 1;
}

int main()
{
  Polyhedron polyhedron;

  // A cube
  make_hexahedron(
        Point_3(-0.5,-0.5,-0.5), Point_3(0.5,-0.5,-0.5), Point_3(0.5,0.5,-0.5),
        Point_3(-0.5,0.5,-0.5), Point_3(-0.5,0.5,0.5), Point_3(-0.5,-0.5,0.5),
        Point_3(0.5,-0.5,0.5), Point_3(0.5,0.5,0.5),
                  polyhedron);

  boost::graph_traits<Polyhedron>::halfedge_descriptor facets[6];
  int i = 0;
  for(boost::graph_traits<Polyhedron>::face_descriptor fd : faces(polyhedron))
    facets[i++] = halfedge(fd, polyhedron);

  for(int i=0; i<6; ++i)
    CGAL::Euler::split_face(facets[i],next(next(facets[i], polyhedron), polyhedron), polyhedron);

  const FT eps = 1e-10;

  int validity =
      test_triangles_2(eps)
      *test_triangles_3(eps)
      *test_volume_mesh(polyhedron, eps)
      *test_T2(eps)
      *test_on_c3t3(polyhedron, eps)
      *test_in_c3t3(polyhedron, eps)
      ;
  assert(validity == 1);

  return 0;
}
