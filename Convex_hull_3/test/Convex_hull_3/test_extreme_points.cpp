#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Vertex_to_point_traits_adapter_3.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <vector>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <iostream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;
typedef K::Point_3 Point_3;


void test_function_overload()
{
  std::vector<Point_3> points;
  points.push_back(Point_3(0,0,0));
  points.push_back(Point_3(10,0,0));
  points.push_back(Point_3(0,10,0));
  points.push_back(Point_3(0,0,10));
  points.push_back(Point_3(5,5,5));
  points.push_back(Point_3(2,5,3));
  points.push_back(Point_3(1,3,2));

  std::vector<Point_3> extreme_points;
  CGAL::extreme_points_3(points, std::back_inserter(extreme_points));
  CGAL_assertion(extreme_points.size() == 5);

  Polyhedron_3 polyhedron;
  CGAL::convex_hull_3(points.begin(), points.end(), polyhedron);
  typedef Polyhedron_3::Point_iterator Point_iterator;
  for (Point_iterator p_it = polyhedron.points_begin(); p_it != polyhedron.points_end(); ++p_it)
  {
    Point_3 p = *p_it;
    CGAL_assertion(std::find(extreme_points.begin(), extreme_points.end(), p) != extreme_points.end());
  }
}

void test_triangulated_cube(const char* fname)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::Point_3 Point_3;
  typedef CGAL::Surface_mesh<Point_3> SurfaceMesh;

  std::ifstream input(fname);
  SurfaceMesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  std::vector<Point_3> mesh_points;
  typedef typename boost::property_map<SurfaceMesh, boost::vertex_point_t>::type Pmap;
  Pmap vpmap = get_property_map(boost::vertex_point, mesh);

  typedef typename boost::graph_traits<SurfaceMesh>::vertex_descriptor vertex_descriptor;
  BOOST_FOREACH(vertex_descriptor v, vertices(mesh))
  {
    Point_3 p = get(vpmap, v);
    mesh_points.push_back(p);
  }

  std::vector<Point_3> extreme_points;
  CGAL::extreme_points_3(mesh_points, std::back_inserter(extreme_points)) ;
  CGAL_assertion(extreme_points.size() == 8);
}

void test_coplanar_points(const char* fname)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::Point_3 Point_3;

  std::ifstream input(fname);
  std::vector<Point_3> points;
  if(input)
  {
    Point_3 p;
    int i = 0;
    while(input >> p)
    {
      if(i % 2 == 0) // avoid normals in .xyz file
      {
        points.push_back(p);
      }
      ++i;
    }
  }
  else
    std::cerr << "error loading file\n";

  CGAL_assertion(points.size() == 289);
  std::vector<Point_3> convex_hull;
  CGAL::extreme_points_3(points, std::back_inserter(convex_hull));
  CGAL_assertion(convex_hull.size() == 4);
}

void test_3_points()
{
  std::vector<Point_3> points;
  points.push_back(Point_3(0,0,0));
  points.push_back(Point_3(10,0,0));
  points.push_back(Point_3(0,10,0));

  Polyhedron_3 polyhedron;
  CGAL::convex_hull_3(points.begin(), points.end(), polyhedron);
  std::vector<Point_3> convex_polyhedron(polyhedron.points_begin(), polyhedron.points_end());
  CGAL_assertion(convex_polyhedron.size() == 3);

  std::vector<Point_3> extreme_points;
  CGAL::extreme_points_3(points, std::back_inserter(extreme_points));
  CGAL_assertion(convex_polyhedron.size() == 3);

  typedef Polyhedron_3::Point_iterator Point_iterator;
  for (Point_iterator p_it = polyhedron.points_begin(); p_it != polyhedron.points_end(); ++p_it)
  {
    Point_3 p = *p_it;
    CGAL_assertion(std::find(extreme_points.begin(), extreme_points.end(), p) != extreme_points.end());
  }
}

void test_3_collinear()
{
  std::vector<Point_3> points;
  points.push_back(Point_3(0,0,0));
  points.push_back(Point_3(1,0,0));
  points.push_back(Point_3(2,0,0));

  Polyhedron_3 polyhedron;
  CGAL::convex_hull_3(points.begin(), points.end(), polyhedron);
  std::vector<Point_3> convex_polyhedron(polyhedron.points_begin(), polyhedron.points_end());
  CGAL_assertion(convex_polyhedron.size() == 2);

  std::vector<Point_3> extreme_points;
  CGAL::extreme_points_3(points, std::back_inserter(extreme_points));
  CGAL_assertion(convex_polyhedron.size() == 2);

  typedef Polyhedron_3::Point_iterator Point_iterator;
  for (Point_iterator p_it = polyhedron.points_begin(); p_it != polyhedron.points_end(); ++p_it)
  {
    Point_3 p = *p_it;
    CGAL_assertion(std::find(extreme_points.begin(), extreme_points.end(), p) != extreme_points.end());
  }
}

void test_up_to_3_extreme_points()
{
  std::vector<Point_3> points;
  std::vector<Point_3> extreme_points;
  CGAL::extreme_points_3(points, std::back_inserter(extreme_points));
  CGAL_assertion(extreme_points.empty());

  Point_3 p1(0, 0, 0);
  points.push_back(p1);
  extreme_points.clear();
  CGAL::extreme_points_3(points, std::back_inserter(extreme_points));
  CGAL_assertion(extreme_points.size() == 1);

  Point_3 p2(1, 0, 0);
  points.push_back(p2);
  extreme_points.clear();
  CGAL::extreme_points_3(points, std::back_inserter(extreme_points));
  CGAL_assertion(extreme_points.size() == 2);

  Point_3 p3(1, 1, 0);
  points.push_back(p3);
  extreme_points.clear();
  CGAL::extreme_points_3(points, std::back_inserter(extreme_points));
  CGAL_assertion(extreme_points.size() == 3);
}

void test_equal_points()
{
  std::vector<Point_3> points;
  std::vector<Point_3> extreme_points;

  // test two equal
  Point_3 p1(0, 0, 0);
  points.push_back(p1);
  points.push_back(p1);
  CGAL::extreme_points_3(points, std::back_inserter(extreme_points));
  CGAL_assertion(extreme_points.size() == 1);

  // test many equal
  extreme_points.clear();
  for(int i = 0; i < 10; ++i)
    points.push_back(p1);
  CGAL::extreme_points_3(points, std::back_inserter(extreme_points));
  CGAL_assertion(extreme_points.size() == 1);

  // test with only 2 different
  extreme_points.clear();
  Point_3 p3(0.1, 0, 0);
  points.push_back(p3);
  CGAL::extreme_points_3(points, std::back_inserter(extreme_points));
  CGAL_assertion(extreme_points.size() == 2);
}



template <class InputRange,
          class OutputIterator,
          class VertexPointMap,
          class Traits>
OutputIterator
extreme_vertices(const InputRange& range, 
                 OutputIterator out,
                 VertexPointMap vpm,
                 const Traits& traits)
{
  CGAL::Vertex_to_point_traits_adapter<Traits, VertexPointMap> traits_adapter(vpm, traits);
  CGAL::extreme_points_3(range, out,traits_adapter);
  return out;
}

void test_extreme_vertices(const char* fname)
{
  std::ifstream input(fname);
  Polyhedron_3 P;
  if (!input || !(input >> P) || P.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }  
  std::vector<boost::graph_traits<Polyhedron_3>::vertex_descriptor> verts;
  extreme_vertices(vertices(P), std::back_inserter(verts) , get(CGAL::vertex_point, P),
                   CGAL::Convex_hull_traits_3<K, Polyhedron_3, CGAL::Tag_true>());
  std::cout<<verts.size()<<std::endl;
}

int main()
{/*
  test_function_overload();
  test_3_points();
  test_up_to_3_extreme_points();
  test_3_collinear();
  test_triangulated_cube("data/cube_meshed.off");
  test_coplanar_points("data/coplanar_points.xyz");
  test_equal_points();*/
  
  test_extreme_vertices("data/cross.off");

  return 0;
}
