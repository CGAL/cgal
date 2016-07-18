#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>



#include <CGAL/Polygon_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/boost/graph/helpers.h>
#include <iostream>
typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K>                        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>              Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>        CDT;
typedef CGAL::Polygon_2<K>                                        Polygon_2;
using namespace CGAL;
int test_triangles_3()
{
#include <iostream>
#include <fstream>
using namespace CGAL;
typedef Simple_cartesian<double>                           K;
typedef K::Point_3                                         Point;


int main()
{
 // Generated points are in that vector
  std::vector<Point> points;
  // Create input triangles
  std::vector<K::Triangle_3> triangles;
  for(int i=0; i< 5; ++i)
  {
    triangles.push_back(K::Triangle_3(Point(i,0,0.5*i), Point(i+1,0,0.5*i), Point(i+0.5,1,0.5*i)));
  }

  // Create the generator, input is the vector of Triangle_3
  Random_points_on_triangles_3<Point> g(triangles);
  // Get 100 random points in cdt
  CGAL::cpp11::copy_n(g, 1000, std::back_inserter(points));

  // Check that we have really created 100 points.
  assert( points.size() == 1000);

  // print the first point that was generated
  std::cout << points[0] << std::endl;

  return 0;
}



}

int test_T2()
{
  typedef CDT::Point                                              Point_2;
  std::vector<Point_2> points;
//construct two non-intersecting nested polygons
::Polygon_2 polygon1;
polygon1.push_back(Point_2(0,0));
polygon1.push_back(Point_2(2,0));
polygon1.push_back(Point_2(2,2));
polygon1.push_back(Point_2(0,2));

//Insert the polygons into a constrained triangulation
CDT cdt;
cdt.insert_constraint(polygon1.vertices_begin(), polygon1.vertices_end(), true);

Random_points_in_triangle_mesh_2<Point_2, CDT>
    g(cdt);
cpp11::copy_n( g, 300, std::back_inserter(points));
for(std::size_t i = 0; i<points.size(); ++i)
{
  Point_2 p= points[i];
  for(int j = 0; j<2; ++j)
  {
    double coords[2] = {p.x(), p.y()};
    if(coords[j]>2.05 || coords[j]<-0.05)
    {
      std::cerr<<"ERROR : Generated point is not on the cube."<<std::endl;
      return 0;
    }
  }
  return 1;
}

return 1;
}


typedef CGAL::Polyhedron_3<K>                                     Polyhedron;
typedef K::Point_3                                                Point;
typedef K::FT                                                     FT;
bool on_face(int face, double coord[3])
{
    if(CGAL::abs(CGAL::abs(coord[face]) - 0.5) < 0.05
       && CGAL::abs(coord[(face+1)%3]) - 0.5 < 0.05
       && CGAL::abs(coord[(face+2)%3]) - 0.5 < 0.05)
      return true;
  return false;
}

int  test_volume_mesh(Polyhedron& polyhedron)
{


  std::vector<Point> points;
  Random_points_on_triangle_mesh_3<Point, Polyhedron>
      g(polyhedron);
  CGAL::cpp11::copy_n( g, 300, std::back_inserter(points));
  for (std::size_t i = 0; i<points.size(); ++i)
  {
    Point p= points[i];
    double coords[3] = {p.x(), p.y(), p.z()};
    if(!(on_face(0, coords) || on_face(1, coords) || on_face(2,coords)))
    {
      std::cerr<<"ERROR : Generated point is not on the cube."<<std::endl;
      return 0;
    }
  }
  return 1;
}


// Domain
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef C3t3::Point                                   Point_c3t3;

int test_on_c3t3(const Polyhedron& polyhedron)
{
  std::vector<Point_c3t3> points;
  points.clear();
  // Create domain
  Mesh_domain domain(polyhedron);
   using namespace CGAL::parameters;
  // Mesh criteria (no cell_size set)
  Mesh_criteria criteria(facet_angle=25, facet_size=0.15, facet_distance=0.008,
                         cell_radius_edge_ratio=3);
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());
  Random_points_on_tetrahedral_mesh_boundary<C3t3>
      g(c3t3);
  CGAL::cpp11::copy_n( g, 300, std::back_inserter(points));
  for (std::size_t i = 0; i<points.size(); ++i)
  {
    Point p= points[i];
    double coords[3] = {p.x(), p.y(), p.z()};
    if(!(on_face(0, coords) || on_face(1, coords) || on_face(2,coords)))
    {
      std::cerr<<"ERROR : Generated point is not on the cube."<<std::endl;
      return 0;
    }
  }
    return 1;
}

int test_in_c3t3(const Polyhedron& polyhedron)
{
  std::vector<Point> points;
  points.clear();
  // Create domain
  Mesh_domain domain(polyhedron);
   using namespace CGAL::parameters;
  // Mesh criteria (no cell_size set)
  Mesh_criteria criteria(facet_angle=25, facet_size=0.15, facet_distance=0.008,
                         cell_radius_edge_ratio=3);
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

  Random_points_in_tetrahedral_mesh_3<C3t3>
      g(c3t3);
  CGAL::cpp11::copy_n( g, 300, std::back_inserter(points));
  for (std::size_t i = 0; i<points.size(); ++i)
  {
    Point p= points[i];
    double coords[3] = {p.x(), p.y(), p.z()};
    for(int j = 0; j< 3; ++j)
      if(CGAL::abs(coords[j]) > 0.501)
      {
        std::cerr<<"ERROR : Generated point is not in the cube."<<std::endl;
        return 0;
      }
  }
    return 1;
}

int
main( )
{
  Polyhedron polyhedron;

  make_hexahedron(Point(-0.5,-0.5,-0.5), Point(0.5,-0.5,-0.5), Point(0.5,0.5,-0.5), Point(-0.5,0.5,-0.5),
                  Point(-0.5,0.5,0.5), Point(-0.5,-0.5,0.5), Point(0.5,-0.5,0.5), Point(0.5,0.5,0.5),
                  polyhedron);
  boost::graph_traits<Polyhedron>::halfedge_descriptor facets[6];
  int i = 0;
  BOOST_FOREACH(boost::graph_traits<Polyhedron>::face_descriptor fd, faces(polyhedron))
    facets[i++] = halfedge(fd, polyhedron);
  for(int i=0; i<6; ++i)
    CGAL::Euler::split_face(facets[i],next(next(facets[i], polyhedron), polyhedron), polyhedron);


  int validity = test_volume_mesh(polyhedron)
      *test_T2()
      *test_on_c3t3(polyhedron)
      *test_in_c3t3(polyhedron)
      ;
  assert(validity == 1);
  return 0;
}
