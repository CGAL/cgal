#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>



#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>

#include <CGAL/point_generators_3.h>

#include <iostream>
#include <fstream>

using namespace CGAL;
void test_volume_mesh()
{
  typedef Simple_cartesian<double>                           R;
  typedef R::Point_3                                         Point;
  typedef R::FT                                              FT;
  typedef Surface_mesh<Point>                                Surface_mesh;

  std::vector<Point> points;
  Surface_mesh sm;
  std::ifstream in("../../../Polyhedron/demo/Polyhedron/data/star.off");
  in >> sm;
  CGAL_assertion(in && !sm.is_empty());


  Random_points_on_triangle_mesh_3<Point, Surface_mesh>
      g(sm);
  CGAL::cpp11::copy_n( g, 3000, std::back_inserter(points));
  for (std::size_t i = 0; i<points.size(); ++i)
    std::cerr<<points[i]<<std::endl;
}


typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K>                        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>              Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>        CDT;
typedef CGAL::Polygon_2<K>                                        Polygon_2;



void test_T2()
{
  typedef CDT::Point                                                Point;
  std::vector<Point> points;
//construct two non-intersecting nested polygons
::Polygon_2 polygon1;
polygon1.push_back(Point(0,0));
polygon1.push_back(Point(2,0));
polygon1.push_back(Point(2,2));
polygon1.push_back(Point(0,2));
::Polygon_2 polygon2;
polygon2.push_back(Point(4.0,-2.0));
polygon2.push_back(Point(4.0,2.0));
polygon2.push_back(Point(6.0,0.0));

//Insert the polygons into a constrained triangulation
CDT cdt;
cdt.insert_constraint(polygon1.vertices_begin(), polygon1.vertices_end(), true);
cdt.insert_constraint(polygon2.vertices_begin(), polygon2.vertices_end(), true);

Random_points_in_triangle_mesh_2<Point, CDT>
    g(cdt);
CGAL::cpp11::copy_n( g, 300, std::back_inserter(points));
for (std::size_t i = 0; i<points.size(); ++i)
  std::cerr<<points[i]<<std::endl;
}

typedef CGAL::Polyhedron_3<K>                         Polyhedron;
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
typedef C3t3::Point                                   Point;

void test_on_c3t3()
{

  std::vector<Point> points;
  points.clear();
  // Create input polyhedron
  Polyhedron polyhedron;
  std::ifstream input("../../../Polyhedron/demo/Polyhedron/data/star.off");
  input >> polyhedron;
  input.close();
  // Create domain
  Mesh_domain domain(polyhedron);
   using namespace CGAL::parameters;
  // Mesh criteria (no cell_size set)
  Mesh_criteria criteria(facet_angle=25, facet_size=0.15, facet_distance=0.008,
                         cell_radius_edge_ratio=3);
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());
  Random_points_on_tetrahedral_mesh_3<C3t3>
      g(c3t3);
  CGAL::cpp11::copy_n( g, 300, std::back_inserter(points));
  for (std::size_t i = 0; i<points.size(); ++i)
    std::cerr<<points[i].x()<<" "<< points[i].y()<<" "<<points[i].z()<<std::endl;
}

void test_in_c3t3()
{


  std::vector<Point> points;
  points.clear();
  // Create input polyhedron
  Polyhedron polyhedron;
  std::ifstream input("../../../Polyhedron/demo/Polyhedron/data/star.off");
  input >> polyhedron;
  input.close();
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
    std::cerr<<points[i].x()<<" "<< points[i].y()<<" "<<points[i].z()<<std::endl;
}

int
main( )
{
  test_volume_mesh();
  test_T2();
  test_on_c3t3();
  test_in_c3t3();

  return 0;
}
