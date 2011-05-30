#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
 #include <CGAL/Delaunay_triangulation_3.h>
 #include <CGAL/Triangulation_vertex_base_with_info_3.h>
 #include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
 #include <CGAL/Mesh_triangulation_3.h>
 #include <CGAL/Mesh_complex_3_in_triangulation_3.h>
 #include <CGAL/Mesh_criteria_3.h>
 #include <CGAL/Polyhedral_mesh_domain_3.h>
 #include <CGAL/make_mesh_3.h>
 #include <CGAL/refine_mesh_3.h>
 #include <CGAL/IO/Polyhedron_iostream.h>
 #include <CGAL/Mesh_domain_with_polyline_features_3.h>
 
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
 typedef CGAL::Mesh_3::Robust_intersection_traits_3<K> Geom_traits;
 typedef Geom_traits::Point_3 Point_3;
 typedef CGAL::Polyhedron_3<Geom_traits> Polyhedron;
 typedef CGAL::Mesh_domain_with_polyline_features_3<
 CGAL::Polyhedral_mesh_domain_3<Polyhedron, Geom_traits> > MD_IntersectedLines;
 typedef CGAL::Mesh_triangulation_3<MD_IntersectedLines>::type Tr_MDI;
 typedef CGAL::Mesh_criteria_3<Tr_MDI> Mesh_criteria_MDI;
 typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr_MDI> C3t3_MDI;
 typedef C3t3_MDI::Facets_in_complex_iterator Facet_Iter_inComplex;
 typedef std::vector<Point_3>      InterPolyline;
 typedef std::list<InterPolyline>       InterPolylines;
using namespace CGAL::parameters;
 
int main()
 {
 
//constructing the bounding box by triangles
Polyhedron bounding_box;   
  //top face
  bounding_box.make_triangle(Point_3(100,100,100), Point_3(0,100,100), Point_3(0,0,100));
  bounding_box.make_triangle(Point_3(100,100,100), Point_3(0,0,100), Point_3(100,0,100));
  //bottom face
  bounding_box.make_triangle(Point_3(100,100,0), Point_3(100,0,0), Point_3(0,0,0));
  bounding_box.make_triangle(Point_3(100,100,0), Point_3(0,0,0), Point_3(0,100,0));
  //front face
  bounding_box.make_triangle(Point_3(0,0,0), Point_3(100,0,100), Point_3(0,0,100));
  bounding_box.make_triangle(Point_3(0,0,0), Point_3(100,0,0), Point_3(100,0,100));
  //back face
  bounding_box.make_triangle(Point_3(0,100,0), Point_3(0,100,100), Point_3(100,100,100));
  bounding_box.make_triangle(Point_3(100,100,100), Point_3(100,100,0), Point_3(0,100,0));
  //left face
  bounding_box.make_triangle(Point_3(0,0,0), Point_3(0,0,100), Point_3(0,100,100));
  bounding_box.make_triangle(Point_3(0,0,0), Point_3(0,100,100), Point_3(0,100,0));
  //right face
  bounding_box.make_triangle(Point_3(100,100,100), Point_3(100,0,0), Point_3(100,100,0));
  bounding_box.make_triangle(Point_3(100,100,100), Point_3(100,0,100), Point_3(100,0,0));
 
//constructing the input polyhedral face
 std::vector<Polyhedron*> AllSurfaces;  
 Polyhedron NewPoly_1;
  NewPoly_1.make_triangle(Point_3(20,20,20), Point_3(25,20,20), Point_3(25,40,20));
  NewPoly_1.make_triangle(Point_3(20,20,20), Point_3(25,40,20), Point_3(20,40,20));
  AllSurfaces.push_back(&NewPoly_1);
  Polyhedron NewPoly_2;
  NewPoly_2.make_triangle(Point_3(25,40,20), Point_3(25,20,20), Point_3(60,20,20));
  NewPoly_2.make_triangle(Point_3(60,40,20), Point_3(25,40,20), Point_3(60,20,20));
  AllSurfaces.push_back(&NewPoly_2);

  std::ofstream file;
  file.open("input_1.off");
  file << NewPoly_1;
  file.close();
  file.open("input_2.off");
  file << NewPoly_2;
  file.close();
  file.open("bbox.off");
  file << bounding_box;
  file.close();

//constructing the line feature
  InterPolylines interPolySet;
  InterPolyline interLine;
  interLine.push_back(Point_3(25,40,20));
  interLine.push_back(Point_3(25,20,20));
  interPolySet.push_back(interLine);
 
//generating the mesh
 MD_IntersectedLines domain(AllSurfaces.begin(), AllSurfaces.end(), bounding_box);
  domain.add_features(interPolySet.begin(), interPolySet.end());
  Mesh_criteria_MDI criteria(edge_size = 5, facet_angle=20,  cell_size=5, facet_topology = CGAL::FACET_VERTICES_ON_SURFACE);
  C3t3_MDI c3t3 = CGAL::make_mesh_3<C3t3_MDI>(domain, criteria, no_exude(), no_perturb());

 //output
  std::ofstream medit_file("out_1.mesh");
  c3t3.output_to_medit(medit_file);
  medit_file.close();
 }

#if 0
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>

// Domain 
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_segment_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main()
{
  // Create domain
  Mesh_domain domain("data/fandisk.off");
  
  // Get sharp features
  domain.detect_features();

  // Mesh criteria
  Mesh_criteria criteria(edge_size = 0.025,
                         facet_angle = 25, facet_size = 0.05, facet_distance = 0.005,
                         cell_radius_edge_ratio = 3, cell_size = 0.05);
  
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file);
}
#endif
