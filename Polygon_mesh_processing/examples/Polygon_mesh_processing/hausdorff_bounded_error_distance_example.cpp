#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/random_perturbation.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Real_timer.h>

#if defined(CGAL_LINKED_WITH_TBB)
#define TAG CGAL::Parallel_tag
#else
#define TAG CGAL::Sequential_tag
#endif

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Triangle_3 Triangle_3;
typedef K::Vector_3 Vector_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef Mesh::Vertex_index Vertex_index;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  // Objects to hold the meshes
  Mesh tm1, tm2;

  // Timer to measure the runtime of h-distance Distance_computation
  CGAL::Real_timer time;

  // Set an error bound
  double error_bound = 0.0001;

// ----------------------------------------------------------------------------

  // Easy Example of a tetrahedron and a remeshed version of itself

  // Create the Tetrahedron
  CGAL::make_tetrahedron(Point_3(.0,.0,.0),
                         Point_3(2,.0,.0),
                         Point_3(1,1,1),
                         Point_3(1,.0,2),
                         tm1);
  // Copy it and remesh it
  tm2=tm1;
  PMP::isotropic_remeshing(tm2.faces(),.05, tm2);
  // Compute the Hausdorff distance
  time.reset();
  time.start();
  std::cout << "Approximated Hausdorff distance: "
            << PMP::bounded_error_Hausdorff_distance
                  <TAG>(tm1, tm2, error_bound)
            << std::endl;
  time.stop();
  std::cout << "Processing took " << time.time() << "s." << std::endl;

// ----------------------------------------------------------------------------

  // Example with point realizing the Hausdorff distance strictly lying in the

  // interior of a triangle
  tm1 = Mesh();
  tm1.add_vertex( Point_3(-1.,1.,1.) );
  tm1.add_vertex( Point_3(0.,-1.,1.) );
  tm1.add_vertex( Point_3(1.,1.,1.) );
  tm1.add_face( tm1.vertices() );

  tm2 = Mesh();
  Vertex_index w0 = tm2.add_vertex( Point_3(-1.,1.,0.) );
  Vertex_index w1 = tm2.add_vertex( Point_3(0.,-1.,0.) );
  Vertex_index w2 = tm2.add_vertex( Point_3(1.,1.,0.) );
  Vertex_index w3 = tm2.add_vertex( Point_3(0.,1.,-1.) );
  Vertex_index w4 = tm2.add_vertex( Point_3(-0.5,0.,-1.) );
  Vertex_index w5 = tm2.add_vertex( Point_3(0.5,0.,-1.) );
  tm2.add_face( w0, w3, w4 );
  tm2.add_face( w1, w4, w5 );
  tm2.add_face( w2, w5, w3 );

  // Compute the Hausdorff distance
  time.reset();
  time.start();
  std::cout << "Approximated Hausdorff distance: "
            << PMP::bounded_error_Hausdorff_distance
                  <TAG>(tm1, tm2, error_bound)
            << std::endl;
  time.stop();
  std::cout << "Processing took " << time.time() << "s." << std::endl;

// ----------------------------------------------------------------------------

  // Read a real meshes given by the user, perturb it slightly and compute the
  // Hausdorff distance between the original mesh and its pertubation

  std::ifstream input( argv[1] );
  input >> tm1;
  std::cout << "Read a mesh with " << tm1.number_of_faces() << " triangles." << std::endl;

  // Copy the mesh and perturb it slightly
  tm2 = tm1;
  bool do_project = false;
  PMP::random_perturbation( tm2.vertices(), tm2, 0.1, do_project );
  std::cout << "Perturbed the input mesh, now computing the Hausdorff distance." << std::endl;

  // Compute the Hausdorff distance
  time.reset();
  time.start();
  std::cout << "Approximated Hausdorff distance: "
            << PMP::bounded_error_Hausdorff_distance
                  <TAG>(tm1, tm2, error_bound)
            << std::endl;
  time.stop();
  std::cout << "Processing took " << time.time() << "s." << std::endl;

// ----------------------------------------------------------------------------

  // Read two meshes given by the user, initially place them at their originally
  // given position. Move the second mesh in 300 steps away from the first one.
  // Print how the Hausdorff distance changes.

  std::ifstream input1( argv[1] );
  input1 >> tm1;
  std::cout << "Read a mesh with " << tm1.number_of_faces() << " triangles." << std::endl;

  std::ifstream input2( argv[2] );
  input2 >> tm2;
  std::cout << "Read a mesh with " << tm2.number_of_faces() << " triangles." << std::endl;

  CGAL::Bbox_3 bb = PMP::bbox(tm2);
  double dist = CGAL::approximate_sqrt( Vector_3(bb.xmax() - bb.xmin(), bb.ymax() - bb.ymin(), bb.zmax() - bb.zmin()).squared_length() );

  for (int i=0; i<300; i++) {
    PMP::transform(
       CGAL::Aff_transformation_3<K> ( CGAL::Translation(), Vector_3( 0.01*dist, 0.01*dist, 0.01*dist ) ),
       tm2
     );
    time.reset();
    time.start();
    std::cout << "Position: " << i << std::endl;
    std::cout << "Approximated Hausdorff distance: "
              << PMP::bounded_error_Hausdorff_distance
                    <TAG>(tm1, tm2, error_bound)
              << std::endl;
    time.stop();
    std::cout << "Processing took " << time.time() << "s." << std::endl;
  }
}
