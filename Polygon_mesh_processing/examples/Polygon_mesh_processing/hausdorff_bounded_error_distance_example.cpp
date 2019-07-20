#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/random_perturbation.h>
#include <CGAL/Real_timer.h>

#if defined(CGAL_LINKED_WITH_TBB)
#define TAG CGAL::Parallel_tag
#else
#define TAG CGAL::Sequential_tag
#endif

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Triangle_3 Triangle_3;
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
  double error_bound = 0.01;

/*
  // Easy Example
  CGAL::make_tetrahedron(Point_3(.0,.0,.0),
                         Point_3(2,.0,.0),
                         Point_3(1,1,1),
                         Point_3(1,.0,2),
                         tm1);
  tm2=tm1;
  CGAL::Polygon_mesh_processing::isotropic_remeshing(tm2.faces(),.05, tm2);
*/
// ----------------------------------------------------------------------------
/*
  // Example with point realizing the Hausdorff distance strictly lying in the
  // interior of a triangle
  tm1 = Mesh();
  tm1.add_vertex( Point_3(-1.,1.,1.) );
  tm1.add_vertex( Point_3(0.,-1.,1.) );
  tm1.add_vertex( Point_3(1.,1.,1.) );
  tm1.add_face( tm1.vertices() );

  std::cout << "TM1 is valid: " << (tm1.is_valid() ? "true" : "false") << std::endl;

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

  std::cout << "TM2 is valid: " << (tm2.is_valid() ? "true" : "false") << std::endl;
*/
// ----------------------------------------------------------------------------
/*
  // Read a real mesh given by the user
  std::ifstream input( argv[1] );
  input >> tm1;
  std::cout << "Read a mesh with " << tm1.number_of_faces() << " triangles." << std::endl;

  std::ifstream input2( argv[2] );
  input2 >> tm2;
  std::cout << "Read a mesh with " << tm2.number_of_faces() << " triangles." << std::endl;

  // Copy the mesh and perturb it slightly
  tm2 = tm1;
  bool do_project = false;
  CGAL::Polygon_mesh_processing::random_perturbation( tm2.vertices(), tm2, 0.1, do_project );
  std::cout << "Perturbed the input mesh, now computing the Hausdorff distance." << std::endl;
*/
// ----------------------------------------------------------------------------

  // Pairwise computation on a set of benchmark models
  std::vector<std::string> models = {"80","102","128","162","204","256","320","402","504","632","792","992","1242"};
  int num_models = models.size();
  std::string prefix = "/home/martin/Downloads/bunnies/bunny_";
  std::string postfix = ".off";
  std::vector<double> error_bounds = { 0.1, 0.01, 0.001, 0.0001 };

  std::ifstream input(prefix + models[0] + postfix);
  // input >> tm1;
  // input >> tm2;
  // std::cout << "Initialized with a mesh at " << (prefix + models[0] + postfix) << " with " << tm1.number_of_faces() << " triangles." << std::endl;

  for(int i=0; i<num_models; i++) {

    input.close();
    input.open(prefix + models[i] + postfix);
    tm1 = Mesh();
    input >> tm1;
    // std::cout << "Read a mesh at " << (prefix + models[i] + postfix) << " with " << tm1.number_of_faces() << " faces." << std::endl;

    for(int j=0; j<num_models; j++) {
      if (i == j) continue;
      input.close();
      input.open(prefix + models[j] + postfix);
      tm2 = Mesh();
      input >> tm2;
      // std::cout << "Read a mesh at " << (prefix + models[j] + postfix) << " with " << tm2.number_of_faces() << " faces." << std::endl;
      //
      // std::cout << "Read two meshes with " << tm1.number_of_faces() << ", " << tm2.number_of_faces() << " triangles respectively." << std::endl;

      for (int k=0; k<error_bounds.size(); k++) {
        time.reset();
        time.start();
        double h_dist = CGAL::Polygon_mesh_processing::bounded_error_Hausdorff_distance<TAG>(tm1, tm2, error_bounds[k]);
        time.stop();
        std::cout << models[i] << " " << models[j] << " " << time.time() << " " << error_bounds[k] << " " << h_dist << std::endl;
      }
    }
  }


/*
  time.start();
  std::cout << "Approximated Hausdorff distance: "
            << CGAL::Polygon_mesh_processing::bounded_error_Hausdorff_distance
                  <TAG>(tm1, tm2, error_bound)
            << std::endl;
  time.stop();
  std::cout << "Processing took " << time.time() << "s." << std::endl;
*/
/*
  std::cout << "Approximated Hausdorff distance (naive): "
            << CGAL::Polygon_mesh_processing::bounded_error_Hausdorff_distance_naive
                  <TAG>(tm1, tm2, error_bound)
            << ", the actual distance is at most " << error_bound << " larger."
            << std::endl;
*/
}
