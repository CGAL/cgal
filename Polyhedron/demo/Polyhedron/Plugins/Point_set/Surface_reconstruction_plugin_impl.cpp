//----------------------------------------------------------
// Poisson reconstruction method:
// Reconstructs a surface mesh from a point set and returns it as a polyhedron.
//----------------------------------------------------------


// CGAL
#include <CGAL/AABB_tree.h> // must be included before kernel
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Timer.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/compute_average_spacing.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_solver_traits.h>
#endif

#include <math.h>

#include "Kernel_type.h"
#include "SMesh_type.h"
#include "Scene_points_with_normal_item.h"



// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif


// Poisson reconstruction method:
// Reconstructs a surface mesh from a point set and returns it as a polyhedron.
template<class FaceGraph, typename Traits>
bool poisson_reconstruct(FaceGraph* graph,
                         Point_set& points,
                         typename Traits::FT sm_angle, // Min triangle angle (degrees).
                         typename Traits::FT sm_radius, // Max triangle size w.r.t. point set average spacing.
                         typename Traits::FT sm_distance, // Approximation error w.r.t. point set average spacing.
                         const QString& solver_name, // solver name
                         bool use_two_passes,
                         bool do_not_fill_holes)
{
  // Poisson implicit function
  typedef CGAL::Poisson_reconstruction_function<Traits> Poisson_reconstruction_function;

  // Surface mesher
  typedef CGAL::Surface_mesh_default_triangulation_3 STr;
  typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
  typedef CGAL::Implicit_surface_3<Traits, Poisson_reconstruction_function> Surface_3;

  // AABB tree
  typedef CGAL::AABB_face_graph_triangle_primitive<FaceGraph> Primitive;
  typedef CGAL::AABB_traits<Traits, Primitive> AABB_traits;
  typedef CGAL::AABB_tree<AABB_traits> AABB_tree;
    CGAL::Timer task_timer; task_timer.start();

    //***************************************
    // Checks requirements
    //***************************************

    if (points.size() == 0)
    {
      std::cerr << "Error: empty point set" << std::endl;
     return false;
    }

    bool points_have_normals = points.has_normal_map();
    if ( ! points_have_normals )
    {
      std::cerr << "Input point set not supported: this reconstruction method requires oriented normals" << std::endl;
     return false;
    }

    CGAL::Timer reconstruction_timer; reconstruction_timer.start();

    //***************************************
    // Computes implicit function
    //***************************************

 
    std::cerr << "Computes Poisson implicit function "
              << "using " << solver_name.toLatin1().data() << " solver...\n";
              
    
    // Creates implicit function from the point set.
    // Note: this method requires an iterator over points
    // + property maps to access each point's position and normal.
    Poisson_reconstruction_function function(points.begin_or_selection_begin(), points.end(),
                                             points.point_map(), points.normal_map());

    bool ok = false;    
    #ifdef CGAL_EIGEN3_ENABLED
    if(solver_name=="Eigen - built-in simplicial LDLt")
    {
      CGAL::Eigen_solver_traits<Eigen::SimplicialCholesky<CGAL::Eigen_sparse_matrix<double>::EigenType> > solver;
      ok = function.compute_implicit_function(solver, use_two_passes);
    }
    if(solver_name=="Eigen - built-in CG")
    {
      CGAL::Eigen_solver_traits<Eigen::ConjugateGradient<CGAL::Eigen_sparse_matrix<double>::EigenType> > solver;
      solver.solver().setTolerance(1e-6);
      solver.solver().setMaxIterations(1000);
      ok = function.compute_implicit_function(solver, use_two_passes);
    }
    #endif

    // Computes the Poisson indicator function f()
    // at each vertex of the triangulation.
    if ( ! ok )
    {
      std::cerr << "Error: cannot compute implicit function" << std::endl;
     return false;
    }

    // Prints status
    std::cerr << "Total implicit function (triangulation+refinement+solver): " << task_timer.time() << " seconds\n";
    task_timer.reset();

    //***************************************
    // Surface mesh generation
    //***************************************

    std::cerr << "Surface meshing...\n";

    // Computes average spacing
    Kernel::FT average_spacing = CGAL::compute_average_spacing<Concurrency_tag>(points.all_or_selection_if_not_empty(),
                                                                                6 /* knn = 1 ring */,
                                                                                points.parameters());

    // Gets one point inside the implicit surface
    Kernel::Point_3 inner_point = function.get_inner_point();
    Kernel::FT inner_point_value = function(inner_point);
    if(inner_point_value >= 0.0)
    {
      std::cerr << "Error: unable to seed (" << inner_point_value << " at inner_point)" << std::endl;
     return false;
    }

    // Gets implicit function's radius
    Kernel::Sphere_3 bsphere = function.bounding_sphere();
    Kernel::FT radius = std::sqrt(bsphere.squared_radius());

    // Defines the implicit surface: requires defining a
  	// conservative bounding sphere centered at inner point.
    Kernel::FT sm_sphere_radius = 5.0 * radius;
    Kernel::FT sm_dichotomy_error = sm_distance*average_spacing/1000.0; // Dichotomy error must be << sm_distance
    Surface_3 surface(function,
                      Kernel::Sphere_3(inner_point,sm_sphere_radius*sm_sphere_radius),
                      sm_dichotomy_error/sm_sphere_radius);

    // Defines surface mesh generation criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
                                                        sm_radius*average_spacing,  // Max triangle size
                                                        sm_distance*average_spacing); // Approximation error

    CGAL_TRACE_STREAM << "  make_surface_mesh(sphere center=("<<inner_point << "),\n"
                      << "                    sphere radius="<<sm_sphere_radius<<",\n"
                      << "                    angle="<<sm_angle << " degrees,\n"
                      << "                    triangle size="<<sm_radius<<" * average spacing="<<sm_radius*average_spacing<<",\n"
                      << "                    distance="<<sm_distance<<" * average spacing="<<sm_distance*average_spacing<<",\n"
                      << "                    dichotomy error=distance/"<<sm_distance*average_spacing/sm_dichotomy_error<<",\n"
                      << "                    Manifold_with_boundary_tag)\n";

    // Generates surface mesh with manifold option
    STr tr; // 3D Delaunay triangulation for surface mesh generation
    C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
    CGAL::make_surface_mesh(c2t3,                                 // reconstructed mesh
                            surface,                              // implicit surface
                            criteria,                             // meshing criteria
                            CGAL::Manifold_with_boundary_tag());  // require manifold mesh

    // Prints status
    std::cerr << "Surface meshing: " << task_timer.time() << " seconds, "
                                     << tr.number_of_vertices() << " output vertices"
                                     << std::endl;
    task_timer.reset();

    if(tr.number_of_vertices() == 0)
     return false;

    // Converts to polyhedron
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, *graph);

    // Prints total reconstruction duration
    std::cerr << "Total reconstruction (implicit function + meshing): " << reconstruction_timer.time() << " seconds\n";

    //***************************************
    // Computes reconstruction error
    //***************************************

    // Constructs AABB tree and computes internal KD-tree
    // data structure to accelerate distance queries
    AABB_tree tree(faces(*graph).first, faces(*graph).second, *graph);
    tree.accelerate_distance_queries();

    // Computes distance from each input point to reconstructed mesh
    double max_distance = DBL_MIN;
    double avg_distance = 0;

    std::set<typename boost::graph_traits<FaceGraph>::face_descriptor> faces_to_keep;
    
    for (Point_set::const_iterator p=points.begin_or_selection_begin(); p!=points.end(); p++)
    {
      typename AABB_traits::Point_and_primitive_id pap = tree.closest_point_and_primitive (points.point (*p));
      double distance = std::sqrt(CGAL::squared_distance (pap.first, points.point(*p)));
      
      max_distance = (std::max)(max_distance, distance);
      avg_distance += distance;

      typename boost::graph_traits<FaceGraph>::face_descriptor f = pap.second;
      faces_to_keep.insert (f);  
    }
    avg_distance /= double(points.size());

    std::cerr << "Reconstruction error:\n"
              << "  max = " << max_distance << " = " << max_distance/average_spacing << " * average spacing\n"
              << "  avg = " << avg_distance << " = " << avg_distance/average_spacing << " * average spacing\n";

    if (do_not_fill_holes)
      {
        typename boost::graph_traits<FaceGraph>::face_iterator it = faces(*graph).begin ();
        while (it != faces(*graph).end ())
	  {
            typename boost::graph_traits<FaceGraph>::face_iterator current = it ++;

            if (faces_to_keep.find (*current) == faces_to_keep.end ())
            {
              CGAL::Euler::remove_face(halfedge (*current, *graph), *graph);
            }

	  }

      }
    return true;
}

SMesh* poisson_reconstruct_sm(Point_set& points,
                              Kernel::FT sm_angle, // Min triangle angle (degrees).
                              Kernel::FT sm_radius, // Max triangle size w.r.t. point set average spacing.
                              Kernel::FT sm_distance, // Approximation error w.r.t. point set average spacing.
                              const QString& solver_name, // solver name
                              bool use_two_passes,
                              bool do_not_fill_holes)
{
  SMesh* res = new SMesh;
  if(poisson_reconstruct<SMesh, Kernel>(res, points, sm_angle, sm_radius, sm_distance,
                         solver_name, use_two_passes, do_not_fill_holes))
    return res;
  else
  {
    delete res;
    return NULL;
  }
}
