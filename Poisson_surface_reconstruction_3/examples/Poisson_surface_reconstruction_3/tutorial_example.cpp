#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Point_set_3/Point_set_processing_3.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Poisson_implicit_surface_3.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/Scale_space_surface_reconstruction_3.h>

#include <cstdlib>
#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Sphere_3 Sphere;
typedef CGAL::Point_set_3<Point, Vector> Point_set;

int main(int argc, char*argv[])
{
  ///////////////////////////////////////////////////////////////////
  //! [Reading input]
  
  Point_set points;

  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " [input.xyz/off/ply] (output.off)" << std::endl;
    return -1;
  }

  const char* input_file = argv[1];
  std::ifstream stream (input_file);
  if (!stream)
  {
    std::cerr << "Error: cannot read file " << input_file << std::endl;
    return -1;
  }

  stream >> points;
  
  std::cout << "Read " << points.size () << " point(s)" << std::endl;
  if (points.empty())
    return -1;

  //! [Reading input]
  ///////////////////////////////////////////////////////////////////
  
  ///////////////////////////////////////////////////////////////////
  //! [Outlier removal]

  CGAL::remove_outliers (points,
                         24, // Number of neighbors considered for evaluation
                         5.0); // Percentage of points to remove

  // Applying point set processing algorithm to a CGAL::Point_set_3
  // object does not erase the points from memory but place them in
  // the garbage of the object: memory can be freeed by the user.
  std::cout << points.number_of_removed_points()
	    << " point(s) are inliers." << std::endl;

  points.collect_garbage();

  //! [Outlier removal]
  ///////////////////////////////////////////////////////////////////
  
  ///////////////////////////////////////////////////////////////////
  //! [Simplification]
  
  // Compute average spacing using neighborhood of 6 points
  double spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag> (points, 6);

  // Simplify using a grid of size 2 * average spacing
  CGAL::grid_simplify_point_set (points, 2. * spacing);

  std::cout << points.number_of_removed_points()
	    << " point(s) remaining after simplification." << std::endl;

  points.collect_garbage();

  //! [Simplification]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Smoothing]

  CGAL::jet_smooth_point_set<CGAL::Sequential_tag> (points, 24);

  //! [Smoothing]
  ///////////////////////////////////////////////////////////////////

  unsigned int reconstruction_choice
    = (argc < 3 ? 0 : atoi(argv[2]));

  if (reconstruction_choice == 0) // Poisson
  {
    ///////////////////////////////////////////////////////////////////
    //! [Normal estimation]

    CGAL::jet_estimate_normals<CGAL::Sequential_tag>
      (points, 24); // Use 24 neighbors

    // Orientation of normals, returns iterator to first unoriented point
    typename Point_set::iterator unoriented_points_begin =
      CGAL::mst_orient_normals(points, 24); // Use 24 neighbors

    points.remove (unoriented_points_begin, points.end());

    //! [Normal estimation]
    ///////////////////////////////////////////////////////////////////
      
    ///////////////////////////////////////////////////////////////////
    //! [Poisson reconstruction]
      
    typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;      
    Poisson_reconstruction_function function
      (points.begin(), points.end(), points.point_map(), points.normal_map());

    if ( ! function.compute_implicit_function() )
    {
      std::cerr << "Error: cannot compute implicit function.";
      return EXIT_FAILURE;
    }

    //! [Poisson reconstruction]
    ///////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////
    //! [Surface mesh generation]
      
    typedef CGAL::Surface_mesh_default_triangulation_3 STr;
    typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
    typedef CGAL::Poisson_implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;

    // Gets one point inside the implicit surface
    // and computes implicit function bounding sphere radius.
    Point inner_point = function.get_inner_point();
    Sphere bsphere = function.bounding_sphere();
    FT radius = std::sqrt(bsphere.squared_radius());

    FT sm_angle = 20.0; // Min triangle angle (degrees).
    FT sm_radius = 100; // Max triangle size w.r.t. point set average spacing.
    FT sm_distance = 0.25; // Approximation error w.r.t. point set average spacing.
      
    // Defines the implicit surface: requires defining a
    // conservative bounding sphere centered at inner point.
    FT sm_sphere_radius = 1.5 * radius;
    FT sm_dichotomy_error = sm_distance * spacing/1000.0;
    // Dichotomy error must be << sm_distance
      
    Surface_3 surface(function,
                      Sphere(inner_point,sm_sphere_radius*sm_sphere_radius),
                      sm_dichotomy_error/sm_sphere_radius);

    // Defines surface mesh generation criteria
    CGAL::Surface_mesh_default_criteria_3<STr>
      criteria(sm_angle,  // Min triangle angle (degrees)
               sm_radius * spacing,  // Max triangle size
               sm_distance * spacing); // Approximation error

    // Generates surface mesh with manifold option
    STr tr; // 3D Delaunay triangulation for surface mesh generation
    C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
    CGAL::make_surface_mesh(c2t3,                                 // reconstructed mesh
                            surface,                              // implicit surface
                            criteria,                             // meshing criteria
                            CGAL::Manifold_with_boundary_tag());  // require manifold mesh

    //! [Surface mesh generation]
    ///////////////////////////////////////////////////////////////////
      
    ///////////////////////////////////////////////////////////////////      
    //! [Output poisson]

    CGAL::Polyhedron_3<Kernel> output_mesh;
    // Convert to Polyhedron
    CGAL::output_surface_facets_to_polyhedron(c2t3, output_mesh);
      
    std::ofstream f ("out.off");
    f << output_mesh;
    f.close ();

    //! [Output poisson]
    ///////////////////////////////////////////////////////////////////      
  }
  else if (reconstruction_choice == 1) // Advancing front
  {
    ///////////////////////////////////////////////////////////////////      
    //! [Advancing front reconstruction]
      
    typedef CGAL::cpp11::array<std::size_t,3> Facet; // Triple of indices

    std::vector<Facet> facets;

    // The function is called using directly the points raw iterators
    CGAL::advancing_front_surface_reconstruction(points.points().begin(),
                                                 points.points().end(),
                                                 std::back_inserter(facets));
    std::cout << facets.size ()
              << " facet(s) generated by reconstruction." << std::endl;

    //! [Advancing front reconstruction]
    ///////////////////////////////////////////////////////////////////      
      
    ///////////////////////////////////////////////////////////////////      
    //! [Output advancing front]
      
    std::ofstream f ("out.off");
  
    f << "OFF" << std::endl << points.size () << " " << facets.size () << " 0" << std::endl;

    for (typename Point_set::iterator it = points.begin(); it != points.end(); ++ it)
      f << points.point(*it) << std::endl;

    for (std::size_t i = 0; i < facets.size (); ++ i)
    {
      f << "3";
      for (std::size_t j = 0; j < 3; ++ j)
        f << " " << facets[i][j];
      f << std::endl;
    }

    f.close ();

    //! [Output advancing front]
    ///////////////////////////////////////////////////////////////////      
  }
  else if (reconstruction_choice == 2) // Scale space
  {
    ///////////////////////////////////////////////////////////////////
    //! [Scale space reconstruction]

    CGAL::Scale_space_surface_reconstruction_3<Kernel> reconstruct
      (10, // Number of neighborhood points
       300); // Number of samples used to estimate neighborhood radius

    reconstruct.reconstruct_surface(points.points().begin (), points.points().end (),
                                    4, // Number of iterations
                                    false, // Do not separate connected components
                                    true); // Force manifold output

    //! [Scale space reconstruction]
    ///////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////
    //! [Output scale space]
      
    std::ofstream f ("out.off");
    f << "OFF" << std::endl << points.size () << " "
      << reconstruct.number_of_triangles() << " 0" << std::endl;

    for (typename Point_set::iterator it = points.begin(); it != points.end(); ++ it)
      f << points.point(*it) << std::endl;

    typedef typename CGAL::Scale_space_surface_reconstruction_3<Kernel>::Triple_iterator Triple_iterator;
    for (Triple_iterator it = reconstruct.surface_begin();
         it != reconstruct.surface_end(); ++it)
      f << "3 "<< *it << std::endl;

    f.close ();

    //! [Output scale space]
    ///////////////////////////////////////////////////////////////////
  }
  else // Handle error
  {
    std::cerr << "Error: invalid reconstruction id: " << reconstruction_choice << std::endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}
