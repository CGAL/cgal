#include <CGAL/AABB_polyhedral_oracle.h>
#include <CGAL/AABB_tree.h>

#include "Polyhedron_type.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_polyhedron_builder.h>

#include <CGAL/Timer.h>

#include <CGAL/assertions_behaviour.h>

#include <algorithm>

Polyhedron* cgal_code_remesh(const Polyhedron* pMesh,
                             const double angle,
                             const double sizing,
                             const double approx) {
  CGAL::set_error_behaviour(CGAL::ABORT);
  if(!pMesh) return 0;

  // remesh

  typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
  typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
  typedef Tr::Geom_traits GT;

  Tr triangulation; // 3D-Delaunay triangulation
  C2t3 c2t3(triangulation); // 2D-complex in 3D-Delaunay triangulation

  // meshing parameters
  CGAL::Surface_mesh_default_criteria_3<Tr> facets_criteria(angle,sizing,approx);

  // AABB tree
  CGAL::Timer timer;
  timer.start();
  std::cerr << "Build AABB tree...";
  typedef CGAL::Simple_cartesian<double> Simple_cartesian_kernel;
  // input surface
  typedef CGAL::AABB_polyhedral_oracle<Polyhedron,GT,Simple_cartesian_kernel> Input_surface;
  Input_surface input(*pMesh);
  std::cerr << "done (" << timer.time() << " ms)" << std::endl;

  // initial point set
  timer.reset();
  std::cerr << "Insert initial point set...";
  unsigned int nb_initial_points = 10;
  typedef CGAL::Cartesian_converter<Kernel,GT> Converter;
  Converter convert;

  { // new scope for the initialization, so that the vector
    // polyhedron_points is destroyed as soon as the initialization is
    // finished
    std::vector<Point> polyhedron_points;
    polyhedron_points.reserve(pMesh->size_of_vertices());
    std::copy(pMesh->points_begin(), pMesh->points_end(), 
              std::back_inserter(polyhedron_points));

    for(int n = 0;
        n < nb_initial_points || (n < 10 * nb_initial_points && 
                                  triangulation.dimension() < 3 );
        n = triangulation.number_of_vertices())
    {
      const int pos = CGAL::default_random.get_int(0, polyhedron_points.size());
      triangulation.insert(convert(polyhedron_points[pos]));
    }
  }
  if(triangulation.dimension() < 3)
    return 0;

  std::cerr << "done (" << timer.time() << " ms)" << std::endl;

  // remesh
  timer.reset();
  std::cerr << "Remesh...";
  CGAL::make_surface_mesh(c2t3, input, input, facets_criteria, CGAL::Manifold_with_boundary_tag());
  std::cerr << "done (" << timer.time() << " ms, " << triangulation.number_of_vertices() << " vertices)" << std::endl;

  if(triangulation.number_of_vertices() > 0)
  {
    // add remesh as new polyhedron
    Polyhedron *pRemesh = new Polyhedron;
    CGAL::Complex_2_in_triangulation_3_polyhedron_builder<C2t3, Polyhedron> builder(c2t3);
    pRemesh->delegate(builder);
    return pRemesh;
  }
  else
    return 0;
}
