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

Polyhedron* cgal_code_remesh(const Polyhedron* pMesh,
                             const double angle,
                             const double sizing,
                             const double approx) {
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
  Polyhedron::Point_const_iterator it;
  typedef CGAL::Cartesian_converter<Kernel,GT> Converter;
  Converter convert;
  unsigned int i = 0;
  for(it = pMesh->points_begin();
      it != pMesh->points_end(), i < nb_initial_points;
      it++, i++)
    triangulation.insert(convert(*it));
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
