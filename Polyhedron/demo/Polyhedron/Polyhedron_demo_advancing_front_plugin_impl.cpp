#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"
#include <Scene_polyhedron_item.h>

#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include "Scene_points_with_normal_item.h"

#include <CGAL/Advancing_front_surface_reconstruction.h>

Polyhedron*
advancing_front_reconstruct(const Point_set& points,
                            double sm_perimeter,
                            double sm_area)
{
  typedef CGAL::Advancing_front_surface_reconstruction<Kernel> Reconstruction;
  typedef Reconstruction::Triangulation_3 Triangulation_3;
  Polyhedron* output_mesh = new Polyhedron;

  Triangulation_3 dt(points.begin(), points.end());
  
  Reconstruction reconstruction(dt);

  reconstruction();

  CGAL::AFSR::construct_polyhedron(*output_mesh, reconstruction);

  return output_mesh;
}

