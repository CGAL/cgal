#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/draw_polyhedron.h>
#include <CGAL/draw_point_set_3.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Qt/Basic_viewer.h>

#include <vector>
#include <iostream>

using Kernel=CGAL::Exact_predicates_inexact_constructions_kernel;
using Point=Kernel::Point_3;
using Vector=Kernel::Vector_3;
using Pwn=std::pair<Point, Vector>;
using Polyhedron=CGAL::Polyhedron_3<Kernel>;
using PS3=CGAL::Point_set_3<Point>;

struct Graphics_scene_options_green_points:
  public CGAL::Graphics_scene_options<PS3, typename PS3::const_iterator,
                                      typename PS3::const_iterator,
                                      typename PS3::const_iterator>
{
  bool colored_vertex(const PS3&, typename PS3::const_iterator) const
  { return true; }
  CGAL::IO::Color vertex_color(const PS3&, typename PS3::const_iterator) const
  { return CGAL::IO::Color(0,220,0); }
};

int main(void)
{
  std::vector<Pwn> points;

  if(!CGAL::IO::read_points(CGAL::data_file_path("points_3/kitten.xyz"), std::back_inserter(points),
                            CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>())
                                             .normal_map(CGAL::Second_of_pair_property_map<Pwn>())))
  {
    std::cerr << "Error: cannot read input file " << CGAL::data_file_path("points_3/kitten.xyz") << std::endl;
    return EXIT_FAILURE;
  }

  Polyhedron output_mesh;

  double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
    (points, 6, CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>()));

  if (CGAL::poisson_surface_reconstruction_delaunay
      (points.begin(), points.end(),
       CGAL::First_of_pair_property_map<Pwn>(),
       CGAL::Second_of_pair_property_map<Pwn>(),
       output_mesh, average_spacing))
  {
    PS3 point_set;
    for(Pwn& it: points)
    { point_set.insert(it.first); }

    CGAL::Graphics_scene scene;
    CGAL::add_to_graphics_scene(point_set, scene, Graphics_scene_options_green_points());
    CGAL::add_to_graphics_scene(output_mesh, scene);
    CGAL::draw_graphics_scene(scene);
  }
  else
  { return EXIT_FAILURE; }

  return EXIT_SUCCESS;
}
