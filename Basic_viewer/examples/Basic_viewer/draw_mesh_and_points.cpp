#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/draw_polyhedron.h>
#include <CGAL/draw_point_set_3.h>
#include <CGAL/Graphics_scene_options.h>

#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point, Vector> Pwn;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
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
    std::cerr << "Error: cannot read input file!" << std::endl;
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

    CGAL::Graphics_scene<float> graphic_buffer;
    CGAL::add_in_graphic_storage(point_set, graphic_buffer, Graphics_scene_options_green_points());
    CGAL::add_in_graphic_storage(output_mesh, graphic_buffer);
    CGAL::draw_graphic_storage(graphic_buffer);
  }
  else
  { return EXIT_FAILURE; }

  return EXIT_SUCCESS;
}
