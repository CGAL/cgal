#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>
#include <CGAL/Graphics_scene_options.h>
#include <fstream>

using Kernel=CGAL::Simple_cartesian<double>;
using Point=Kernel::Point_3;
using Mesh=CGAL::Surface_mesh<Point>;

// Inherit from CGAL::Graphics_scene_options to get all the default values.
struct My_graphics_scene_options:
  public CGAL::Graphics_scene_options<Mesh,
                                      typename boost::graph_traits<Mesh>::vertex_descriptor,
                                      typename boost::graph_traits<Mesh>::edge_descriptor,
                                      typename boost::graph_traits<Mesh>::face_descriptor>
{
  // All vertices are colored.
  bool colored_vertex(const Mesh&,
                      typename boost::graph_traits<Mesh>::vertex_descriptor) const
  { return true; }

  // Change the color of vertices randomly.
  CGAL::IO::Color vertex_color(const Mesh&,
                               typename boost::graph_traits<Mesh>::vertex_descriptor) const
  {
    static bool v_green=true;
    v_green=!v_green;
    if(v_green) // 1 vertex out of two green (randomly)
    { return CGAL::IO::Color(0,220,0); }
    else // the others are blue
    { return CGAL::IO::Color(0,0,220); }
  }
};


int main(int argc, char* argv[])
{
  const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");

  Mesh sm;
  if(!CGAL::IO::read_polygon_mesh(filename, sm))
  {
    std::cerr << "Invalid input file: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  // Draw the mesh using the new graphics scene option.
  CGAL::draw(sm, My_graphics_scene_options());
  return EXIT_SUCCESS;
}
