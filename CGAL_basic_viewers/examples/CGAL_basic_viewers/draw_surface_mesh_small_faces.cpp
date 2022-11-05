#include <CGAL/license/Surface_mesh.h>
#include <CGAL/Graphic_buffer.h>
#include <CGAL/Drawing_functor.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

struct Drawing_functor_small_face: public CGAL::Drawing_functor
<Mesh, Mesh::vertex_descriptor, Mesh::edge_descriptor, Mesh::face_descriptor>
{
   drawing_functor.colored_face = [](const FG&,
             typename boost::graph_traits<FG>::face_descriptor fg) -> bool
   { return is_small(fg); };

  drawing_functor.face_color = [] (const FG&,
             typename boost::graph_traits<FG>::face_descriptor) -> CGAL::IO::Color
  { return CGAL::IO::Color(200, 60, 60); }; // Red

  unsigned int m_threshold;
  FT m_min_size, m_max_size;
  bool m_draw_small_faces;
  bool m_draw_big_faces;
};

int main(int argc, char* argv[])
{
  const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");

  Mesh sm;
  if(!CGAL::IO::read_polygon_mesh(filename, sm))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  Drawing_functor_small_face drawing_functor;
  CGAL::draw(sm, drawing_functor);

  // setKeyDescription(Qt::Key_I, "Increment threshold for small faces");
  // setKeyDescription(Qt::Key_D, "Decrement threshold for small faces");
  // setKeyDescription(Qt::Key_S, "Draw small faces only , big faces only, both");

  return EXIT_SUCCESS;
}
