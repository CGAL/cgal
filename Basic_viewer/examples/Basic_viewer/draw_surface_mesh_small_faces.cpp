#include <CGAL/license/Surface_mesh.h>
#include <CGAL/Graphic_buffer.h>
#include <CGAL/Drawing_functor.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>
#include "CGAL/draw_surface_mesh_small_faces.h"

#include <fstream>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;
typedef Kernel::FT                                             FT;

struct Drawing_functor_small_face: public CGAL::Drawing_functor
<Mesh, Mesh::Vertex_index, Mesh::Edge_index, Mesh::Face_index>
{
 Drawing_functor_small_face() {

  // TODO: change threshold in realtime.
  m_threshold = 2000;

  // Return false if face value less than threshold.
  is_small = [=] (typename boost::graph_traits<Mesh>::face_descriptor fg) -> bool {
    return fg <= m_threshold;
  };

  colored_face = [=](const Mesh&,
             typename boost::graph_traits<Mesh>::face_descriptor fg) -> bool
   { return is_small(fg); };

  face_color = [] (const Mesh&,
             typename boost::graph_traits<Mesh>::face_descriptor) -> CGAL::IO::Color
  { return CGAL::IO::Color(200, 60, 60); }; // Red

  }

  std::function<bool(typename boost::graph_traits<Mesh>::face_descriptor)> is_small;

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
