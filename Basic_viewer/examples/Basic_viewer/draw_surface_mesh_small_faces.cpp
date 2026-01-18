#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <iostream>
#include <string>
#include <cassert>
#include <CGAL/draw_surface_mesh.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef Mesh::Vertex_index             vertex_descriptor;
typedef Mesh::Face_index               face_descriptor;
typedef K::FT                          FT;

template<class SM>
struct Graphics_scene_options_small_faces:
  public CGAL::Graphics_scene_options<SM,
                               typename boost::graph_traits<SM>::vertex_descriptor,
                               typename boost::graph_traits<SM>::edge_descriptor,
                               typename boost::graph_traits<SM>::face_descriptor>
{
  using Base=CGAL::Graphics_scene_options<SM,
                                          typename boost::graph_traits<SM>::vertex_descriptor,
                                          typename boost::graph_traits<SM>::edge_descriptor,
                                          typename boost::graph_traits<SM>::face_descriptor>;

  Graphics_scene_options_small_faces(const SM& sm): Base(), m_sm(sm)
  {
    std::optional<typename SM::template Property_map<face_descriptor, FT>> faces_size
      = sm.template property_map<face_descriptor, FT>("f:size");
    m_with_size = faces_size.has_value();
    if(!m_with_size)
    { return; }

    m_min_size=faces_size.value()[*(sm.faces().begin())];
    m_max_size=m_min_size;
    FT cur_size;
    for (typename SM::Face_range::iterator f=sm.faces().begin(); f!=sm.faces().end(); ++f)
    {
      cur_size=faces_size.value()[*f];
      if (cur_size<m_min_size) m_min_size=cur_size;
      if (cur_size>m_max_size) m_max_size=cur_size;
    }

    this->face_color=[this](const SM& sm,
                         typename boost::graph_traits<SM>::face_descriptor fh) -> CGAL::IO::Color
    { return this->get_face_color(sm, fh); };

    this->colored_face = [](const SM&,
                            typename boost::graph_traits<SM>::face_descriptor) -> bool
    { return true; };
  }

  CGAL::IO::Color get_face_color(const SM& sm,
                                 typename boost::graph_traits<SM>::face_descriptor fh)
  {
    // Default color of faces
    CGAL::IO::Color c(75,160,255);
    if(!m_with_size) { return c; }

    // Compare the size of the face with the % m_threshold
    std::optional<typename SM::template Property_map<face_descriptor, FT>> faces_size
      = sm.template property_map<face_descriptor, FT>("f:size");
    assert(faces_size.has_value());

    // If the face is small, color it in red.
    if (get(faces_size.value(), fh)<m_min_size + ((m_max_size - m_min_size) / (100 - m_threshold)))
    { return CGAL::IO::Color(255,20,20); }

    return c; // Default color
  }

  const SM& m_sm;
  bool m_with_size;
  FT m_min_size, m_max_size;
  unsigned int m_threshold=85; // 85%
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

  CGAL::Polygon_mesh_processing::triangulate_faces(sm);

  Mesh::Property_map<face_descriptor, FT> faces_size;
  bool created;
  std::tie(faces_size, created)=sm.add_property_map<face_descriptor, FT>("f:size",0.);
  assert(created);

  for(face_descriptor fd : sm.faces())
  { faces_size[fd]=CGAL::Polygon_mesh_processing::face_area(fd, sm); }

  Graphics_scene_options_small_faces gsosm(sm);
  CGAL::Graphics_scene gs;

  add_to_graphics_scene(sm, gs, gsosm);

#ifdef CGAL_USE_BASIC_VIEWER

  CGAL::Qt::QApplication_and_basic_viewer app(gs, "Small faces");
  if(app)
  {
    app.basic_viewer().on_key_pressed=
      [&sm, &gsosm, &gs] (QKeyEvent* e, CGAL::Qt::Basic_viewer* basic_viewer) -> bool
      {
        const ::Qt::KeyboardModifiers modifiers = e->modifiers();
        if ((e->key() == ::Qt::Key_I) && (modifiers == ::Qt::NoButton))
        {
          gsosm.m_threshold+=5;
          if(gsosm.m_threshold>100) { gsosm.m_threshold=100; }
          basic_viewer->displayMessage
            (QString("Small faces threshold=%1.").arg(gsosm.m_threshold));

          gs.clear();
          add_to_graphics_scene(sm, gs, gsosm);
          basic_viewer->redraw();
        }
        else if ((e->key() == ::Qt::Key_D) && (modifiers == ::Qt::NoButton))
        {
          if(gsosm.m_threshold<5) { gsosm.m_threshold=0; }
          else  { gsosm.m_threshold-=5; }
          basic_viewer->displayMessage
            (QString("Small faces threshold=%1.").arg(gsosm.m_threshold));

          gs.clear();
          add_to_graphics_scene(sm, gs, gsosm);
          basic_viewer->redraw();
        }
        else
        {
          // Return false will call the base method to process others/classicals key
          return false;
        }
        return true;
      };

    // Here we add shortcut descriptions
    app.basic_viewer().setKeyDescription(::Qt::Key_I, "Increase threshold for small faces");
    app.basic_viewer().setKeyDescription(::Qt::Key_D, "Decrease threshold for small faces");

    // Then we run the app
    app.run();
  }

#endif

  sm.remove_property_map(faces_size);

  return EXIT_SUCCESS;
}
