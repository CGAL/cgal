#include <CGAL/Qt/Basic_viewer_qt.h>
#include <QApplication>
#include <QColor>
#include <QVariantAnimation>
#include <QSpinBox>
#include <QSlider>
#include <QShortcut>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <algorithm>
#include <thread>
#include <cstdlib>
#include <CGAL/Surface_mesh.h>
#include <CGAL/facets_in_complex_3_to_triangle_mesh.h>
#include <boost/iterator/filter_iterator.hpp>
#include "ui_draw_c3t3_surface.h"

template <typename Mesh>
int draw_c3t3_surface_from_surface_mesh(const Mesh& mesh)
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (cgal_test_suite) return EXIT_SUCCESS;

  int argc=1;
  const char* argv[2]={"surface_mesh_viewer","\0"};
  QApplication app(argc,const_cast<char**>(argv));

  using Point =
    typename boost::property_map_value<Mesh, CGAL::vertex_point_t>::type;
  using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
  using Vector = typename Kernel::Vector_3;

  auto vnormals = get(CGAL::dynamic_vertex_property_t<Vector>(), mesh);
  {
    // temporary face property map needed by `compute_normals`
    auto fpm = get(CGAL::dynamic_face_property_t<Vector>(), mesh);

    CGAL::Polygon_mesh_processing::compute_normals(mesh, vnormals, fpm);
  }

  Ui_Window ui_window;
  QWidget window(app.activeWindow());
  ui_window.setupUi(&window);
  using Viewer = CGAL::Basic_viewer_qt;
  Viewer* viewer = ui_window.viewer;
  auto threadBox = ui_window.threadBox;
  auto timelineSlider = ui_window.timelineSlider;
  auto timelineBox = ui_window.timelineBox;
  new QShortcut(QKeySequence(Qt::CTRL + Qt::Key_Q), &window,
                SLOT(close()));

  viewer->set_draw_edges(false);
  viewer->set_draw_vertices(true);
  using graph_traits = boost::graph_traits<Mesh>;
  using v_descriptor_t = typename graph_traits::vertex_descriptor;
  std::cerr << "Properties:\n";
  for(auto p: mesh.template properties<v_descriptor_t>()) {
    std::cerr << "  \"" << p << "\"\n";
  }
  auto v_thread_id_pair =
      mesh.template property_map<v_descriptor_t, int>("v:thread_id");
  if(! v_thread_id_pair.second ) {
    std::cerr << "Missing property map v:thread_id\n";
    return EXIT_FAILURE;
  }
  auto v_thread_id = v_thread_id_pair.first;
  auto v_timestamp_pair =
      mesh.template property_map<v_descriptor_t, int>("v:timestamp");
  if(! v_timestamp_pair.second ) {
    std::cerr << "Missing property map v:timestamp\n";
    return EXIT_FAILURE;
  }
  auto v_timestamp = v_timestamp_pair.first;
  {
    auto it_pair = std::minmax_element(vertices(mesh).first,
                                       vertices(mesh).second,
                                       [&](auto v1, auto v2) {
                                         return get(v_timestamp, v1) < get(v_timestamp, v2);
                                       });
    auto it_min = it_pair.first;
    auto it_max = it_pair.second;
    if(it_min != it_max) {
      auto min_ = get(v_timestamp, *it_min);
      auto max_ = get(v_timestamp, *it_max);
      timelineSlider->setMinimum(min_);
      timelineSlider->setMaximum(max_);
      timelineBox->setMinimum(min_);
      timelineBox->setMaximum(max_);
      timelineSlider->setValue(max_);
      timelineBox->setValue(max_);
    }
  }
  {
    auto it_pair = std::minmax_element(vertices(mesh).first,
                                       vertices(mesh).second,
                                       [&](auto v1, auto v2) {
                                         return get(v_thread_id, v1) < get(v_thread_id, v2);
                                       });
    auto it_min = it_pair.first;
    auto it_max = it_pair.second;
    if(it_min != it_max) {
      threadBox->setMinimum(get(v_thread_id, *it_min));
      threadBox->setMaximum(get(v_thread_id, *it_max));
      threadBox->setValue(get(v_thread_id, *it_min));
      std::cerr << "Threads: [ #" << threadBox->minimum()
                << " - #" << threadBox->maximum() << "]\n";
    }
  }
  {
    auto first = vertices(mesh).first;
    auto end = vertices(mesh).second;
    auto filter = [&](auto v) {
      return get(v_thread_id, v)  > 0;
    };
    auto filter_begin = boost::make_filter_iterator(std::cref(filter), first, end);
    auto filter_end  = boost::make_filter_iterator(std::cref(filter), end, end);
    auto it = std::min_element(filter_begin,
                               filter_end,
                               [&](auto v1, auto v2) {
                                 return get(v_timestamp, v1) < get(v_timestamp, v2);
                               });
    if(it != filter_end) {
      const auto min_timestamp_of_non_0_thread = *it;
      timelineSlider->setMinimum(min_timestamp_of_non_0_thread);
      timelineBox->setMinimum(min_timestamp_of_non_0_thread);
      std::cerr << "Timeline: [ " << min_timestamp_of_non_0_thread
                << " - " << timelineBox->maximum() << " ]\n";
    }
  }


  QVariantAnimation rainbow_colors;

  rainbow_colors.setKeyValueAt(0.0f, QVariant{ QColor::fromRgbF(0.75, 0.75, 1.  ) });
  rainbow_colors.setKeyValueAt(0.2f, QVariant{ QColor::fromRgbF(0.  , 0.  , 1.  ) });
  rainbow_colors.setKeyValueAt(0.4f, QVariant{ QColor::fromRgbF(0.  , 1.  , 0.  ) });
  rainbow_colors.setKeyValueAt(0.6f, QVariant{ QColor::fromRgbF(1.  , 1.  , 0.  ) });
  rainbow_colors.setKeyValueAt(0.8f, QVariant{ QColor::fromRgbF(1.  , 0.  , 0.  ) });
  rainbow_colors.setKeyValueAt(1.0f, QVariant{ QColor::fromRgbF(0.5 , 0.  , 0.  ) });
  rainbow_colors.setDuration(1000);

  auto fromQColor = [](QColor c) -> CGAL::IO::Color {
    return {
      (unsigned char)c.red(),
      (unsigned char)c.green(),
      (unsigned char)c.blue()
    };
  };

  auto get_color = [&](double key) {
    rainbow_colors.setCurrentTime(1000 * key);
    return fromQColor( rainbow_colors.currentValue().value<QColor>() );
  };

  auto draw_mesh = [&]() {
    for(auto fh : faces(mesh)) {
      viewer->face_begin(CGAL::IO::gray());
      auto hd=mesh.halfedge(fh);
      do
        {
          auto v = mesh.source(hd);
          viewer->add_point_in_face(mesh.point(v), get(vnormals, v));
          hd=mesh.next(hd);
        }
      while(hd!=mesh.halfedge(fh));
      viewer->face_end();
    }
  };
  bool call_redraw = false;
  std::function<void()> recompute_elements = [&]() {
    viewer->clear();
    draw_mesh();
    const auto thread_id_to_display = threadBox->value();
    const auto max_timestamp = timelineSlider->value();
    std::cerr << "Display thread #" << thread_id_to_display << '\n';
    std::cerr << "Max timestamp: " << max_timestamp << '\n';
    auto min = threadBox->minimum();
    auto delta = threadBox->maximum()-threadBox->minimum();
    for(auto v : vertices(mesh)) {
      const auto ts = get(v_timestamp, v);
      const auto thread_id = get(v_thread_id, v);
      if((thread_id_to_display == 0 || thread_id == thread_id_to_display) &&
         ts < max_timestamp)
      {
        CGAL::IO::Color c = get_color(((0.+thread_id) - min) / delta);
        viewer->add_point(mesh.point(v), c);
      }
    }
    if(call_redraw) viewer->redraw();
  };

  QObject::connect(threadBox, QOverload<int>::of(&QSpinBox::valueChanged),
                   recompute_elements);
  QObject::connect(timelineSlider, &QAbstractSlider::valueChanged,
                   recompute_elements);
  recompute_elements();
  window.show();
  call_redraw = true;
  return app.exec();
}


template <typename C3t3>
auto draw_c3t3_surface(const C3t3& c3t3)
{
  std::map<std::thread::id, int> mapping;
#if CGAL_MESH_3_STATS_THREADS
  {
    int patch_id_nb = 0;
    for(auto facet: CGAL::make_range(c3t3.facets_in_complex_begin(),
                                     c3t3.facets_in_complex_end()))
      {
        auto& cell = facet.first;
        auto& index = facet.second;
        int& patch_id = mapping[cell->thread_ids[index]];
        if(patch_id == 0) patch_id = ++patch_id_nb;
        cell->set_surface_patch_index(index, patch_id);
        facet = c3t3.triangulation().mirror_facet(facet);
        cell->set_surface_patch_index(index, patch_id);
      }
  }
#endif

  using Tr_geom_traits = typename C3t3::Triangulation::Geom_traits;
  using Point_3 = typename Tr_geom_traits::Point_3;
  using Mesh = CGAL::Surface_mesh<Point_3>;
  Mesh surface_mesh;
  CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, surface_mesh);

#if CGAL_MESH_3_STATS_THREADS
  auto point_pmap = get(CGAL::vertex_point, surface_mesh);
  auto thread_pmap =
      surface_mesh
          .template add_property_map<typename Mesh::Vertex_index, int>(
              "v:thread_id")
          .first;
  auto time_pmap =
      surface_mesh
          .template add_property_map<typename Mesh::Vertex_index, int>(
              "v:timestamp")
          .first;

  auto vit = c3t3.triangulation().finite_vertices_begin();
  auto vend = c3t3.triangulation().finite_vertices_end();
  std::unordered_map<Point_3, typename C3t3::Triangulation::Vertex_handle>
      point_to_vertices;
  for (auto v : c3t3.triangulation().all_vertex_handles()) {
    point_to_vertices[v->point().point()] = v;
  }
  for(auto v: vertices(surface_mesh)) {
    CGAL_assertion(vit != vend);
    const auto point = get(point_pmap, v);
    auto vh = point_to_vertices[point];
    put(thread_pmap, v, mapping[vh->thread_id]);
    put(time_pmap, v, vh->time_stamp());
    ++vit;
  }
#endif
  draw_c3t3_surface_from_surface_mesh(surface_mesh);
  return surface_mesh;
}
