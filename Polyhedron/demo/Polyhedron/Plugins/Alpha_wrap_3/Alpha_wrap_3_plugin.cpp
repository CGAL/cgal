#include <QtCore/qglobal.h>

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Scene_polylines_item.h"
#include "Scene_points_with_normal_item.h"

#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>

#include <QElapsedTimer>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QString>
#include <QDialog>
#include <QMessageBox>
#include <QtPlugin>

#include <vector>
#include <algorithm>
#include <queue>
#include <sstream>
#include <stdexcept>

#include "ui_alpha_wrap_3_dialog.h"

template <typename Scene, typename WrapItem>
struct AW3_visu_visitor
{
  AW3_visu_visitor(Scene* scene,
                   WrapItem* wrap_item)
    : m_scene(scene),
      m_wrap_item(wrap_item)
  { }

public:
  template <typename AlphaWrapper, typename Point>
  void before_Steiner_point_insertion(const AlphaWrapper& wrapper,
                                      const Point& p)
  {
    // @todo something nicer & smarter
    if(wrapper.triangulation().number_of_vertices() < 500)
    {
      // Do all iterations
    }
    else if(wrapper.triangulation().number_of_vertices() < 5000)
    {
      if(wrapper.triangulation().number_of_vertices() % 10 != 0)
        return;
    }
    else if(wrapper.triangulation().number_of_vertices() < 10000)
    {
      if(wrapper.triangulation().number_of_vertices() % 100 != 0)
        return;
    }
    else if(wrapper.triangulation().number_of_vertices() < 1000000)
    {
      if(wrapper.triangulation().number_of_vertices() % 1000 != 0)
        return;
    }

    if(!m_wrap_item)
      return;

//    for(auto cit=wrapper.triangulation().finite_cells_begin(), cend=wrapper.triangulation().finite_cells_end(); cit!=cend; ++cit)
//      std::cout << cit->time_stamp() << std::endl;

    std::cout << "on_Steiner_point_insert (" << wrapper.triangulation().number_of_vertices() << ")" << std::endl;

    // EXTRACT ---

    using Dt = typename std::decay<decltype(wrapper.triangulation())>::type;
    using Vertex_handle = typename Dt::Vertex_handle;
    using Facet = typename Dt::Facet;
    using Cell_handle = typename Dt::Cell_handle;

    std::vector<Kernel::Point_3> points;
    std::vector<std::vector<std::size_t> > faces;

    std::unordered_map<Vertex_handle, std::size_t> vertex_to_id;
    std::size_t nv = 0;

    std::size_t min_time_stamp = -1, max_time_stamp = 0;
    for(auto cit=wrapper.triangulation().finite_cells_begin(), cend=wrapper.triangulation().finite_cells_end(); cit!=cend; ++cit)
    {
      if(cit->time_stamp() > max_time_stamp)
        max_time_stamp = cit->time_stamp();
      if(cit->time_stamp() < min_time_stamp)
        min_time_stamp = cit->time_stamp();
    }

    std::vector<CGAL::IO::Color> vcolors;
    std::vector<CGAL::IO::Color> fcolors;

    std::size_t vi = 0, fi = 0;
    for(auto fit=wrapper.triangulation().finite_facets_begin(), fend=wrapper.triangulation().finite_facets_end(); fit!=fend; ++fit)
    {
      Facet f = *fit;
      if(!f.first->info().is_outside)
        f = wrapper.triangulation().mirror_facet(f);

      const Cell_handle c = f.first;
      const int s = f.second;
      const Cell_handle nh = c->neighbor(s);
      if(c->info().is_outside == nh->info().is_outside)
        continue;

      std::array<std::size_t, 3> ids;
      for(int pos=0; pos<3; ++pos)
      {
        Vertex_handle vh = c->vertex(Dt::vertex_triple_index(s, pos));
        auto insertion_res = vertex_to_id.emplace(vh, nv);
        if(insertion_res.second) // successful insertion, never-seen-before vertex
        {
          points.push_back(wrapper.triangulation().point(vh));
          vcolors.push_back(CGAL::IO::Color(0, 0, 0));
          ++nv;
        }

        ids[pos] = insertion_res.first->second;
      }

      faces.emplace_back(std::vector<std::size_t>{ids[0], ids[1], ids[2]});
      double color_val = double(c->time_stamp() - min_time_stamp) / double(max_time_stamp - min_time_stamp);
      color_val = int(256. * color_val);
      // fcolors.push_back(CGAL::IO::Color(color_val, 10, 150)); // young is red, old is blue
      // fcolors.push_back(CGAL::IO::Color(256 - color_val, 256 - color_val, 256 - color_val)); // young is light, old is dark
      fcolors.push_back(CGAL::IO::Color(100, 100, 100)); // darkish gray
    }

//    std::cout << "DT NV " << wrapper.triangulation().number_of_vertices() << " vs " << nv << " vs " << points.size() << std::endl;
//    std::cout << "DT NF " << wrapper.triangulation().number_of_finite_facets() << " vs " << fi << " vs " << faces.size() << std::endl;

    // --- EXTRACT

    m_wrap_item->load(points, faces, fcolors, vcolors);
    m_wrap_item->setAlpha(255 / 2);

    m_wrap_item->invalidateOpenGLBuffers();
    m_wrap_item->redraw();
    m_wrap_item->itemChanged();

//    SMesh wrap;
//    auto vpm = get(CGAL::vertex_point, wrap);
//    wrapper.extract_surface(wrap, vpm, true /*tolerate non manifoldness*/);

//    Scene_surface_mesh_item* new_wrap_item = new Scene_surface_mesh_item(wrap);
//    new_wrap_item->setName(QString("Wrap"));
//    new_wrap_item->setColor(Qt::gray);

//    new_wrap_item->setName(m_wrap_item->name());
//    new_wrap_item->setColor(m_wrap_item->color());
//    new_wrap_item->setRenderingMode(m_wrap_item->renderingMode());
//    new_wrap_item->setVisible(m_wrap_item->visible());
//    Scene_item_with_properties *property_item = dynamic_cast<Scene_item_with_properties*>(new_wrap_item);
//    m_scene->replaceItem(m_scene->item_id(m_wrap_item), new_wrap_item, true);
//    if(property_item)
//      property_item->copyProperties(m_wrap_item);
//    new_wrap_item->invalidateOpenGLBuffers();
//    m_wrap_item->deleteLater();

    QApplication::processEvents();
  }

  template <typename AlphaWrapper, typename VertexHandle>
  void after_Steiner_point_insertion(const AlphaWrapper& wrapper,
                                     const VertexHandle vh)
  {

  }

  void after_alpha_wrapping()
  {
    m_wrap_item->setAlpha(255);
    m_wrap_item->invalidateOpenGLBuffers();
    m_wrap_item->redraw();
    m_wrap_item->itemChanged();
    QApplication::processEvents();
  }

private:
  Scene* m_scene;
  WrapItem* m_wrap_item;
};

class Polyhedron_demo_alpha_wrap_3_plugin
  : public QObject,
    public CGAL::Three::Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

private:
  Ui::alpha_wrap_3_dialog ui;

public:
  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;

    actionAlpha_wrap_3_ = new QAction("3D Alpha Wrapping", this->mw);
    if(actionAlpha_wrap_3_)
      connect(actionAlpha_wrap_3_, SIGNAL(triggered()), this, SLOT(on_actionAlpha_wrap_3_triggered()));
  }

  bool applicable(QAction*) const
  {
    Q_FOREACH(int index, scene->selectionIndices())
    {
      if(qobject_cast<Scene_polygon_soup_item*>(scene->item(index)))
        return true;
      if(qobject_cast<Scene_surface_mesh_item*>(scene->item(index)))
        return true;
      if(qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index)))
        return true;
      if(qobject_cast<Scene_polylines_item*>(scene->item(index)))
        return true;
      if(qobject_cast<Scene_points_with_normal_item*>(scene->item(index)))
        return true;
    }

    return false;
  }

  QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionAlpha_wrap_3_;
  }

private:
  void print_message(QString message) const
  {
    CGAL::Three::Three::information(message);
  }

  Ui::alpha_wrap_3_dialog create_dialog(QDialog* dialog)
  {
    ui.setupUi(dialog);

    connect(ui.wrapEdges, SIGNAL(clicked(bool)), this, SLOT(toggle_wrap_faces()));
    connect(ui.wrapFaces, SIGNAL(clicked(bool)), this, SLOT(toggle_wrap_edges()));

    connect(ui.buttonBox, SIGNAL(accepted()), dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()), dialog, SLOT(reject()));

    return ui;
  }

public Q_SLOTS:
  void toggle_wrap_faces()
  {
    if(!ui.wrapEdges->isChecked()) // if edges are disabled, so are faces
      ui.wrapFaces->setChecked(false);
  }

  void toggle_wrap_edges()
  {
    if(ui.wrapFaces->isChecked()) // if faces are enabled, so are edges
      ui.wrapEdges->setChecked(true);
  }

  void on_actionAlpha_wrap_3_triggered()
  {
    using Triangles = std::vector<Kernel::Triangle_3>;
    using Segments = std::vector<Kernel::Segment_3>;
    using Points = std::vector<Kernel::Point_3>;

    using TS_Oracle = CGAL::Alpha_wraps_3::internal::Triangle_soup_oracle<Kernel>;
    using SS_Oracle = CGAL::Alpha_wraps_3::internal::Segment_soup_oracle<Kernel, TS_Oracle>;
    using Oracle = CGAL::Alpha_wraps_3::internal::Point_set_oracle<Kernel, SS_Oracle>;

    QDialog dialog(mw);
    ui = create_dialog(&dialog);
    dialog.setWindowFlags(Qt::Dialog|Qt::CustomizeWindowHint|Qt::WindowCloseButtonHint);

    int i = dialog.exec();
    if(i == QDialog::Rejected)
      return;

    const bool is_relative_alpha = ui.relativeAlpha->isChecked();
    const bool is_relative_offset = ui.relativeOffset->isChecked();
    const bool enforce_manifoldness = ui.runManifoldness->isChecked();
    double alpha = ui.alphaValue->value();
    double offset = ui.offsetValue->value();

    if(alpha <= 0. || offset <= 0.)
      return;

//    QApplication::setOverrideCursor(Qt::WaitCursor);

    TS_Oracle ts_oracle;
    SS_Oracle ss_oracle(ts_oracle);
    Oracle oracle(ss_oracle);

    Triangles triangles;
    Segments segments;
    Points points;

    Q_FOREACH(int index, scene->selectionIndices())
    {
      // ---
      Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
      if(sm_item != nullptr)
      {
        SMesh* pMesh = sm_item->polyhedron();
        auto vpm = get(CGAL::vertex_point, *pMesh);

        triangles.reserve(triangles.size() + num_faces(*pMesh));
        for(boost::graph_traits<SMesh>::face_descriptor f : faces(*pMesh))
        {
          boost::graph_traits<SMesh>::halfedge_descriptor h = halfedge(f, *pMesh);
          if(!is_triangle(h, *pMesh))
          {
            print_message("Warning: non-triangular face in input");
            continue;
          }

          triangles.emplace_back(get(vpm, target(h, *pMesh)),
                                 get(vpm, target(next(h, *pMesh), *pMesh)),
                                 get(vpm, source(h, *pMesh)));
        }

        sm_item->setRenderingMode(Flat);

        continue;
      }

      // ---
      Scene_polygon_soup_item* soup_item = qobject_cast<Scene_polygon_soup_item*>(scene->item(index));
      if(soup_item != nullptr)
      {
        triangles.reserve(triangles.size() + soup_item->polygons().size());

        for(const auto& p : soup_item->polygons())
        {
          if(p.size() != 3)
          {
            print_message("Warning: non-triangular face in input");
            continue;
          }

          triangles.emplace_back(soup_item->points()[p[0]],
                                 soup_item->points()[p[1]],
                                 soup_item->points()[p[2]]);
        }

        soup_item->setRenderingMode(Flat);

        continue;
      }

      // ---
      Scene_polyhedron_selection_item* selection_item =
          qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));
      if(selection_item != nullptr)
      {
        SMesh* pMesh = selection_item->polyhedron();
        auto vpm = get(CGAL::vertex_point, *pMesh);

        triangles.reserve(triangles.size() + selection_item->selected_facets.size());
        for(const auto& f : selection_item->selected_facets)
        {
          boost::graph_traits<SMesh>::halfedge_descriptor h = halfedge(f, *pMesh);
          if(!is_triangle(h, *pMesh))
          {
            print_message("Warning: non-triangular face in input");
            continue;
          }

          triangles.emplace_back(get(vpm, target(h, *pMesh)),
                                 get(vpm, target(next(h, *pMesh), *pMesh)),
                                 get(vpm, source(h, *pMesh)));
        }

        segments.reserve(segments.size() + selection_item->selected_edges.size());
        for(const auto& e : selection_item->selected_edges)
        {
          segments.emplace_back(get(vpm, target(halfedge(e, *pMesh), *pMesh)),
                                get(vpm, target(opposite(halfedge(e, *pMesh), *pMesh), *pMesh)));
        }

        points.reserve(points.size() + selection_item->selected_vertices.size());
        for(const auto& v : selection_item->selected_vertices)
        {
          points.push_back(get(vpm, v));
        }

        continue;
      }

      // ---
      Scene_polylines_item* lines_item = qobject_cast<Scene_polylines_item*>(scene->item(index));
      if(lines_item != nullptr)
      {
        for(auto& polyline : lines_item->polylines)
        {
          for(auto it=std::begin(polyline), end=std::end(polyline); it!=end; ++it)
          {
            auto nit = std::next(it);
            if(nit != end)
              segments.emplace_back(*it, *nit);
          }
        }

        continue;
      }

      // ---
      Scene_points_with_normal_item* pts_item = qobject_cast<Scene_points_with_normal_item*>(scene->item(index));
      if(pts_item != nullptr)
      {
        points.insert(std::cend(points),
                      std::cbegin(pts_item->point_set()->points()),
                      std::cend(pts_item->point_set()->points()));

        continue;
      }
    }

    const bool wrap_triangles = ui.wrapFaces->isChecked();
    const bool wrap_segments = ui.wrapEdges->isChecked();

    if(!wrap_triangles)
    {
      segments.reserve(segments.size() + 3 * triangles.size());
      for(const auto& tr : triangles)
      {
        segments.emplace_back(tr[0], tr[1]);
        segments.emplace_back(tr[1], tr[2]);
        segments.emplace_back(tr[2], tr[0]);
      }
    }

    if(!wrap_segments)
    {
      points.reserve(points.size() + 2 * segments.size());
      for(const auto& s : segments)
      {
        points.push_back(s[0]);
        points.push_back(s[1]);
      }
    }

    std::cout << triangles.size() << " triangles" << std::endl;
    std::cout << segments.size() << " edges" << std::endl;
    std::cout << points.size() << " points" << std::endl;
    std::cout << "wrap s/tr: " << wrap_segments << " " << wrap_triangles << std::endl;

    if(wrap_triangles)
      oracle.add_triangle_soup(triangles);
    if(wrap_segments)
      oracle.add_segment_soup(segments);
    oracle.add_point_set(points);

    if(!oracle.do_call())
    {
      print_message("Warning: empty input - nothing to wrap");
//      QApplication::restoreOverrideCursor();
      return;
    }

    // Oracles set up, time to wrap

    CGAL::Bbox_3 bbox = oracle.bbox();
    const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                         CGAL::square(bbox.ymax() - bbox.ymin()) +
                                         CGAL::square(bbox.zmax() - bbox.zmin()));

    if(is_relative_alpha)
      alpha = diag_length / alpha;
    if(is_relative_offset)
      offset = diag_length / offset;

    CGAL::Alpha_wraps_3::internal::Alpha_wrap_3<Oracle> aw3(oracle);

    Scene_polygon_soup_item* iterative_wrap_item = new Scene_polygon_soup_item();
    iterative_wrap_item->setName(tr("Iterative wrap").arg(alpha).arg(offset));
    scene->addItem(iterative_wrap_item);

    AW3_visu_visitor<std::remove_reference<decltype(*(this->scene))>::type,
                     Scene_polygon_soup_item> visitor(scene, iterative_wrap_item);

    SMesh wrap;
    aw3(alpha, offset, wrap,
        CGAL::parameters::do_enforce_manifoldness(enforce_manifoldness)
                         .visitor(visitor));

//    Scene_surface_mesh_item* wrap_item = new Scene_surface_mesh_item(wrap);
//    wrap_item->setName(tr("Wrap alpha %2 offset %3").arg(alpha).arg(offset));
//    wrap_item->setColor(Qt::gray);
//    scene->addItem(wrap_item);

//    QApplication::restoreOverrideCursor();
  }

private:
  QAction* actionAlpha_wrap_3_;
};

#include "Alpha_wrap_3_plugin.moc"
