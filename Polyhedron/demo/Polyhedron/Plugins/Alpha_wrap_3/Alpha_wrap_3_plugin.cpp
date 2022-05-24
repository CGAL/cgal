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

    QApplication::setOverrideCursor(Qt::WaitCursor);

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

    if(wrap_triangles)
      oracle.add_triangle_soup(triangles);
    if(wrap_segments)
      oracle.add_segment_soup(segments);
    oracle.add_point_set(points);

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

    SMesh wrap;
    aw3(alpha, offset, wrap,
        CGAL::parameters::do_enforce_manifoldness(enforce_manifoldness));

    Scene_surface_mesh_item* wrap_item = new Scene_surface_mesh_item(wrap);
    wrap_item->setName(tr("Wrap alpha %2 offset %3").arg(alpha).arg(offset));
    wrap_item->setColor(Qt::cyan);
    scene->addItem(wrap_item);

    QApplication::restoreOverrideCursor();
  }

private:
  QAction* actionAlpha_wrap_3_;
};

#include "Alpha_wrap_3_plugin.moc"
