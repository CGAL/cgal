#include <QtCore/qglobal.h>

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"

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
    // Ok if there's at least one mesh
    Q_FOREACH(int index, scene->selectionIndices())
    {
      if(qobject_cast<Scene_surface_mesh_item*>(scene->item(index)))
        return true;
      if(qobject_cast<Scene_polygon_soup_item*>(scene->item(index)))
        return true;
    }

    return false;
  }

  QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionAlpha_wrap_3_;
  }

private:
  Ui::alpha_wrap_3_dialog create_dialog(QDialog* dialog)
  {
    Ui::alpha_wrap_3_dialog ui;
    ui.setupUi(dialog);
    connect(ui.buttonBox, SIGNAL(accepted()), dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()), dialog, SLOT(reject()));

    return ui;
  }

public Q_SLOTS:
  void on_actionAlpha_wrap_3_triggered()
  {
    using Points = typename Scene_polygon_soup_item::Points;
    using Polygons = typename Scene_polygon_soup_item::Polygons;
    using Oracle = CGAL::Alpha_wraps_3::internal::Triangle_soup_oracle<Points, Polygons>;

    QDialog dialog(mw);
    Ui::alpha_wrap_3_dialog ui = create_dialog(&dialog);
    dialog.setWindowFlags(Qt::Dialog|Qt::CustomizeWindowHint|Qt::WindowCloseButtonHint);

    int i = dialog.exec();
    if(i == QDialog::Rejected)
      return;

    const bool is_relative_alpha = ui.relativeAlpha->isChecked();
    const bool is_relative_offset = ui.relativeOffset->isChecked();
    const bool enforce_manifoldenss = ui.runManifoldness->isChecked();
    double alpha = ui.alphaValue->value();
    double offset = ui.offsetValue->value();

    if(alpha <= 0. || offset <= 0.)
      return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    Oracle oracle;
    Q_FOREACH(int index, scene->selectionIndices())
    {
      Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
      if(sm_item != nullptr)
      {
        if(!is_triangle_mesh(*(sm_item->polyhedron())))
          continue;

        Points points;
        Polygons polygons;
        CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup(*(sm_item->polyhedron()),
                                                                    points, polygons);

        oracle.add_triangle_soup(points, polygons);
      }
      else
      {
        Scene_polygon_soup_item* soup_item = qobject_cast<Scene_polygon_soup_item*>(scene->item(index));
        if(soup_item != nullptr)
        {
          bool is_triangle_soup = true;
          for(auto p : soup_item->polygons())
          {
            if(p.size() != 3)
            {
              is_triangle_soup = false;
              break;
            }
          }

          if(is_triangle_soup)
            oracle.add_triangle_soup(soup_item->points(), soup_item->polygons());
        }
      }
    }

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
        CGAL::parameters::do_enforce_manifoldness(enforce_manifoldenss));

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
