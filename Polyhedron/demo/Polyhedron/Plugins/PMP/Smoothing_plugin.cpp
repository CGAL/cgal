#include <QtCore/qglobal.h>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>

#include <CGAL/iterator.h>
#include <CGAL/utility.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>

#include "Scene.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_selection_item.h"

#include "ui_Smoothing_plugin.h"

#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QDockWidget>
#include <QString>
#include <QInputDialog>
#include <QtPlugin>
#include <QMessageBox>

typedef Scene_surface_mesh_item                                     Scene_face_graph_item;
typedef Scene_face_graph_item::Face_graph                           Face_graph;

using namespace CGAL::Polygon_mesh_processing;
using namespace CGAL::Three;

class Polyhedron_demo_smothing_plugin
  : public QObject, public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*)
  {
    scene = scene_interface;
    mw = mainWindow;

    actionSmoothing_ = new QAction(tr("Smoothing"), mw);
    actionSmoothing_->setProperty("subMenuName", "Polygon Mesh Processing");

    connect(actionSmoothing_, SIGNAL(triggered()), this, SLOT(smoothing_action()));

    dock_widget = new QDockWidget("Smoothing", mw);
    dock_widget->setVisible(false);

    ui_widget.setupUi(dock_widget);
    addDockWidget(dock_widget);

    connect(ui_widget.angle_button,  SIGNAL(clicked()), this, SLOT(on_angle_smoothing_clicked()));
    connect(ui_widget.area_button,  SIGNAL(clicked()), this, SLOT(on_area_smoothing_clicked()));
    connect(ui_widget.curvature_flow_button,  SIGNAL(clicked()), this, SLOT(on_curvature_flow_clicked()));
  }

  QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionSmoothing_;
  }

  bool applicable(QAction*) const
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    if(qobject_cast<Scene_face_graph_item*>(scene->item(index)))
      return true;
    else if(qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index)))
      return true;
    else
      return false;
  }

  virtual void closure()
  {
    dock_widget->hide();
  }

  void init_ui()
  {
    ui_widget.angles_iter_spinBox ->setValue(1);
    ui_widget.angles_iter_spinBox->setSingleStep(1);
    ui_widget.angles_iter_spinBox->setMinimum(1);

    ui_widget.areas_iter_spinBox->setValue(1);
    ui_widget.areas_iter_spinBox->setSingleStep(1);
    ui_widget.areas_iter_spinBox->setMinimum(1);

    ui_widget.time_step_spinBox->setValue(0.00001);
    ui_widget.time_step_spinBox->setSingleStep(0.00001);
    ui_widget.time_step_spinBox->setMinimum(1e-6);

    // todo: replace this spinbox with a sliding bar
    ui_widget.precision_spinBox->setValue(0.00001);
    ui_widget.precision_spinBox->setSingleStep(0.0000001);
    ui_widget.precision_spinBox->setMinimum(1e-7);

    ui_widget.projection_checkBox->setChecked(true);

    ui_widget.explicit_checkBox->setChecked(false);

    ui_widget.curvature_iter_spinBox->setValue(1);
  }

public Q_SLOTS:
  void smoothing_action()
  {
    dock_widget->show();
    dock_widget->raise();

    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_face_graph_item* poly_item = qobject_cast<Scene_face_graph_item*>(scene->item(index));

    Scene_polyhedron_selection_item* selection_item =
        qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

    if(poly_item || selection_item)
    {
      init_ui();
    }
  }

  void on_angle_smoothing_clicked()
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_face_graph_item* poly_item = qobject_cast<Scene_face_graph_item*>(scene->item(index));
    Scene_polyhedron_selection_item* selection_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));
    Face_graph& pmesh = (poly_item != nullptr) ? * poly_item->polyhedron() : * selection_item->polyhedron();

    const unsigned int nb_iter = ui_widget.angles_iter_spinBox->value();
    bool projection = ui_widget.projection_checkBox->isChecked();
    QApplication::setOverrideCursor(Qt::WaitCursor);

    if(poly_item)
    {
      smooth_angles(pmesh,
                    parameters::number_of_iterations(nb_iter).do_project(projection));

      poly_item->invalidateOpenGLBuffers();
      poly_item->itemChanged();
    }

    else if(selection_item)
    {
      smooth_angles(selection_item->selected_facets, pmesh,
                    parameters::number_of_iterations(nb_iter)
                    .do_project(projection));

      selection_item->poly_item_changed();
      selection_item->changed_with_poly_item();
    }

    QApplication::restoreOverrideCursor();
  }

  void on_area_smoothing_clicked()
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_face_graph_item* poly_item = qobject_cast<Scene_face_graph_item*>(scene->item(index));
    Scene_polyhedron_selection_item* selection_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));
    Face_graph& pmesh = (poly_item != nullptr) ? * poly_item->polyhedron() : * selection_item->polyhedron();

    unsigned int nb_iter = ui_widget.areas_iter_spinBox->value();
    bool projection = ui_widget.projection_checkBox->isChecked();
    QApplication::setOverrideCursor(Qt::WaitCursor);

    if(poly_item)
    {
      smooth_areas(pmesh, parameters::number_of_iterations(nb_iter).do_project(projection));

      poly_item->invalidateOpenGLBuffers();
      poly_item->itemChanged();
    }
    else if(selection_item)
    {
      smooth_areas(selection_item->selected_facets, pmesh, parameters::number_of_iterations(nb_iter)
                                                                      .do_project(projection));

      selection_item->poly_item_changed();
      selection_item->changed_with_poly_item();
    }
    else
    {
      std::cerr << "Something's gone wrong.\n";
      CGAL_assertion(false);
    }

    QApplication::restoreOverrideCursor();
  }

  void on_curvature_flow_clicked()
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_face_graph_item* poly_item = qobject_cast<Scene_face_graph_item*>(scene->item(index));
    Scene_polyhedron_selection_item* selection_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));
    Face_graph& pmesh = (poly_item != nullptr) ? * poly_item->polyhedron() : * selection_item->polyhedron();

    int index_id = scene->item_id(poly_item);

    const double time_step = ui_widget.time_step_spinBox->value();
    const unsigned int nb_iter = ui_widget.curvature_iter_spinBox->value();

    const bool use_constrained_vertex_map = ui_widget.border_button->isChecked() && !CGAL::is_closed(pmesh);
    const bool use_explicit = ui_widget.explicit_checkBox->isChecked();

    QApplication::setOverrideCursor(Qt::WaitCursor);

    if(poly_item)
    {
      if(use_constrained_vertex_map)
      {
        smooth_along_curvature_flow(pmesh, time_step, parameters::use_explicit_scheme(use_explicit)
                                                                 .number_of_iterations(nb_iter)
                                                                 .vertex_is_constrained_map(get_border_constrained_map(pmesh)));
      }
      else
      {
        smooth_along_curvature_flow(pmesh, time_step, parameters::use_explicit_scheme(use_explicit)
                                                                 .number_of_iterations(nb_iter));
      }

      poly_item->invalidateOpenGLBuffers();
      poly_item->itemChanged();
    }
    else if(selection_item)
    {
      if(use_constrained_vertex_map)
      {
        smooth_along_curvature_flow(selection_item->selected_facets,
                                    pmesh, time_step, parameters::use_explicit_scheme(use_explicit)
                                                                 .number_of_iterations(nb_iter)
                                                                 .vertex_is_constrained_map(get_border_constrained_map(pmesh)));
      }
      else
      {
        smooth_along_curvature_flow(selection_item->selected_facets,
                                    pmesh, time_step, parameters::use_explicit_scheme(use_explicit)
                                                                 .number_of_iterations(nb_iter));
      }

      selection_item->poly_item_changed();
      selection_item->changed_with_poly_item();
    }
    else
    {
      std::cerr << "Something's gone wrong.\n";
      CGAL_assertion(false);
    }

    // recenter scene
    //poly_item->compute_bbox();
    //static_cast<Scene*>(scene)->updated_bbox(true);
    QApplication::restoreOverrideCursor();
  }

private:
  typedef boost::graph_traits<Face_graph>::vertex_descriptor vertex_descriptor;

  template<class TriangleMesh>
  typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<bool> >::type
  get_border_constrained_map(TriangleMesh& tm)
  {
    typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor         vertex_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor       halfedge_descriptor;

    typedef CGAL::dynamic_vertex_property_t<bool>                                 Vertex_bool_property;
    typename boost::property_map<TriangleMesh, Vertex_bool_property>::type vcm = get(Vertex_bool_property(), tm);

    for(vertex_descriptor v : vertices(tm))
      put(vcm, v, false);

    for(halfedge_descriptor h : halfedges(tm))
    {
      if(CGAL::is_border(h, tm))
        put(vcm, target(h, tm), true);
    }

    return vcm;
  }

private:
  QAction* actionSmoothing_;
  QDockWidget* dock_widget;
  Ui::Smoothing ui_widget;
};

#include "Smoothing_plugin.moc"
