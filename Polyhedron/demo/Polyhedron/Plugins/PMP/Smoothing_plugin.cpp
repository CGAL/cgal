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

#include <QElapsedTimer>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QDockWidget>
#include <QString>
#include <QInputDialog>
#include <QtPlugin>
#include <QMessageBox>

using namespace CGAL::Polygon_mesh_processing;
using namespace CGAL::Three;

typedef Scene_surface_mesh_item                                     Scene_face_graph_item;
typedef Scene_face_graph_item::Face_graph                           Face_graph;

typedef boost::graph_traits<Face_graph>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Face_graph>::edge_descriptor            edge_descriptor;
typedef boost::graph_traits<Face_graph>::face_descriptor            face_descriptor;

typedef CGAL::dynamic_vertex_property_t<bool>                       Vertex_bool_property;
typedef typename boost::property_map<Face_graph, Vertex_bool_property>::type VCMap;

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

    actionSmoothing_ = new QAction(tr("Mesh and Shape Smoothing"), mw);
    actionSmoothing_->setProperty("subMenuName", "Polygon Mesh Processing");

    connect(actionSmoothing_, SIGNAL(triggered()), this, SLOT(smoothing_action()));

    dock_widget = new QDockWidget("Smoothing", mw);
    dock_widget->setVisible(false);

    ui_widget.setupUi(dock_widget);
    addDockWidget(dock_widget);

    connect(ui_widget.area_smoothing_checkBox, SIGNAL(toggled(bool)),
            ui_widget.flip_checkBox, SLOT(setEnabled(bool)));

    connect(ui_widget.mesh_smoothing_button,  SIGNAL(clicked()), this, SLOT(on_mesh_smoothing_clicked()));
    connect(ui_widget.shape_smoothing_button,  SIGNAL(clicked()), this, SLOT(on_shape_smoothing_clicked()));
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
    ui_widget.time_step_spinBox->setValue(0.00001);
    ui_widget.time_step_spinBox->setMinimum(1e-6);

    ui_widget.smooth_iter_spinBox->setValue(1);
    ui_widget.projection_checkBox->setChecked(true);

    ui_widget.area_smoothing_checkBox->setChecked(false);
    ui_widget.flip_checkBox->setDisabled(true);
  }

  void mark_border_vertices(const VCMap vcmap, const Face_graph& pmesh) const
  {
    for(halfedge_descriptor h : halfedges(pmesh))
    {
      if(CGAL::is_border(h, pmesh))
        put(vcmap, target(h, pmesh), true);
    }
  }

  void mark_selected_vertices(const VCMap vcmap,
                              const Face_graph& pmesh,
                              Scene_polyhedron_selection_item* selection_item) const
  {
    for(vertex_descriptor v : selection_item->selected_vertices)
      put(vcmap, v, true);

    for(edge_descriptor e : selection_item->selected_edges)
    {
      put(vcmap, source(e, pmesh), true);
      put(vcmap, target(e, pmesh), true);
    }
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

  void on_mesh_smoothing_clicked()
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_face_graph_item* poly_item = qobject_cast<Scene_face_graph_item*>(scene->item(index));
    Scene_polyhedron_selection_item* selection_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

    if(!poly_item && !selection_item)
      return;

    Face_graph& pmesh = (poly_item != nullptr) ? * poly_item->polyhedron() : * selection_item->polyhedron();

    const unsigned int nb_iter = ui_widget.smooth_iter_spinBox->value();
    const bool projection = ui_widget.projection_checkBox->isChecked();
    const bool use_safety_measures = ui_widget.sanity_checkBox->isChecked();
    const bool constrain_border_vertices = ui_widget.border_button->isChecked() && !CGAL::is_closed(pmesh);
    const bool use_angle_smoothing = ui_widget.angle_smoothing_checkBox->isChecked();
    const bool use_area_smoothing = ui_widget.area_smoothing_checkBox->isChecked();
    const bool use_Delaunay_flips = ui_widget.flip_checkBox->isChecked();

    QApplication::setOverrideCursor(Qt::WaitCursor);

    VCMap vcmap = get(Vertex_bool_property(), pmesh);
    for(vertex_descriptor v : vertices(pmesh))
      put(vcmap, v, false);

    if(constrain_border_vertices)
      mark_border_vertices(vcmap, pmesh);

    if(poly_item)
    {
      smooth_mesh(pmesh, parameters::do_project(projection)
                                    .number_of_iterations(nb_iter)
                                    .vertex_is_constrained_map(vcmap)
                                    .use_safety_constraints(use_safety_measures)
                                    .use_angle_smoothing(use_angle_smoothing)
                                    .use_area_smoothing(use_area_smoothing)
                                    .use_Delaunay_flips(use_Delaunay_flips));

      poly_item->invalidateOpenGLBuffers();
      poly_item->itemChanged();
    }
    else if(selection_item)
    {
      mark_selected_vertices(vcmap, pmesh, selection_item);

      // No faces selected --> use all faces
      if(std::begin(selection_item->selected_facets) == std::end(selection_item->selected_facets))
      {
        smooth_mesh(pmesh, parameters::do_project(projection)
                                      .number_of_iterations(nb_iter)
                                      .vertex_is_constrained_map(vcmap)
                                      .edge_is_constrained_map(selection_item->constrained_edges_pmap())
                                      .use_safety_constraints(use_safety_measures)
                                      .use_angle_smoothing(use_angle_smoothing)
                                      .use_area_smoothing(use_area_smoothing)
                                      .use_Delaunay_flips(use_Delaunay_flips));
      }
      else // some faces exist in the selection
      {
        smooth_mesh(selection_item->selected_facets, pmesh, parameters::do_project(projection)
                                                                       .number_of_iterations(nb_iter)
                                                                       .vertex_is_constrained_map(vcmap)
                                                                       .edge_is_constrained_map(selection_item->constrained_edges_pmap())
                                                                       .use_safety_constraints(use_safety_measures)
                                                                       .use_angle_smoothing(use_angle_smoothing)
                                                                       .use_area_smoothing(use_area_smoothing)
                                                                       .use_Delaunay_flips(use_Delaunay_flips));
      }

      selection_item->poly_item_changed();
      selection_item->changed_with_poly_item();
    }

    QApplication::restoreOverrideCursor();
  }

  void on_shape_smoothing_clicked()
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_face_graph_item* poly_item = qobject_cast<Scene_face_graph_item*>(scene->item(index));
    Scene_polyhedron_selection_item* selection_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

    if(!poly_item && !selection_item)
      return;

    Face_graph& pmesh = (poly_item != nullptr) ? * poly_item->polyhedron() : * selection_item->polyhedron();

    const double time_step = ui_widget.time_step_spinBox->value();
    const unsigned int nb_iter = ui_widget.smooth_iter_spinBox->value();

    const bool constrain_border_vertices = ui_widget.border_button->isChecked() && !CGAL::is_closed(pmesh);

    QApplication::setOverrideCursor(Qt::WaitCursor);

    VCMap vcmap = get(Vertex_bool_property(), pmesh);
    for(vertex_descriptor v : vertices(pmesh))
      put(vcmap, v, false);

    if(constrain_border_vertices)
      mark_border_vertices(vcmap, pmesh);

    if(poly_item)
    {
      smooth_shape(pmesh, time_step, parameters::number_of_iterations(nb_iter)
                                                .vertex_is_constrained_map(vcmap));

      poly_item->invalidateOpenGLBuffers();
      poly_item->itemChanged();
    }
    else if(selection_item)
    {
      mark_selected_vertices(vcmap, pmesh, selection_item);

      if(std::begin(selection_item->selected_facets) == std::end(selection_item->selected_facets))
      {
        smooth_shape(pmesh, time_step, parameters::number_of_iterations(nb_iter)
                                                  .vertex_is_constrained_map(vcmap));
      }
      else
      {
        smooth_shape(selection_item->selected_facets, pmesh, time_step,
                     parameters::number_of_iterations(nb_iter)
                                .vertex_is_constrained_map(vcmap));
      }

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

private:
  QAction* actionSmoothing_;
  QDockWidget* dock_widget;
  Ui::Smoothing ui_widget;
};

#include "Smoothing_plugin.moc"
