#include <QtCore/qglobal.h>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>

#include <CGAL/iterator.h>
#include <CGAL/utility.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/Polygon_mesh_processing/tangential_relaxation.h>

#include "Scene.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_selection_item.h"

#include "ui_Smoothing_plugin.h"
#include "ui_Smoothing_tangential_relaxation.h"

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

    actionRelax_ = new QAction(tr("Tangential Relaxation"), mw);
    actionRelax_->setProperty("subMenuName", "Polygon Mesh Processing");
    connect(actionRelax_, SIGNAL(triggered()), this, SLOT(tangential_relaxation_action()));
  }

  QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionSmoothing_ << actionRelax_;
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

  Ui::Tangential_relaxation_dialog
  relaxation_dialog(QDialog* dialog)
  {
    Ui::Tangential_relaxation_dialog ui;
    ui.setupUi(dialog);
    connect(ui.buttonBox, SIGNAL(accepted()), dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()), dialog, SLOT(reject()));

    ui.nbIterations_spinbox->setSingleStep(1);
    ui.nbIterations_spinbox->setRange(1/*min*/, 1000/*max*/);
    ui.nbIterations_spinbox->setValue(1);

    ui.smooth1D_checkbox->setChecked(true);

    return ui;
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

  void tangential_relaxation_action()
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();

    Scene_facegraph_item* poly_item =
      qobject_cast<Scene_facegraph_item*>(scene->item(index));
    Scene_polyhedron_selection_item* selection_item =
      qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

    if (poly_item || selection_item)
    {
      if (selection_item && selection_item->selected_facets.empty())
      {
        QMessageBox::warning(mw, "Empty Facets", "There are no selected facets. Aborting.");
        return;
      }
      // Create dialog box
      QDialog dialog(mw);
      Ui::Tangential_relaxation_dialog ui = relaxation_dialog(&dialog);

      // Get values
      int i = dialog.exec();
      if (i == QDialog::Rejected)
      {
        std::cout << "Tangential relaxation aborted" << std::endl;
        return;
      }

      unsigned int nb_iter = ui.nbIterations_spinbox->value();
      bool smooth_features = ui.smooth1D_checkbox->isChecked();

      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);
      QElapsedTimer time;
      time.start();

      FaceGraph& pmesh = (poly_item != nullptr)
        ? *poly_item->polyhedron()
        : *selection_item->polyhedron();

      if (selection_item)
      {
        boost::unordered_set<vertex_descriptor> vset;
        for (face_descriptor f : selection_item->selected_facets)
        {
          for(vertex_descriptor fv : CGAL::vertices_around_face(halfedge(f, pmesh), pmesh))
            vset.insert(fv);
        }

        CGAL::Polygon_mesh_processing::tangential_relaxation(
          vset,
          pmesh,
          CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
          .edge_is_constrained_map(selection_item->constrained_edges_pmap())
          .vertex_is_constrained_map(selection_item->constrained_vertices_pmap())
          .relax_constraints(smooth_features));
        selection_item->polyhedron_item()->invalidateOpenGLBuffers();
        Q_EMIT selection_item->polyhedron_item()->itemChanged();
        selection_item->invalidateOpenGLBuffers();
        selection_item->setKeepSelectionValid(Scene_polyhedron_selection_item::None);
      }
      else if (poly_item)
      {
        CGAL::Polygon_mesh_processing::tangential_relaxation(
          vertices(pmesh),
          pmesh,
          CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter));

        poly_item->invalidateOpenGLBuffers();
        Q_EMIT poly_item->itemChanged();
      }

      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
    }

    // default cursor
    QApplication::restoreOverrideCursor();
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
      angle_and_area_smoothing(pmesh, parameters::do_project(projection)
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
        angle_and_area_smoothing(pmesh, parameters::do_project(projection)
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
        angle_and_area_smoothing(selection_item->selected_facets, pmesh, parameters::do_project(projection)
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
  QAction* actionRelax_;
  QDockWidget* dock_widget;
  Ui::Smoothing ui_widget;
};

#include "Smoothing_plugin.moc"
