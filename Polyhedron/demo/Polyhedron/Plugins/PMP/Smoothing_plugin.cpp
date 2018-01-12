
#include <QtCore/qglobal.h>
#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QDockWidget>
#include <QString>
#include <QInputDialog>
#include <QtPlugin>
#include <QMessageBox>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>

#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Polyhedron_type.h"
#include "Scene.h"

#include <CGAL/iterator.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/utility.h>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/property_map.h>

#include <CGAL/Polygon_mesh_processing/mesh_smoothing.h>
#include <CGAL/Polygon_mesh_processing/shape_smoothing.h>

#include "ui_Smoothing_plugin.h"


using namespace CGAL::Polygon_mesh_processing;
using namespace CGAL::Three;
class Polyhedron_demo_smothing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
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

    connect(ui_widget.angle_button,  SIGNAL(clicked()), this, SLOT(on_angles_clicked()));
    connect(ui_widget.area_button,  SIGNAL(clicked()), this, SLOT(on_areas_clicked()));
    connect(ui_widget.curvature_button,  SIGNAL(clicked()), this, SLOT(on_curvature_clicked()));
    connect(ui_widget.modified_curv_button,  SIGNAL(clicked()), this, SLOT(on_modified_curv_clicked()));

}

    QList<QAction*> actions() const
    {
        return QList<QAction*>() << actionSmoothing_;
    }

    bool applicable(QAction*) const
    {
      const Scene_interface::Item_id index = scene->mainSelectionIndex();
      if (qobject_cast<Scene_polyhedron_item*>(scene->item(index)))
        return true;
      else if (qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index)))
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
      // todo: better time step
      ui_widget.time_step_spinBox->setSingleStep(0.00001);
      ui_widget.time_step_spinBox->setMinimum(1e-6);

      QObject* scene_object = dynamic_cast<QObject*>(scene);
      connect(scene_object, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)),
              this, SLOT(on_actionItemAboutToBeDestroyed()));
    }

public Q_SLOTS:
    void smoothing_action()
    {
        dock_widget->show();
        dock_widget->raise();

        // needed?
        const Scene_interface::Item_id index = scene->mainSelectionIndex();
        Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));

        Scene_polyhedron_selection_item* selection_item =
          qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

        if(poly_item || selection_item)
        {
            init_ui();
        }
    }

    void on_angles_clicked()
    {
        const Scene_interface::Item_id index = scene->mainSelectionIndex();
        Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
        Scene_polyhedron_selection_item* selection_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));
        Polyhedron& pmesh = (poly_item != NULL) ? * poly_item->polyhedron() : * selection_item->polyhedron();

        QApplication::setOverrideCursor(Qt::WaitCursor);

        unsigned int nb_iter = ui_widget.angles_iter_spinBox->value();
        smooth_angles(pmesh, parameters::number_of_iterations(nb_iter));

        poly_item->invalidateOpenGLBuffers();
        poly_item->itemChanged();
        QApplication::restoreOverrideCursor();

     }

    void on_areas_clicked()
    {
        const Scene_interface::Item_id index = scene->mainSelectionIndex();
        Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
        Scene_polyhedron_selection_item* selection_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));
        Polyhedron& pmesh = (poly_item != NULL) ? * poly_item->polyhedron() : * selection_item->polyhedron();

        QApplication::setOverrideCursor(Qt::WaitCursor);

        unsigned int nb_iter = ui_widget.areas_iter_spinBox->value();
        smooth_areas(pmesh, parameters::number_of_iterations(nb_iter));

        poly_item->invalidateOpenGLBuffers();
        poly_item->itemChanged();
        QApplication::restoreOverrideCursor();
     }

    void on_curvature_clicked()
    {
        const Scene_interface::Item_id index = scene->mainSelectionIndex();
        Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
        Scene_polyhedron_selection_item* selection_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));
        Polyhedron& pmesh = (poly_item != NULL) ? * poly_item->polyhedron() : * selection_item->polyhedron();

        QApplication::setOverrideCursor(Qt::WaitCursor);

        unsigned int nb_iter = ui_widget.curvature_iter_spinBox ->value();
        smooth_curvature_flow_explicit(faces(pmesh), pmesh, parameters::number_of_iterations(nb_iter));

        poly_item->invalidateOpenGLBuffers();
        poly_item->itemChanged();
        QApplication::restoreOverrideCursor();
    }

    void on_modified_curv_clicked()
    {
        const Scene_interface::Item_id index = scene->mainSelectionIndex();
        Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
        Scene_polyhedron_selection_item* selection_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));
        Polyhedron& pmesh = (poly_item != NULL) ? * poly_item->polyhedron() : * selection_item->polyhedron();

        int index_id = scene->item_id(poly_item);
        const double time_step = ui_widget.time_step_spinBox->value();

        QApplication::setOverrideCursor(Qt::WaitCursor);

        if(index_id != last_index_id)
        {
            setup_mcf_system(faces(pmesh), pmesh, stiffness_matrix_, parameters::all_default());
            last_index_id = index_id;
        }

        solve_mcf_system(faces(pmesh), pmesh, time_step, stiffness_matrix_, parameters::all_default());

        // todo: update scene for selection item
        poly_item->compute_bbox();
        static_cast<Scene*>(scene)->updated_bbox(true);
        poly_item->invalidateOpenGLBuffers();
        poly_item->itemChanged();

        QApplication::restoreOverrideCursor();

    }

    void on_actionItemAboutToBeDestroyed()
    {
        last_index_id= -1;
    }

private:
    QAction* actionSmoothing_;
    QDockWidget* dock_widget;
    Ui::Smoothing ui_widget;

    int last_index_id = -1;
    Eigen::SparseMatrix<double> stiffness_matrix_;


};



#include "Smoothing_plugin.moc"
