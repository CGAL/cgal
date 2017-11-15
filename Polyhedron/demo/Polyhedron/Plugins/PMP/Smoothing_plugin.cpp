
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

#include <CGAL/iterator.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/utility.h>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/property_map.h>

#include <CGAL/Polygon_mesh_processing/smoothing.h>
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

        connect(ui_widget.Apply_button,  SIGNAL(clicked()), this, SLOT(on_Apply_by_type_clicked()));
        connect(ui_widget.MCF_Button,  SIGNAL(clicked()), this, SLOT(on_Apply_MCF_clicked()));
        connect(ui_widget.modified_MCF_button ,  SIGNAL(clicked()), this, SLOT(on_Apply_mMCF_clicked()));
        connect(ui_widget.Run_convergence_button,  SIGNAL(clicked()), this, SLOT(on_Run_convergence_clicked()));


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
        ui_widget.Angle_spinBox->setValue(1);
        ui_widget.Angle_spinBox->setSingleStep(1);
        ui_widget.Angle_spinBox->setMinimum(1);

        ui_widget.Area_spinBox->setValue(1);
        ui_widget.Area_spinBox->setSingleStep(1);
        ui_widget.Area_spinBox->setMinimum(1);

        /*
        ui_widget.gd_dSpinBox->setSingleStep(0.0001);
        ui_widget.gd_dSpinBox->setDecimals(4);
        ui_widget.gd_dSpinBox->setMinimum(0.0001);
        ui_widget.gd_dSpinBox->setValue(0.001);

        ui_widget.use_weights_checkBox->setChecked(true);

        ui_widget.dist_dSpinBox->setValue(0.01);
        ui_widget.dist_dSpinBox->setSingleStep(0.0001);
        ui_widget.dist_dSpinBox->setDecimals(4);
        ui_widget.dist_dSpinBox->setMinimum(0.0001);

        ui_widget.gd_precision_label->setToolTip("Tradeoff between precision and speed. Less is more precise.");
        ui_widget.distance_label->setToolTip("Tradeoff between precision and speed. Less is more precise.");
        */

        ui_widget.iterations_spinBox->setValue(20);
        ui_widget.iterations_spinBox->setSingleStep(1);
        ui_widget.iterations_spinBox->setMinimum(1);

        /*
        ui_widget.curv_iterations_spinBox->setValue(1);
        ui_widget.curv_iterations_spinBox->setSingleStep(1);
        ui_widget.curv_iterations_spinBox->setMinimum(1);

        ui_widget.curv_iterations_spinBox_2->setValue(1);
        ui_widget.curv_iterations_spinBox_2->setSingleStep(1);
        ui_widget.curv_iterations_spinBox_2->setMinimum(1);
        */
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

    void on_Apply_by_type_clicked()
    {
        const Scene_interface::Item_id index = scene->mainSelectionIndex();

        Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));

        Scene_polyhedron_selection_item* selection_item =
          qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

        Polyhedron& pmesh = (poly_item != NULL) ?
                    * poly_item->polyhedron() :
                    * selection_item->polyhedron();

        QApplication::setOverrideCursor(Qt::WaitCursor);

        if(ui_widget.Angle_checkBox->isChecked())
        {
            unsigned int nb_iter = ui_widget.Angle_spinBox->value();
            //bool use_weights = ui_widget.use_weights_checkBox->isChecked();
            angle_smoothing(pmesh, parameters::number_of_iterations(nb_iter));

            poly_item->invalidateOpenGLBuffers();
            Q_EMIT poly_item->itemChanged();
        }

        if(ui_widget.Area_checkBox->isChecked())
        {
            unsigned int nb_iter = ui_widget.Area_spinBox->value();
            //double gd_precision = ui_widget.gd_dSpinBox->value();
            area_smoothing(pmesh, parameters::number_of_iterations(nb_iter));

            poly_item->invalidateOpenGLBuffers();
            Q_EMIT poly_item->itemChanged();
        }

        QApplication::restoreOverrideCursor();
    }

    void on_Apply_MCF_clicked()
    {
        const Scene_interface::Item_id index = scene->mainSelectionIndex();
        Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
        Scene_polyhedron_selection_item* selection_item =
          qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

        Polyhedron& pmesh = (poly_item != NULL) ?
                    * poly_item->polyhedron() :
                    * selection_item->polyhedron();

        QApplication::setOverrideCursor(Qt::WaitCursor);

        //unsigned int nb_iter = ui_widget.curv_iterations_spinBox->value();
        unsigned int nb_iter = 1;
		unsigned int itime = ui_widget.time_spinBox->value();
		const double time = itime * 1e-6;
		std::cout << "time: " << time << std::endl;
		curvature_flow_smoothing(pmesh, time, parameters::number_of_iterations(nb_iter));

        poly_item->invalidateOpenGLBuffers();
        Q_EMIT poly_item->itemChanged();

        QApplication::restoreOverrideCursor();
    }

    void on_Apply_mMCF_clicked()
    {
        const Scene_interface::Item_id index = scene->mainSelectionIndex();
        Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
        Scene_polyhedron_selection_item* selection_item =
          qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

        Polyhedron& pmesh = (poly_item != NULL) ?
                    * poly_item->polyhedron() :
                    * selection_item->polyhedron();

		unsigned int itime = ui_widget.time_spinBox->value();
		const double time = itime * 1e-6;

        QApplication::setOverrideCursor(Qt::WaitCursor);

        //unsigned int nb_iter = ui_widget.curv_iterations_spinBox_2->value();
        unsigned int nb_iter = 1;

        if(!is_stiffness_matrix_setup)
        {
            setup_mcf_system(pmesh, nb_iter, stiffness_matrix);
            is_stiffness_matrix_setup = true;
        }
        // todo: pass nb_iter as named parameter
        solve_mcf_system(pmesh, time, nb_iter, stiffness_matrix);


        poly_item->invalidateOpenGLBuffers();
        Q_EMIT poly_item->itemChanged();

        QApplication::restoreOverrideCursor();
    }

    void on_Run_convergence_clicked()
    {
        const Scene_interface::Item_id index = scene->mainSelectionIndex();
        Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
        Scene_polyhedron_selection_item* selection_item =
          qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

        Polyhedron& pmesh = (poly_item != NULL) ?
                    * poly_item->polyhedron() :
                    * selection_item->polyhedron();

        QApplication::setOverrideCursor(Qt::WaitCursor);

        unsigned int nb_iter = ui_widget.iterations_spinBox->value();
        double dist = ui_widget.dist_dSpinBox->value();
        //double gd_precision = ui_widget.gd_dSpinBox->value();
        //bool use_weights = ui_widget.use_weights_checkBox->isChecked();

        compatible_smoothing(pmesh,
                             parameters::number_of_iterations(nb_iter).
                             distance_precision(dist));

        poly_item->invalidateOpenGLBuffers();
        Q_EMIT poly_item->itemChanged();

        QApplication::restoreOverrideCursor();
    }



private:
    QAction* actionSmoothing_;
    QDockWidget* dock_widget;
    Ui::Smoothing ui_widget;

    Eigen::SparseMatrix<double> stiffness_matrix;
    bool is_stiffness_matrix_setup = false;



};



#include "Smoothing_plugin.moc"


















