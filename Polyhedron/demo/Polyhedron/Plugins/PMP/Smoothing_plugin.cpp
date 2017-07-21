
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

#include "ui_Smoothing_plugin.h"


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

        connect(ui_widget.Compatible_button,  SIGNAL(clicked()), this, SLOT(on_Compatible_button_clicked()));
        connect(ui_widget.Curvature_button,  SIGNAL(clicked()), this, SLOT(on_Curvature_button_clicked()));

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

public Q_SLOTS:
    void smoothing_action()
    {
        dock_widget->show();
        dock_widget->raise();
    }

    void on_Compatible_button_clicked()
    {
        const Scene_interface::Item_id index = scene->mainSelectionIndex();
        Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));

        Polyhedron& pmesh = *poly_item->polyhedron();
        CGAL::Polygon_mesh_processing::compatible_remeshing(pmesh);

        poly_item->invalidateOpenGLBuffers();
        Q_EMIT poly_item->itemChanged();
        QApplication::restoreOverrideCursor();
    }

    void on_Curvature_button_clicked()
    {
        const Scene_interface::Item_id index = scene->mainSelectionIndex();
        Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));

        Polyhedron& pmesh = *poly_item->polyhedron();
        CGAL::Polygon_mesh_processing::curvature_flow(pmesh);

        poly_item->invalidateOpenGLBuffers();
        Q_EMIT poly_item->itemChanged();
        QApplication::restoreOverrideCursor();
    }


private:
    QAction* actionSmoothing_;
    QDockWidget* dock_widget;
    Ui::Smoothing ui_widget;



};



#include "Smoothing_plugin.moc"


















