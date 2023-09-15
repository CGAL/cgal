#undef NDEBUG
#include <QtCore/qglobal.h>

#include "Messages_interface.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>
#include "Scene_polyhedron_selection_item.h"
#include "ui_ACVD_simplification_widget.h"

#include "SMesh_type.h"
typedef Scene_surface_mesh_item Scene_facegraph_item;

#include <CGAL/iterator.h>
#include <CGAL/Polygon_mesh_processing/acvd/acvd.h>

#include <QElapsedTimer>
#include <QRadioButton>
#include <QCheckBox>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QDockWidget>
#include <QEvent>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QFileDialog>

#include <vector>
#include <algorithm>
#include <queue>

typedef Scene_facegraph_item::Face_graph FaceGraph;

using namespace CGAL::Three;
class Polyhedron_demo_acvd_simplification_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "acvd_simplification_plugin.json")
public:
  bool applicable(QAction*) const {
    return qobject_cast<Scene_facegraph_item*>(scene->item(scene->mainSelectionIndex()));
  }
  void print_message(QString message) { CGAL::Three::Three::information(message);}
  QList<QAction*> actions() const { return QList<QAction*>() << actionACVD; }


  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface* m) {
    mw = mainWindow;
    scene = scene_interface;
    messages = m;
    actionACVD = new QAction(tr(
                                "ACVD Simplification"
                                ), mw);
    //actionACVD->setProperty("subMenuName", "");
    if (actionACVD)
      connect(actionACVD, SIGNAL(triggered()), this, SLOT(acvd_action()));

    dock_widget = new QDockWidget(
          "ACVD Simplification", mw);

    dock_widget->setVisible(false);

    ui_widget.setupUi(dock_widget);
    ui_widget.nb_clusters_spin_box->setMinimum(10);
    ui_widget.nb_clusters_spin_box->setMaximum(1000000);
    // get number of vertices of the mesh
    Scene_facegraph_item* item = getSelectedItem<Scene_facegraph_item>();
    if(!item) { return; }
    // set default number of clusters to 5% of the number of vertices
    ui_widget.nb_clusters_spin_box->setValue(max(item->face_graph()->number_of_faces() / 20, (unsigned int)10));

    addDockWidget(dock_widget);
    dock_widget->setWindowTitle(tr(
                                  "ACVD Simplification"
                                  ));

    connect(ui_widget.ACVD_button,  SIGNAL(clicked()), this, SLOT(on_ACVD_button_clicked()));
  }
  virtual void closure()
  {
    dock_widget->hide();
  }

public Q_SLOTS:
  void acvd_action() {
    dock_widget->show();
    dock_widget->raise();
  }

  void on_ACVD_button_clicked() {
    /*typename FaceGraphItem::Face_graph* graph = item->face_graph();
    if(!graph) return;
    QElapsedTimer time;
    time.start();
    CGAL::Three::Three::information("Upsample subdivision...");
    QApplication::setOverrideCursor(Qt::WaitCursor);
    CGAL::Subdivision_method_3::Upsample_subdivision(
        *graph,
        CGAL::params::number_of_iterations(ui_widget.nb_clusters_spin_box->value()/10)
    );

    CGAL::Three::Three::information(QString("ok (%1 ms)").arg(time.elapsed()));
    QApplication::restoreOverrideCursor();
    item->invalidateOpenGLBuffers();
    scene->itemChanged(item);*/
  }

private:
  Messages_interface* messages;
  QAction* actionACVD;

  QDockWidget* dock_widget;
  Ui::ACVD_Simplification ui_widget;

}; // end Polyhedron_demo_acvd_simplification_plugin

// Q_EXPORT_PLUGIN2(Polyhedron_demo_acvd_simplification_plugin, Polyhedron_demo_acvd_simplification_plugin)

#include "ACVD_simplification_plugin.moc"
