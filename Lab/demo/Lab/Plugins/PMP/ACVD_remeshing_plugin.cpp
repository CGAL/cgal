#include <QtCore/qglobal.h>

#include "Messages_interface.h"
#include <CGAL/Three/CGAL_Lab_plugin_helper.h>
#include <CGAL/Three/Three.h>
#include "Scene_polyhedron_selection_item.h"
#include "ui_ACVD_remeshing_dialog.h"

#include "SMesh_type.h"
typedef Scene_surface_mesh_item Scene_facegraph_item;

#include <CGAL/iterator.h>
#include <CGAL/Polygon_mesh_processing/approximated_centroidal_Voronoi_diagram_remeshing.h>

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
namespace PMP = CGAL::Polygon_mesh_processing;

class CGAL_Lab_acvd_remeshing_plugin :
  public QObject,
  public CGAL_Lab_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "acvd_remeshing_plugin.json")
public:
  bool applicable(QAction*) const {
    return qobject_cast<Scene_facegraph_item*>(scene->item(scene->mainSelectionIndex()));
  }
  void print_message(QString message) { CGAL::Three::Three::information(message);}
  QList<QAction*> actions() const { return QList<QAction*>() << actionACVD; }


  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    mw = mainWindow;
    scene = scene_interface;
    actionACVD = new QAction(tr(
                                "ACVD Remeshing"
                                ), mw);
    //actionACVD->setProperty("subMenuName", "");
    if (actionACVD)
      connect(actionACVD, SIGNAL(triggered()), this, SLOT(acvd_action()));
  }

public Q_SLOTS:
  void acvd_action() {

    // Create dialog box
    QDialog dialog(mw);
    ui.setupUi(&dialog);

    ui.nb_clusters_spin_box->setMinimum(10);
    ui.nb_clusters_spin_box->setMaximum(1000000);
    // get number of vertices of the mesh
    Scene_facegraph_item* item = getSelectedItem<Scene_facegraph_item>();
    if(!item) { return; }
    // set default number of clusters to 5% of the number of vertices
    ui.nb_clusters_spin_box->setValue((std::max)(item->face_graph()->number_of_faces() / 20, (unsigned int)10));

    int i = dialog.exec();
    if (i == QDialog::Rejected)
    {
      std::cout << "Remeshing aborted" << std::endl;
      return;
    }


    Face_graph* graph = item->face_graph();
    if(!graph) return;
    QElapsedTimer time;
    time.start();
    CGAL::Three::Three::information("ACVD Remeshing...");


    QApplication::setOverrideCursor(Qt::WaitCursor);
    FaceGraph remeshed = *graph;

    // TODO add post-processing

    if (ui.IsotropicClustering->isChecked())
      PMP::approximated_centroidal_Voronoi_diagram_remeshing(remeshed, ui.nb_clusters_spin_box->value());
    else if (ui.QEMClustering->isChecked())
      PMP::approximated_centroidal_Voronoi_diagram_remeshing(remeshed, ui.nb_clusters_spin_box->value(), CGAL::parameters::use_qem_based_energy(true));
    else
      PMP::approximated_centroidal_Voronoi_diagram_remeshing(remeshed, ui.nb_clusters_spin_box->value(), CGAL::parameters::gradation_factor(0.8));

    // update the scene
    CGAL::Three::Three::information(QString("ok (%1 ms)").arg(time.elapsed()));
    QApplication::restoreOverrideCursor();
    item->invalidateOpenGLBuffers();

    // add the simplified mesh to the scene
    Scene_facegraph_item* remeshed_item = new Scene_facegraph_item(remeshed);
    remeshed_item->setName(QString(item->name() + "_simplified"));
    scene->addItem(remeshed_item);
    remeshed_item->invalidateOpenGLBuffers();
    scene->itemChanged(remeshed_item);
    item->setVisible(false);
  }

private:
  QAction* actionACVD;
  Ui::ACVDDialog ui;
}; // end CGAL_Lab_acvd_remeshing_plugin


#include "ACVD_remeshing_plugin.moc"
