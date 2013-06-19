#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "ui_MCFSkeleton_dialog.h"

#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"

#include <QApplication>
#include <QMainWindow>
#include <QInputDialog>
#include <QTime>

#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>

class Polyhedron_demo_mean_curvature_flow_skeleton_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionMCFSkeleton";
  }

  bool applicable() const {
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }
public slots:
  void on_actionMCFSkeleton_triggered();

}; // end Polyhedron_demo_mean_curvature_flow_skeleton_plugin

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionMCFSkeleton_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    Polyhedron* pMesh = item->polyhedron();

    if(!pMesh) return;

    QDialog dialog(mw);
    Ui::MCFSkeleton_dialog ui;
    ui.setupUi(&dialog);
    connect(ui.buttonBox, SIGNAL(accepted()),
            &dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()),
            &dialog, SLOT(reject()));
    double diag = scene->len_diagonal();

    ui.omega_L->setValue(1);
    ui.omega_H->setValue(0.1);
    ui.edgelength_TH->setDecimals(7);
    ui.edgelength_TH->setValue(0.002 * diag);
    ui.alpha->setValue(0.15);
    ui.zero_TH->setDecimals(8);
    ui.zero_TH->setValue(1e-07);

    int i = dialog.exec();
    if(i == QDialog::Rejected)
      return;

    const double omega_L = ui.omega_L->value();
    const double omega_H = ui.omega_H->value();
    const double edgelength_TH = ui.edgelength_TH->value();
    const double alpha = ui.alpha->value();
    const double zero_TH = ui.zero_TH->value();

    // skeletonize
    QTime time;
    time.start();
    std::cout << "Skeletonize...";
    QApplication::setOverrideCursor(Qt::WaitCursor);



    std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

    // update scene
    scene->itemChanged(index);
    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_mean_curvature_flow_skeleton_plugin, Polyhedron_demo_mean_curvature_flow_skeleton_plugin)

#include "Polyhedron_demo_mean_curvature_flow_skeleton_plugin.moc"
