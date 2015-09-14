#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/hierarchical_clustering.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>

#include "ui_Polyhedron_demo_point_set_hierarchical_clustering_plugin.h"

class Polyhedron_demo_point_set_hierarchical_clustering_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  QAction* actionHierarchicalCluster;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    actionHierarchicalCluster = new QAction(tr("Point set hierarchical clustering"), mainWindow);
    actionHierarchicalCluster->setObjectName("actionHierarchicalCluster");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  bool applicable(QAction*) const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionHierarchicalCluster;
  }

public Q_SLOTS:
  void on_actionHierarchicalCluster_triggered();

}; // end Polyhedron_demo_point_set_hierarchical_clustering_plugin

class Point_set_demo_point_set_hierarchical_clustering_dialog : public QDialog, private Ui::PointSetHierarchicalClusteringDialog
{
  Q_OBJECT
public:
  Point_set_demo_point_set_hierarchical_clustering_dialog(QWidget * /*parent*/ = 0)
  {
    setupUi(this);
  }

  int size() const { return m_maximumClusterSize->value(); }
  double var_max() const { return m_maximumSurfaceVariation->value(); }
};

void Polyhedron_demo_point_set_hierarchical_clustering_plugin::on_actionHierarchicalCluster_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
    {
      // Gets point set
      Point_set* points = item->point_set();
      if(points == NULL)
        return;

      // Gets options
      Point_set_demo_point_set_hierarchical_clustering_dialog dialog;
      if(!dialog.exec())
	return;

      QApplication::setOverrideCursor(Qt::WaitCursor);

      CGAL::Timer task_timer; task_timer.start();

      std::cerr << "Hierarchical clustering (cluster size = " << dialog.size()
		<< ", maximum variation = " << dialog.var_max() << ")" << std::endl;



      Scene_points_with_normal_item* new_item = new Scene_points_with_normal_item();
      CGAL::hierarchical_clustering(points->begin(), points->end(),
				    std::back_inserter (*(new_item->point_set())),
				    dialog.size(), dialog.var_max());

      new_item->setName(QString("%1 (hierarchical clustering)").arg(item->name()));
      new_item->set_has_normals (false);
      new_item->setColor(item->color());
      new_item->setRenderingMode(item->renderingMode());
      new_item->setVisible(item->visible());
      scene->addItem(new_item);
      
      std::size_t memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << "Clustering: " << new_item->point_set()->size () << " point(s) generated ("
		<< task_timer.time() << " seconds, "
		<< (memory>>20) << " Mb allocated)"
		<< std::endl;

      QApplication::restoreOverrideCursor();

    }
}

#include "Polyhedron_demo_point_set_hierarchical_clustering_plugin.moc"
