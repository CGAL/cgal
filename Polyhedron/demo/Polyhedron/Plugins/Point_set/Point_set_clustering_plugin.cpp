#include "config.h"
#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Scene_group_item.h>

#include <CGAL/cluster_point_set.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Random.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>
#include <QSpinBox>
#include <QDoubleSpinBox>

#include <QMultipleInputDialog.h>

#include "run_with_qprogressdialog.h"

struct Clustering_functor
  : public Functor_with_signal_callback
{
  Point_set* points;
  Point_set::Property_map<std::size_t> cluster_map;
  const int nb_neighbors;
  const double neighbor_radius;
  boost::shared_ptr<std::size_t> result; 

  Clustering_functor (Point_set* points, const int nb_neighbors,
                      const double neighbor_radius,
                      Point_set::Property_map<std::size_t> cluster_map)
    : points (points), cluster_map (cluster_map),
      nb_neighbors (nb_neighbors), neighbor_radius (neighbor_radius),
      result (new std::size_t(0)) { }

  void operator()()
  {
    *result = CGAL::cluster_point_set (*points, cluster_map, nb_neighbors,
                                       points->parameters().neighbor_radius(neighbor_radius).
                                       callback(*(this->callback())));
  }
};

using namespace CGAL::Three;

class Polyhedron_demo_point_set_clustering_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  QAction* actionCluster;

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {

    scene = scene_interface;
    mw = mainWindow;
    actionCluster = new QAction(tr("Cluster Point Set"), mainWindow);
    actionCluster->setObjectName("actionCluster");
    actionCluster->setProperty("subMenuName","Point Set Processing");
    autoConnectActions();
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionCluster;
  }

  bool applicable(QAction* action) const {
    Scene_points_with_normal_item* item = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    return item;
  }

public Q_SLOTS:
  void on_actionCluster_triggered();

}; // end 

void Polyhedron_demo_point_set_clustering_plugin::on_actionCluster_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
  {
    // Gets point set
    Point_set* points = item->point_set();
    if(points == NULL)
        return;

    QMultipleInputDialog dialog ("Clustering", mw);
    QSpinBox* nb_neighbors = dialog.add<QSpinBox> ("Number of neighbors (0 = use radius):");
    nb_neighbors->setRange (0, 10000000);
    nb_neighbors->setValue (12);
    QDoubleSpinBox* neighbor_radius = dialog.add<QDoubleSpinBox> ("Neighbor radius (0 = use number):");
    neighbor_radius->setRange (0, 10000000);
    neighbor_radius->setValue (0);
    QSpinBox* min_nb = dialog.add<QSpinBox> ("Minimum number of points per cluster:");
    min_nb->setRange (1, 10000000);
    min_nb->setValue (1);

    if (!dialog.exec())
      return;
    
    QApplication::setOverrideCursor(Qt::BusyCursor);
    QApplication::processEvents();
    CGAL::Real_timer task_timer; task_timer.start();

    Point_set::Property_map<std::size_t>
      cluster_map = points->add_property_map<std::size_t> ("cluster_point_set_property_map").first;
    
    // Computes average spacing
    Clustering_functor functor (points, nb_neighbors->value(), neighbor_radius->value(), cluster_map);
    run_with_qprogressdialog (functor, "Clustering...", mw);
    
    std::size_t nb_clusters = *functor.result;
 
    CGAL::Random rand(static_cast<unsigned int>(time(0)));

    Scene_group_item* group = new Scene_group_item(QString("%1 (clusters)").arg(item->name()));
    scene->addItem(group);
    
    std::vector<Scene_points_with_normal_item*> new_items;
    new_items.reserve (nb_clusters);
    for (std::size_t i = 0; i < nb_clusters; ++ i)
    {
      Scene_points_with_normal_item* new_item = new Scene_points_with_normal_item;
      new_item->point_set()->copy_properties (*points);
      unsigned char r, g, b;
      r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      b = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      new_item->setRgbColor(r, g, b);
      new_item->setName (QString("Cluster %1 of %2").arg(i).arg(item->name()));
      new_items.push_back (new_item);
    }

    for (Point_set::Index idx : *points)
      new_items[cluster_map[idx]]->point_set()->insert (*points, idx);

    for (Scene_points_with_normal_item* new_item : new_items)
    {
      if (new_item->point_set()->size() >= min_nb->value())
      {
        scene->addItem(new_item);
        scene->changeGroup (new_item, group);
      }
      else
        delete new_item;
    }
    
    std::size_t memory = CGAL::Memory_sizer().virtual_size();
    std::cerr << "Number of clusters = " << nb_clusters << " ("
              << task_timer.time() << " seconds, "
              << (memory>>20) << " Mb allocated)"
              << std::endl;
    QApplication::restoreOverrideCursor();

    item->setVisible (false);
    item->invalidateOpenGLBuffers();
    scene->itemChanged(item);
  }
}



#include "Point_set_clustering_plugin.moc"
