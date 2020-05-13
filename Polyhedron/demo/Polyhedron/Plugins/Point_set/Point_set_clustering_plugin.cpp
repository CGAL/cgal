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
#include <QCheckBox>

#include <QMultipleInputDialog.h>

#include "run_with_qprogressdialog.h"

struct Clustering_functor
  : public Functor_with_signal_callback
{
  Point_set* points;
  Point_set::Property_map<std::size_t> cluster_map;
  const double neighbor_radius;
  boost::shared_ptr<std::size_t> result;

  Clustering_functor (Point_set* points,
                      const double neighbor_radius,
                      Point_set::Property_map<std::size_t> cluster_map)
    : points (points), cluster_map (cluster_map),
      neighbor_radius (neighbor_radius),
      result (new std::size_t(0)) { }

  void operator()()
  {
    *result = CGAL::cluster_point_set (*points, cluster_map,
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

  bool applicable(QAction*) const {
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
    QDoubleSpinBox* neighbor_radius = dialog.add<QDoubleSpinBox> ("Neighbor radius (0 = automatic):");
    neighbor_radius->setRange (0, 10000000);
    neighbor_radius->setValue (0);
    QSpinBox* min_nb = dialog.add<QSpinBox> ("Minimum number of points per cluster:");
    min_nb->setRange (1, 10000000);
    min_nb->setValue (1);

    QCheckBox* add_property = dialog.add<QCheckBox> ("Add a \"cluster\" property to the input item");
    add_property->setChecked (true);

    QCheckBox* gen_color = dialog.add<QCheckBox> ("Generate one colored point set");
    gen_color->setChecked (true);

    QCheckBox* gen_sub = dialog.add<QCheckBox> ("Generate N point subsets");
    gen_sub->setChecked (false);

    if (!dialog.exec())
      return;

    QApplication::setOverrideCursor(Qt::BusyCursor);
    QApplication::processEvents();
    CGAL::Real_timer task_timer; task_timer.start();

    Point_set::Property_map<std::size_t> cluster_map;

    if (add_property->isChecked())
      cluster_map = points->add_property_map<std::size_t> ("cluster_map").first;
    else
      // Use long name to avoid overwriting potentially existing map
      cluster_map = points->add_property_map<std::size_t> ("cluster_point_set_property_map").first;

    // Default value
    if (neighbor_radius->value() == 0)
    {
      neighbor_radius->setRange (-1, 10000000);
      neighbor_radius->setValue(-1);
    }

    // Computes average spacing
    Clustering_functor functor (points, neighbor_radius->value(), cluster_map);
    run_with_qprogressdialog (functor, "Clustering...", mw);

    std::size_t nb_clusters = *functor.result;

    Scene_group_item* group = nullptr;
    std::vector<Scene_points_with_normal_item*> new_items;

    if (gen_sub->isChecked())
    {
      group = new Scene_group_item(QString("%1 (clusters)").arg(item->name()));
      scene->addItem(group);
      new_items.reserve (nb_clusters);
      for (std::size_t i = 0; i < nb_clusters; ++ i)
      {
        Scene_points_with_normal_item* new_item = new Scene_points_with_normal_item;
        new_item->point_set()->copy_properties (*points);
        CGAL::Random rand((unsigned int)(i));
        unsigned char r, g, b;
        r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
        g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
        b = static_cast<unsigned char>(64 + rand.get_int(0, 192));
        new_item->setRgbColor(r, g, b);
        new_item->setName (QString("Cluster %1 of %2").arg(i).arg(item->name()));
        new_items.push_back (new_item);
      }
    }

    std::vector<std::size_t> cluster_size (nb_clusters, 0);
    for (Point_set::Index idx : *points)
    {
      if (gen_sub->isChecked())
        new_items[cluster_map[idx]]->point_set()->insert (*points, idx);
      cluster_size[cluster_map[idx]] ++;
    }

    if (gen_color->isChecked())
    {
      Scene_points_with_normal_item* colored;
      Point_set::Property_map<unsigned char> red, green, blue;

      colored = new Scene_points_with_normal_item;
      colored->setName (QString("%1 (clustering)").arg(item->name()));

      red = colored->point_set()->add_property_map<unsigned char>("red", 0).first;
      green = colored->point_set()->add_property_map<unsigned char>("green", 0).first;
      blue = colored->point_set()->add_property_map<unsigned char>("blue", 0).first;
      colored->point_set()->check_colors();

      colored->point_set()->reserve (points->size());

      for (Point_set::Index idx : *points)
      {
        Point_set::Index iidx = *(colored->point_set()->insert (points->point(idx)));
        if (cluster_size[cluster_map[idx]] >= std::size_t(min_nb->value()))
        {
          CGAL::Random rand((unsigned int)(cluster_map[idx] + 1));
          unsigned char r, g, b;
          r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
          g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
          b = static_cast<unsigned char>(64 + rand.get_int(0, 192));
          red[iidx] = r;
          green[iidx] = g;
          blue[iidx] = b;
        }
      }
      scene->addItem(colored);
    }

    if (gen_sub->isChecked())
    {
      for (Scene_points_with_normal_item* new_item : new_items)
      {
        if (new_item->point_set()->size() >= std::size_t(min_nb->value()))
        {
          scene->addItem(new_item);
          scene->changeGroup (new_item, group);
        }
        else
          delete new_item;
      }
    }

    if (!add_property->isChecked())
      points->remove_property_map (cluster_map);

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
