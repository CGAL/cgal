#include <QApplication>
#include <QAction>
#include <QList>
#include <QMainWindow>
#include <QMessageBox>
#include <QtDebug>

#include "Scene_points_with_normal_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include "Messages_interface.h"

#include <CGAL/Polygon_mesh_processing/distance.h>

using namespace CGAL::Three;
class Polyhedron_demo_point_set_from_sampling_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface, Messages_interface*);

  bool applicable(QAction*) const {
    const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

    return qobject_cast<Scene_surface_mesh_item*>(scene->item(index))
      || qobject_cast<Scene_polygon_soup_item*>(scene->item(index));
  }

  QList<QAction*> actions() const;

public Q_SLOTS:
  void createPointSet();

private:
  CGAL::Three::Scene_interface* scene;
  QAction* actionPointSetFromSampling;


}; // end Polyhedron_demo_point_set_from_sampling_plugin

void Polyhedron_demo_point_set_from_sampling_plugin::init(QMainWindow* mainWindow,
                                                          CGAL::Three::Scene_interface* scene_interface,
                                                          Messages_interface*)
{
  scene = scene_interface;
  actionPointSetFromSampling = new QAction(tr("&Create Point Set from Sampling"), mainWindow);
  actionPointSetFromSampling->setObjectName("actionPointSetFromSampling");
  connect(actionPointSetFromSampling, SIGNAL(triggered()),
          this, SLOT(createPointSet()));
}

QList<QAction*> Polyhedron_demo_point_set_from_sampling_plugin::actions() const {
  return QList<QAction*>() << actionPointSetFromSampling;
}



void Polyhedron_demo_point_set_from_sampling_plugin::createPointSet()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* points = new Scene_points_with_normal_item();
  std::vector<Kernel::Point_3> pts;

  if (points){
    points->setColor(Qt::blue);
  }else{
    return;
  }
  Scene_surface_mesh_item* sm_item =
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

  if (sm_item){
    points->setName(QString("%1 (sampled)").arg(sm_item->name()));
    CGAL::Polygon_mesh_processing::sample_triangle_mesh(*sm_item->polyhedron(),
                                                        std::back_inserter(pts));
  }

  Scene_polygon_soup_item* soup_item =
    qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

  if (soup_item){
    points->setName(QString("%1 (sampled)").arg(soup_item->name()));
    CGAL::Polygon_mesh_processing::sample_triangle_soup(soup_item->points(),
                                                        soup_item->polygons(),
                                                        std::back_inserter(pts));
  }
  for (std::size_t i = 0; i < pts.size(); ++i){
    points->point_set()->insert(pts[i]);
  }
    scene->addItem(points);
  QApplication::restoreOverrideCursor();
}


#include "Point_set_from_sampling_plugin.moc"
