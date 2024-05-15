
#include <CGAL/Three/Three.h>
#include <QApplication>
#include <QAction>
#include <QList>
#include <QMainWindow>
#include <QMessageBox>
#include <QtDebug>
#include <QInputDialog>

#include "Scene_points_with_normal_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include "Messages_interface.h"

#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/boost/graph/helpers.h>

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
  actionPointSetFromSampling = new QAction(tr("Create Point Set from Sampling"), mainWindow);
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


  if (points){
    points->setColor(Qt::blue);
  }else{
    QApplication::restoreOverrideCursor();
    return;
  }
  Scene_surface_mesh_item* sm_item =
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

  if (sm_item){
    if(! CGAL::is_triangle_mesh(*sm_item->polyhedron())){
      CGAL::Three::Three::error(QString("The mesh must have triangle faces"));
      QApplication::restoreOverrideCursor();
      return;
    }
    int nf = num_faces(*sm_item->polyhedron());

    bool ok;
    int nb = 0;
    nb = QInputDialog::getInt(QApplication::activeWindow(), "Sampling",
                              "Number of sample points:",
                              nf , 0, (std::numeric_limits<int>::max)(), 1, &ok);

    points->setName(QString("%1 (sampled)").arg(sm_item->name()));
    if( ok & (nb > 0)){
      points->point_set()->reserve(nb);
      CGAL::Polygon_mesh_processing::sample_triangle_mesh(*sm_item->polyhedron(),
                                                          points->point_set()->point_back_inserter(),
                                                          CGAL::parameters::number_of_points_on_faces(nb)
                                                          .do_sample_vertices(false)
                                                          .do_sample_edges(false));
      scene->addItem(points);
    }
  }

  Scene_polygon_soup_item* soup_item =
    qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

  if (soup_item){
    int nf = static_cast<int>(soup_item->polygons().size());

    for(const auto& f : soup_item->polygons()){
      if(f.size() != 3){
        CGAL::Three::Three::error(QString("The polygons must be triangles"));
        QApplication::restoreOverrideCursor();
        return;
      }
    }

    bool ok;
    int nb = 0;
    nb = QInputDialog::getInt(QApplication::activeWindow(), "Sampling",
                              "Number of sample points:",
                              nf , 0, (std::numeric_limits<int>::max)(), 1, &ok);
    points->setName(QString("%1 (sampled)").arg(soup_item->name()));
    if( ok & (nb > 0)){
      points->point_set()->reserve(nb);
      CGAL::Polygon_mesh_processing::sample_triangle_soup(soup_item->points(),
                                                          soup_item->polygons(),
                                                          points->point_set()->point_back_inserter(),
                                                          CGAL::parameters::number_of_points_on_faces(nb)
                                                          .do_sample_vertices(false)
                                                          .do_sample_edges(false));
      scene->addItem(points);
    }
  }


  QApplication::restoreOverrideCursor();
}


#include "Point_set_from_sampling_plugin.moc"
