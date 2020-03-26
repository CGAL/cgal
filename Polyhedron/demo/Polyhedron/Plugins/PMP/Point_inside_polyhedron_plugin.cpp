#include <QtCore/qglobal.h>


#include "Messages_interface.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>

#include <CGAL/Timer.h>
#include <CGAL/Random.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include "ui_Point_inside_polyhedron_widget.h"

#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QDockWidget>
#include <QInputDialog>

#include <vector>
#include <algorithm>

#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <boost/optional/optional.hpp>

using namespace CGAL::Three;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;



class Polyhedron_demo_point_inside_polyhedron_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  bool applicable(QAction*) const
  {
    for(CGAL::Three::Scene_interface::Item_id i = 0, end = scene->numberOfEntries();
        i < end; ++i)
    {
      if( qobject_cast<Scene_surface_mesh_item*>(scene->item(i)) != NULL)
        return true;
    }

    //if the loop ends without returning true, return false
    return false;
  }
  void print_message(QString message) { CGAL::Three::Three::information(message); }
  QList<QAction*> actions() const { return QList<QAction*>() << actionPointInsidePolyhedron; }


  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface* m)
  {
    mw = mainWindow;
    scene = scene_interface;
    messages = m;

    actionPointInsidePolyhedron = new QAction(tr("Point Inside Polyhedron"), mw);
    actionPointInsidePolyhedron->setProperty("subMenuName", "Polygon Mesh Processing");
    connect(actionPointInsidePolyhedron, SIGNAL(triggered()), this, SLOT(point_inside_polyhedron_action()));

    dock_widget = new QDockWidget("Point Inside Polyhedron", mw);
    dock_widget->setVisible(false);
    ui_widget.setupUi(dock_widget);

    addDockWidget(dock_widget);

    connect(ui_widget.Select_button,  SIGNAL(clicked()), this, SLOT(on_Select_button()));
    connect(ui_widget.Sample_random_points_from_bbox,  SIGNAL(clicked()), this, SLOT(on_Sample_random_points_from_bbox()));

  }

  virtual void closure()
  {
    dock_widget->hide();
  }


public Q_SLOTS:
  void point_inside_polyhedron_action() {
    dock_widget->show();
    dock_widget->raise();
  }

  void on_Select_button()
  {
    bool inside = ui_widget.Inside_check_box->isChecked();
    bool on_boundary = ui_widget.On_boundary_check_box->isChecked();
    bool outside = ui_widget.Outside_check_box->isChecked();

    if(!(inside || on_boundary || outside)) {
      print_message("Error: please check at least one parameter check box.");
      return;
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    // place all selected polyhedron and point items to vectors below
    std::vector<const SMesh*> smeshs;

    typedef CGAL::Side_of_triangle_mesh<SMesh, Kernel> Point_inside_smesh;
        // it does not support copy-construction so let's use pointers
    std::vector<Point_inside_smesh*>inside_smesh_testers;

    std::vector<Point_set*> point_sets;
    Q_FOREACH(CGAL::Three::Scene_interface::Item_id id, scene->selectionIndices()) {
      Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(id));
      if (sm_item){
        inside_smesh_testers.push_back(new Point_inside_smesh(*(sm_item->polyhedron())));
         smeshs.push_back(sm_item->polyhedron());
      }

      Scene_points_with_normal_item* point_item = qobject_cast<Scene_points_with_normal_item*>(scene->item(id));
      if(point_item) { point_sets.push_back(point_item->point_set()); }
    }


    // there should be at least one selected polyhedron and point item
    if(inside_smesh_testers.empty()) { print_message("Error: there is no selected polyhedron item(s)."); }
    if(point_sets.empty()) {
    if(!generated_points.empty())
      point_sets.push_back(generated_points.last()->point_set());
    else
      print_message("Error: there is no selected point set item(s).");
    }
    if((inside_smesh_testers.empty()) || point_sets.empty()) { QApplication::restoreOverrideCursor(); return; }

    // deselect all points
    for(std::vector<Point_set*>::iterator point_set_it = point_sets.begin();
      point_set_it != point_sets.end(); ++point_set_it) {
      (*point_set_it)->unselect_all();
    }

    CGAL::Timer timer; timer.start();

    print_message(
      QString("Constructing with %1 items is done in %2 sec.").arg(inside_smesh_testers.size()).arg(timer.time()));
    timer.reset();

    std::size_t nb_selected = 0;
    for (Point_set* point_set : point_sets)
      point_set->set_first_selected
        (std::partition
         (point_set->begin(), point_set->end(),
          [&](const Point_set::Index& idx) -> bool
          {
            for (const Point_inside_smesh* tester : inside_smesh_testers)
            {
              CGAL::Bounded_side res = (*tester)(point_set->point(idx));
              if ( (inside      && res == CGAL::ON_BOUNDED_SIDE) ||
                   (on_boundary && res == CGAL::ON_BOUNDARY)     ||
                   (outside     && res == CGAL::ON_UNBOUNDED_SIDE) )
              {
                ++ nb_selected;
                return false;
              }
            }
            return true;
          }));

    print_message(QString("%1 points are selected. All Done!").arg(nb_selected));

    // delete testers
  for (std::size_t i = 0; i < inside_smesh_testers.size(); ++i)
      delete inside_smesh_testers[i];

    bool found = false;
    // for repaint
    Q_FOREACH(CGAL::Three::Scene_interface::Item_id id, scene->selectionIndices()) {
      Scene_points_with_normal_item* point_item = qobject_cast<Scene_points_with_normal_item*>(scene->item(id));
      if(point_item) {
        found = true;
        point_item->invalidateOpenGLBuffers();
        scene->itemChanged(point_item);
      }
    }
    if(!found && !generated_points.empty()) {
      generated_points.last()->invalidateOpenGLBuffers();
      scene->itemChanged(generated_points.last());
    }
    QApplication::restoreOverrideCursor();
  }

  void on_Sample_random_points_from_bbox() {

    // calculate bbox of selected polyhedron items
    boost::optional<CGAL::Three::Scene_interface::Bbox> bbox
      = boost::make_optional(false, CGAL::Three::Scene_interface::Bbox());
    // Workaround a bug in g++-4.8.3:
    //   http://stackoverflow.com/a/21755207/1728537
    // Using boost::make_optional to copy-initialize 'bbox' hides the
    //   warning about '*bbox' not being initialized.
    // -- Laurent Rineau, 2014/10/30

    Q_FOREACH(CGAL::Three::Scene_interface::Item_id id, scene->selectionIndices()) {
      Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(id));
      if(sm_item) {
        if(!bbox) {
          bbox = sm_item->bbox();
        }
        else {
          *bbox = *bbox + sm_item->bbox();
        }
      }
    }

    if(!bbox) {
      print_message("Error: there is no selected polyhedron item(s).");
      return;
    }

    // take number of points param
    bool ok;
    const int nb_points =
      QInputDialog::getInt(mw, tr("Number of Points"),
      tr("Number of Points:"),
      100000, // default value
      1, // min
      (int)1.e9, // max
      10, // step for the spinbox
      &ok);

    if(!ok) { return; }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QApplication::processEvents();
    // sample random points and constuct item
    Scene_points_with_normal_item* point_item = new Scene_points_with_normal_item();
    point_item->setName(QString("sample-%1").arg(nb_points));
    CGAL::Random rg(1340818006);

    double grid_dx = bbox->xmax() - bbox->xmin();
    double grid_dy = bbox->ymax() - bbox->ymin();
    double grid_dz = bbox->zmax() - bbox->zmin();

    for(int i=0; i < nb_points; i++){
      point_item->point_set()->insert(
      Epic_kernel::Point_3(bbox->xmin ()+ rg.get_double()* grid_dx,
        bbox->ymin() + rg.get_double()* grid_dy,
        bbox->zmin() + rg.get_double()* grid_dz)
      );
    }

    scene->addItem(point_item);
    scene->itemChanged(point_item);
    generated_points.append(point_item);
    connect(point_item, SIGNAL(destroyed(QObject*)),
            this, SLOT(resetGeneratedPoints(QObject*)));
    QApplication::restoreOverrideCursor();
  }
private Q_SLOTS:
  void resetGeneratedPoints(QObject* o)
  {
    Q_FOREACH(Scene_points_with_normal_item* item , generated_points)
    if(item == o)
    {
      generated_points.removeAll(item);
    }
  }
private:
  Messages_interface* messages;
  QAction* actionPointInsidePolyhedron;
  QList<Scene_points_with_normal_item*> generated_points;
  QDockWidget* dock_widget;
  Ui::Point_inside_polyhedron ui_widget;

}; // end Polyhedron_demo_point_inside_polyhedron_plugin

#include "Point_inside_polyhedron_plugin.moc"
