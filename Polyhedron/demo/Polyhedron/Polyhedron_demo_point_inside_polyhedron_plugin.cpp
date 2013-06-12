#include <QtCore/qglobal.h>

#include "Messages_interface.h"
#include "Scene_polyhedron_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_interface.h"

#include "Polyhedron_demo_plugin_interface.h"
#include "Polyhedron_type.h"

#include <CGAL/Point_inside_polyhedron_3.h>
#include "ui_Point_inside_polyhedron_widget.h"

#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QDockWidget>

#include <vector>
#include <algorithm>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;

class Polyhedron_demo_point_inside_polyhedron_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  bool applicable() const 
  {
    bool poly_item_exists = false;
    bool point_item_exists = false;

    for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); 
      i < end && (!poly_item_exists || !point_item_exists); ++i)
    {
      poly_item_exists |= qobject_cast<Scene_polyhedron_item*>(scene->item(i)) != NULL;
      point_item_exists |= qobject_cast<Scene_points_with_normal_item*>(scene->item(i)) != NULL;
    }

    return poly_item_exists && point_item_exists;
  }
  void print_message(QString message) { messages->information(message); }
  QList<QAction*> actions() const { return QList<QAction*>() << actionPointInsidePolyhedron; }
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* m)
  {
    mw = mainWindow;
    scene = scene_interface;
    messages = m;

    actionPointInsidePolyhedron = new QAction(tr("Point Inside Polyhedron"), mw);
    connect(actionPointInsidePolyhedron, SIGNAL(triggered()), this, SLOT(point_inside_polyhedron_action()));

    dock_widget = new QDockWidget("Point Inside Polyhedron", mw);
    dock_widget->setVisible(false);
    ui_widget = new Ui::Point_inside_polyhedron();
    ui_widget->setupUi(dock_widget);
    mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

    connect(ui_widget->Select_button,  SIGNAL(clicked()), this, SLOT(on_Select_button())); 
  }

public slots:
  void point_inside_polyhedron_action() { dock_widget->show(); }
  void on_Select_button() 
  {
    bool inside = ui_widget->Inside_check_box->isChecked();
    bool on_boundary = ui_widget->On_boundary_check_box->isChecked();
    bool outside = ui_widget->Outside_check_box->isChecked();

    if(!(inside || on_boundary || outside)) {
      print_message("Error: please check at least one parameter check box.");
      return;
    }

    // place all selected polyhedron and point items to vectors below
    std::vector<const Polyhedron*> polys;
    std::vector<Point_set*> point_sets;
    foreach(Scene_interface::Item_id id, scene->selectionIndices()) {
      Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(id));
      if(poly_item) { polys.push_back(poly_item->polyhedron()); }

      Scene_points_with_normal_item* point_item = qobject_cast<Scene_points_with_normal_item*>(scene->item(id));
      if(point_item) { point_sets.push_back(point_item->point_set()); }
    }
    
    // there should be at least one selected polyhedron and point item
    if(polys.empty())      { print_message("Error: there is no selected polyhedron item."); }
    if(point_sets.empty()) { print_message("Error: there is no selected point set item."); }
    if(polys.empty() || point_sets.empty()) { return; }

    // deselect all points
    for(std::vector<Point_set*>::iterator point_set_it = point_sets.begin(); 
      point_set_it != point_sets.end(); ++point_set_it) {
      (*point_set_it)->select((*point_set_it)->begin(), (*point_set_it)->end(), false);
    }

    // for each polyhedron, loop through all points of selected points sets
    for(std::vector<const Polyhedron*>::iterator poly_it = polys.begin(); 
      poly_it != polys.end(); ++poly_it)
    {
      CGAL::Point_inside_polyhedron_3<Polyhedron, Epic_kernel> query_functor(**poly_it);

      for(std::vector<Point_set*>::iterator point_set_it = point_sets.begin(); 
        point_set_it != point_sets.end(); ++point_set_it) 
      {
        Point_set* point_set = *point_set_it;
        for(Point_set::iterator point_it = point_set->begin(); point_it != point_set->end(); ++point_it)
        {
          CGAL::Bounded_side res = query_functor(point_it->position());

          if( (inside      && res == CGAL::ON_BOUNDED_SIDE) ||
              (on_boundary && res == CGAL::ON_BOUNDARY)     ||
              (outside     && res == CGAL::ON_UNBOUNDED_SIDE) )
          {
            point_set->select(&*point_it);
          }
        } // loop on points in point_set
      }// loop on selected point sets
    }// loop on selected polyhedrons

    // for repaint
    foreach(Scene_interface::Item_id id, scene->selectionIndices()) {
      Scene_points_with_normal_item* point_item = qobject_cast<Scene_points_with_normal_item*>(scene->item(id));
      if(point_item) { 
        scene->itemChanged(point_item);
      }
    }
  }

private:
  QMainWindow* mw;
  Scene_interface* scene;
  Messages_interface* messages;
  QAction* actionPointInsidePolyhedron;

  QDockWidget* dock_widget;
  Ui::Point_inside_polyhedron* ui_widget;

}; // end Polyhedron_demo_point_inside_polyhedron_plugin

Q_EXPORT_PLUGIN2(Polyhedron_demo_point_inside_polyhedron_plugin, Polyhedron_demo_point_inside_polyhedron_plugin)

#include "Polyhedron_demo_point_inside_polyhedron_plugin.moc"
