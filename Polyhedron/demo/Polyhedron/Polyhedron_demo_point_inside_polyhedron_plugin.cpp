#include <QtCore/qglobal.h>
#include "opengl_tools.h"

#include "Messages_interface.h"
#include "Scene_polyhedron_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_interface.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_type.h"

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

#include <boost/iterator/transform_iterator.hpp>
#include <boost/optional/optional.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;

class Polyhedron_demo_point_inside_polyhedron_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  bool applicable(QAction*) const 
  {
    bool poly_item_exists = false;
    bool point_item_exists = false;

    for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); 
      i < end && (!poly_item_exists || !point_item_exists); ++i)
    {
      poly_item_exists |= qobject_cast<Scene_polyhedron_item*>(scene->item(i)) != NULL;
      point_item_exists |= qobject_cast<Scene_points_with_normal_item*>(scene->item(i)) != NULL;
    }

    //return poly_item_exists && point_item_exists;
    return poly_item_exists || point_item_exists;
  }
  void print_message(QString message) { messages->information(message); }
  QList<QAction*> actions() const { return QList<QAction*>() << actionPointInsidePolyhedron; }

  using Polyhedron_demo_plugin_helper::init;
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* m)
  {
    mw = mainWindow;
    scene = scene_interface;
    messages = m;

    actionPointInsidePolyhedron = new QAction(tr("Point Inside Polyhedron"), mw);
    connect(actionPointInsidePolyhedron, SIGNAL(triggered()), this, SLOT(point_inside_polyhedron_action()));

    dock_widget = new QDockWidget("Point Inside Polyhedron", mw);
    dock_widget->setVisible(false);
    ui_widget.setupUi(dock_widget);

    add_dock_widget(dock_widget);

    connect(ui_widget.Select_button,  SIGNAL(clicked()), this, SLOT(on_Select_button())); 
    connect(ui_widget.Sample_random_points_from_bbox,  SIGNAL(clicked()), this, SLOT(on_Sample_random_points_from_bbox())); 
    
  }
private:
  // for transform iterator
  struct Get_ref {
    typedef const Polyhedron& result_type;
    result_type operator()(const Polyhedron* poly_ptr) const
    { return *poly_ptr; }
  };


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

    // place all selected polyhedron and point items to vectors below
    std::vector<const Polyhedron*> polys;
    typedef CGAL::Side_of_triangle_mesh<Polyhedron, Kernel> Point_inside;
    std::vector<Point_inside*> inside_testers;// to put all polyhedra to query object
        // it does not support copy-construction so let's use pointers
    std::vector<Point_set*> point_sets;
    Q_FOREACH(Scene_interface::Item_id id, scene->selectionIndices()) {
      Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(id));
      if (poly_item)
        inside_testers.push_back(new Point_inside(*(poly_item->polyhedron())));
      if(poly_item) { polys.push_back(poly_item->polyhedron()); }

      Scene_points_with_normal_item* point_item = qobject_cast<Scene_points_with_normal_item*>(scene->item(id));
      if(point_item) { point_sets.push_back(point_item->point_set()); }
    }

    // there should be at least one selected polyhedron and point item
    if(inside_testers.empty()) { print_message("Error: there is no selected polyhedron item(s)."); }
    if(point_sets.empty()) { print_message("Error: there is no selected point set item(s)."); }
    if(inside_testers.empty() || point_sets.empty()) { return; }

    // deselect all points
    for(std::vector<Point_set*>::iterator point_set_it = point_sets.begin(); 
      point_set_it != point_sets.end(); ++point_set_it) {
      (*point_set_it)->select((*point_set_it)->begin(), (*point_set_it)->end(), false);
    }

    CGAL::Timer timer; timer.start();

    print_message(
      QString("Constructing with %1 items is done in %2 sec.").arg(inside_testers.size()).arg(timer.time()));
    timer.reset();

    std::size_t nb_query = 0, nb_selected = 0;// for print message
    for(std::vector<Point_set*>::iterator point_set_it = point_sets.begin(); 
      point_set_it != point_sets.end(); ++point_set_it) 
    {
      Point_set* point_set = *point_set_it;
      for(Point_set::iterator point_it = point_set->begin();
          point_it != point_set->end();
          ++point_it, ++nb_query)
      {
        for (std::size_t i = 0; i < inside_testers.size(); ++i)
        {
        CGAL::Bounded_side res = (*inside_testers[i])(point_it->position());

        if( (inside      && res == CGAL::ON_BOUNDED_SIDE) ||
            (on_boundary && res == CGAL::ON_BOUNDARY)     ||
            (outside     && res == CGAL::ON_UNBOUNDED_SIDE) )
        {
          point_set->select(&*point_it); ++nb_selected;
          break;//loop on i
        }
        }
      } // loop on points in point_set
    }// loop on selected point sets

    print_message(QString("Querying with %1 points is done in %2 sec.").arg(nb_query).arg(timer.time()));
    print_message(QString("%1 points are selected. All Done!").arg(nb_selected));

    // delete testers
    for (std::size_t i = 0; i < inside_testers.size(); ++i)
      delete inside_testers[i];

    // for repaint
    Q_FOREACH(Scene_interface::Item_id id, scene->selectionIndices()) {
      Scene_points_with_normal_item* point_item = qobject_cast<Scene_points_with_normal_item*>(scene->item(id));
      if(point_item) { 
        point_item->invalidate_buffers();
        scene->itemChanged(point_item);
      }
    }
  }

  void on_Sample_random_points_from_bbox() {
    
    // calculate bbox of selected polyhedron items
    boost::optional<Scene_interface::Bbox> bbox
      = boost::make_optional(false, Scene_interface::Bbox());
    // Workaround a bug in g++-4.8.3:
    //   http://stackoverflow.com/a/21755207/1728537
    // Using boost::make_optional to copy-initialize 'bbox' hides the
    //   warning about '*bbox' not being initialized.
    // -- Laurent Rineau, 2014/10/30

    Q_FOREACH(Scene_interface::Item_id id, scene->selectionIndices()) {
      Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(id));
      if(poly_item) {
        if(!bbox) {
          bbox = poly_item->bbox();
        }
        else {
          *bbox = *bbox + poly_item->bbox();
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

    // sample random points and constuct item
    Scene_points_with_normal_item* point_item = new Scene_points_with_normal_item();
    point_item->setName(QString("sample-%1").arg(nb_points));
    CGAL::Random rg(1340818006);

    double grid_dx = bbox->xmax - bbox->xmin;
    double grid_dy = bbox->ymax - bbox->ymin;
    double grid_dz = bbox->zmax - bbox->zmin;

    for(int i=0; i < nb_points; i++){
      point_item->point_set()->push_back(
      Epic_kernel::Point_3(bbox->xmin + rg.get_double()* grid_dx, 
        bbox->ymin + rg.get_double()* grid_dy,
        bbox->zmin + rg.get_double()* grid_dz)
      );
    }

    scene->addItem(point_item);
    scene->itemChanged(point_item);
  }
private:
  Messages_interface* messages;
  QAction* actionPointInsidePolyhedron;

  QDockWidget* dock_widget;
  Ui::Point_inside_polyhedron ui_widget;

}; // end Polyhedron_demo_point_inside_polyhedron_plugin

#include "Polyhedron_demo_point_inside_polyhedron_plugin.moc"
