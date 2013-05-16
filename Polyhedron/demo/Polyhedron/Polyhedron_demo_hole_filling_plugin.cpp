#include <QtCore/qglobal.h>

#include "Messages_interface.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polylines_item.h"

#include "Polyhedron_demo_plugin_interface.h"
#include "ui_Hole_filling_widget.h"
#include "Polyhedron_type.h"

//#include <CGAL/Fill_hole.h>
#include <CGAL/Fill_hole_Polyhedron_3.h>

#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QDockWidget>

#include <vector>
#include <algorithm>
#include <boost/bimap.hpp>

class Polyhedron_demo_hole_filling_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  bool applicable() const { return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())); }
  void print_message(QString message) { messages->information(message);}
  QList<QAction*> actions() const { return QList<QAction*>() << actionHoleFilling; }
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* m);

public slots:
  void hole_filling_action();
  void on_Create_polyline_item_button();
  void on_Fill_and_update_button();

private:
  Scene_interface* scene;
  Messages_interface* messages;
  QAction* actionHoleFilling;

  QDockWidget* dock_widget;
  Ui::HoleFilling* ui_widget;

  typedef std::map<Scene_polylines_item*,
    std::pair<Scene_polyhedron_item*, Polyhedron::Halfedge_handle> > Polyline_data_map;
  Polyline_data_map polyline_data_map;
}; // end Polyhedron_demo_hole_filling_plugin

void Polyhedron_demo_hole_filling_plugin::init(QMainWindow* mw,
                                      Scene_interface* scene_interface,
                                      Messages_interface* m)
{
  scene = scene_interface;
  messages = m;
  actionHoleFilling = new QAction(tr("Hole Filling"), mw);
  connect(actionHoleFilling, SIGNAL(triggered()),
          this, SLOT(hole_filling_action()));

  dock_widget = new QDockWidget("Hole filling", mw);
  dock_widget->setVisible(false); // do not show at the beginning
  ui_widget = new Ui::HoleFilling();

  ui_widget->setupUi(dock_widget);
  mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

  connect(ui_widget->Create_polyline_item_button,  SIGNAL(clicked()), this, SLOT(on_Create_polyline_item_button()));  
  connect(ui_widget->Fill_and_update_button,  SIGNAL(clicked()), this, SLOT(on_Fill_and_update_button())); 
}

void Polyhedron_demo_hole_filling_plugin::hole_filling_action(){
  if(dock_widget != NULL) { 
    dock_widget->show(); 
  }
}

void Polyhedron_demo_hole_filling_plugin::on_Create_polyline_item_button() {
  typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
  typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;

  int item_id = scene->mainSelectionIndex();
  Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(item_id));
  if(!poly_item) {
    print_message("Error: there is no selected polyhedron item!");
    return;
  }

  Polyhedron& poly = *poly_item->polyhedron();  
  
  for(Halfedge_iterator it = poly.halfedges_begin(); it != poly.halfedges_end(); ++it)
  { it->id() = 0; }

  int counter = 0;
  for(Halfedge_iterator it = poly.halfedges_begin(); it != poly.halfedges_end(); ++it){
    if(it->is_border() && it->id() == 0){
      Scene_polylines_item* new_polylines_item = new Scene_polylines_item();
      polyline_data_map[new_polylines_item] = std::make_pair(poly_item, it);
      new_polylines_item->polylines.push_back(Scene_polylines_item::Polyline());
      Halfedge_around_facet_circulator hf_around_facet = it->facet_begin();
      do {
        CGAL_assertion(hf_around_facet->id() == 0);
        hf_around_facet->id() = 1;
        new_polylines_item->polylines.front().push_back(hf_around_facet->vertex()->point());
      } while(++hf_around_facet != it->facet_begin());
      new_polylines_item->polylines.front().push_back(hf_around_facet->vertex()->point());

      new_polylines_item->setName(QString("Border polyline %1").arg(counter++) );
      new_polylines_item->setColor(Qt::red);
      new_polylines_item->setRenderingMode(Wireframe);
      scene->addItem(new_polylines_item);
    }
  }
}

void Polyhedron_demo_hole_filling_plugin::on_Fill_and_update_button() {
  int item_id = scene->mainSelectionIndex();
  Scene_polylines_item* polyline_item = qobject_cast<Scene_polylines_item*>(scene->item(item_id));
  if(!polyline_item) {
    print_message("Error: there is no selected polyline item!");
    return;
  }
  Polyline_data_map::iterator it = polyline_data_map.find(polyline_item);
  if(it == polyline_data_map.end()) {
    print_message("Error: polyline should be associated with a polyhedron by this plugin!");
    return;
  }
  fill(*(it->second.first->polyhedron()), it->second.second);
  scene->itemChanged(it->second.first);
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_hole_filling_plugin, Polyhedron_demo_hole_filling_plugin)

#include "Polyhedron_demo_hole_filling_plugin.moc"
