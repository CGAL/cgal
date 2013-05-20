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
  void on_Fill_all_holes_button();
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
  dock_widget->setVisible(false);
  ui_widget = new Ui::HoleFilling();

  ui_widget->setupUi(dock_widget);
  mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

  connect(ui_widget->Create_polyline_item_button,  SIGNAL(clicked()), this, SLOT(on_Create_polyline_item_button()));  
  connect(ui_widget->Fill_and_update_button,  SIGNAL(clicked()), this, SLOT(on_Fill_and_update_button())); 
  connect(ui_widget->Fill_all_holes_button,  SIGNAL(clicked()), this, SLOT(on_Fill_all_holes_button()));
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

void Polyhedron_demo_hole_filling_plugin::on_Fill_all_holes_button() {
  typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
  typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;

  int item_id = scene->mainSelectionIndex();
  Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(item_id));
  if(!poly_item) {
    print_message("Error: there is no selected polyhedron item!");
    return;
  }

  // take parameters
  bool fair = ui_widget->Triangulate_refine_fair_radio_button->isChecked();
  bool refine = fair || ui_widget->Triangulate_refine_radio_button->isChecked();
  double alpha = ui_widget->Density_control_factor_spin_box->value();
  bool create_new = ui_widget->Create_new_polyhedron_check_box->checkState() == Qt::Checked;
  CGAL::Fairing_weight_type_tag weight_tag;
  if(ui_widget->Scale_dependent_weight_radio_button->isChecked()) {
    weight_tag = CGAL::SCALE_DEPENDENT_WEIGHTING;
  } else if(ui_widget->Uniform_weight_radio_button->isChecked()) {
    weight_tag = CGAL::UNIFORM_WEIGHTING;
  } else {
    weight_tag = CGAL::COTANGENT_WEIGHTING;
  }
  // create new polyhedron item if required
  Polyhedron* poly_pointer;
  Scene_polyhedron_item* new_item = 0;
  if(create_new) {
    new_item = new Scene_polyhedron_item(*poly_item->polyhedron());
    QString param_exp = fair ? "Refine + Fair" : refine ? "Refine" : "Triangulate";
    new_item->setName(tr("%1-%2-(alpha:%3)").arg(poly_item->name()).arg(param_exp).arg(alpha));
    poly_pointer = new_item->polyhedron();
  }
  else {
    poly_pointer = poly_item->polyhedron();
  }

  // TODO: check whether there is a better way to do this iteration
  bool any_changes = false;
  Polyhedron& poly = *poly_pointer;  
  for(Halfedge_iterator it = poly.halfedges_begin(); it != poly.halfedges_end(); ){
    if(it->is_border()){
      any_changes = true;
      CGAL::fill(poly, it, refine, alpha, fair, weight_tag);
      it = poly.halfedges_begin();
      continue;
    }
    ++it;
  }
  if(!any_changes) {
    print_message("There is no holes in selected polyhedron!");
    delete new_item; // if any, delete new_item
    return;
  }

  // update or add polyhedron item
  if(create_new) {
    Scene_interface::Item_id index = scene->addItem(new_item);
    scene->itemChanged(new_item);
    scene->setSelectedItem(index);
  }
  else {
    scene->itemChanged(poly_item);
  }

}

void Polyhedron_demo_hole_filling_plugin::on_Fill_and_update_button() {
  int item_id = scene->mainSelectionIndex();
  Scene_polylines_item* polyline_item = qobject_cast<Scene_polylines_item*>(scene->item(item_id));
  if(!polyline_item) {
    print_message("Error: there is no selected polyline item!");
    return;
  }

  // take parameters
  bool fair = ui_widget->Triangulate_refine_fair_radio_button->isChecked();
  bool refine = fair || ui_widget->Triangulate_refine_radio_button->isChecked();
  double alpha = ui_widget->Density_control_factor_spin_box->value();
  bool create_new = ui_widget->Create_new_polyhedron_check_box->checkState() == Qt::Checked;

  Polyline_data_map::iterator it = polyline_data_map.find(polyline_item);
  if(it == polyline_data_map.end()) {
    print_message("Error: polyline should be associated with a polyhedron by this plugin!");
    return;
  }

  Scene_polyhedron_item* poly_item = it->second.first;
  //// create new polyhedron item if required
  //Scene_polyhedron_item* new_item = 0;
  //if(create_new) {
  //  new_item = new Scene_polyhedron_item(*poly_item->polyhedron());
  //  QString param_exp = fair ? "Refine + Fair" : refine ? "Refine" : "Triangulate";
  //  new_item->setName(tr("%1-Filled-%2-(alpha:%3)").arg(poly_item->name()).arg(param_exp).arg(alpha));
  //}

  CGAL::fill(*poly_item->polyhedron(), it->second.second, refine, alpha, fair);
  scene->itemChanged(poly_item);
  //if(create_new) {
  //  using std::swap;
  //  swap(*poly_item->polyhedron(), *new_item->polyhedron());
  //  Scene_interface::Item_id index = scene->addItem(new_item);
  //  scene->itemChanged(new_item);
  //  scene->setSelectedItem(index);
  //}  
  //else {
  //  scene->itemChanged(poly_item);
  //}
  // remove polyline
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_hole_filling_plugin, Polyhedron_demo_hole_filling_plugin)

#include "Polyhedron_demo_hole_filling_plugin.moc"
