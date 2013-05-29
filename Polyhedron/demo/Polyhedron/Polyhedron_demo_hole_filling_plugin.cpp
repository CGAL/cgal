#include <QtCore/qglobal.h>

#include "Messages_interface.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polylines_item.h"
#include "Scene.h"

#include "Polyhedron_demo_plugin_interface.h"
#include "ui_Hole_filling_widget.h"
#include "Polyhedron_type.h"

#include <CGAL/Fill_hole_Polyhedron_3.h>

#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QDockWidget>

#include <QEvent>
#include <QKeyEvent>
#include <QMouseEvent>

#include <vector>
#include <algorithm>

#include <QGLViewer/qglviewer.h>
#include <CGAL/gl_render.h>

// Class for visualizing holes in a polyhedron
// provides mouse selection functionality
class Q_DECL_EXPORT Scene_polylines_collection : public Scene_item
{
  Q_OBJECT
public:
  // structs
  struct Polyline_data {
    Scene_polylines_item* polyline;
    Polyhedron::Halfedge_handle halfedge;
    qglviewer::Vec position;
  };
  struct Mouse_keyboard_state
  {
    bool ctrl_pressing, left_button_pressing;
    Mouse_keyboard_state() : ctrl_pressing(false), left_button_pressing(false) { }
  };
public: typedef std::list<Polyline_data> Polyline_data_list;
private:
  struct List_iterator_comparator {
    bool operator()(Polyline_data_list::const_iterator it_1, Polyline_data_list::const_iterator it_2) const 
    { return (&*it_1) < (&*it_2); }
  };
public:
  typedef std::set<Polyline_data_list::const_iterator, List_iterator_comparator> Selected_holes_set;

  Scene_polylines_collection(Scene_polyhedron_item* poly_item)
    : poly_item(poly_item) {
    active_hole = polyline_data_list.end();

    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    viewer->installEventFilter(this);
  }
  ~Scene_polylines_collection() {
    for(Polyline_data_list::const_iterator it = polyline_data_list.begin(); it != polyline_data_list.end(); ++it) {
      delete it->polyline;
    }
  }
  bool isFinite() const { return true; }
  bool isEmpty() const { return polyline_data_list.empty(); }
  Bbox bbox() const {
    if(polyline_data_list.empty()) { return Bbox(); }
    Bbox bbox = polyline_data_list.begin()->polyline->bbox();
    for(Polyline_data_list::const_iterator it = polyline_data_list.begin(); it != polyline_data_list.end(); ++it) {
      bbox = bbox + it->polyline->bbox();
    }
    return bbox;
  }
  Scene_polylines_collection* clone() const {
    return 0;
  }
  QString toolTip() const {
    return tr("%1 with %2 holes").arg(name()).arg(polyline_data_list.size());
  }

  bool supportsRenderingMode(RenderingMode m) const {
    return (m == Wireframe);
  }
  void draw() const {}
  void draw_edges() const {
    ::glLineWidth(3.f);
    for(Polyline_data_list::const_iterator it = polyline_data_list.begin(); it != polyline_data_list.end(); ++it) {
      if(selected_holes.find(it) != selected_holes.end()) 
      { it->polyline->setRbgColor(0, 255, 0); }
      else if(it == active_hole) 
      { it->polyline->setRbgColor(255, 0, 0); }
      else 
      { it->polyline->setRbgColor(0, 0, 255); }
      it->polyline->draw_edges();
    }
  }

  // find holes in polyhedron and construct a internal polyline for each
  void get_holes(Polyhedron& poly) {
    typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
    typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;

    clear();

    for(Halfedge_iterator it = poly.halfedges_begin(); it != poly.halfedges_end(); ++it)
    { it->id() = 0; }

    int counter = 0;
    for(Halfedge_iterator it = poly.halfedges_begin(); it != poly.halfedges_end(); ++it){
      if(it->is_border() && it->id() == 0){
        polyline_data_list.push_back(Polyline_data());
        Scene_polylines_collection::Polyline_data& polyline_data = polyline_data_list.back();
        polyline_data.polyline = new Scene_polylines_item();
        polyline_data.polyline->polylines.push_back(Scene_polylines_item::Polyline());
        polyline_data.halfedge = it;

        qglviewer::Vec center;
        int counter = 0;
        Halfedge_around_facet_circulator hf_around_facet = it->facet_begin();
        do {
          CGAL_assertion(hf_around_facet->id() == 0);
          hf_around_facet->id() = 1;
          const Polyhedron::Traits::Point_3& p = hf_around_facet->vertex()->point();          
          polyline_data.polyline->polylines.front().push_back(p);
          center += qglviewer::Vec(p.x(), p.y(), p.z());
          ++counter;
        } while(++hf_around_facet != it->facet_begin());
        polyline_data.polyline->polylines.front().push_back(hf_around_facet->vertex()->point());
        polyline_data.position = center / counter;
      }
    }
  }
  // filter events for selecting / activating holes with mouse input
  bool eventFilter(QObject* /*target*/, QEvent *event)
  {
    // This filter is both filtering events from 'viewer' and 'main window'
    Mouse_keyboard_state old_state = state;
    // key events
    if(event->type() == QEvent::KeyPress || event->type() == QEvent::KeyRelease)  {
      QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
      Qt::KeyboardModifiers modifiers = keyEvent->modifiers();

      state.ctrl_pressing = modifiers.testFlag(Qt::ControlModifier);
    }
    // mouse events
    if(event->type() == QEvent::MouseButtonPress || event->type() == QEvent::MouseButtonRelease) {
      QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
      if(mouse_event->button() == Qt::LeftButton) {
        state.left_button_pressing = event->type() == QEvent::MouseButtonPress;
      }   
    }

    if(!visible()) { return false; } // if not visible just update event state but don't do any action

    // activate closest hole
    if(event->type() == QEvent::HoverMove) 
    {
      QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
      const QPoint& p = viewer->mapFromGlobal(QCursor::pos());
      bool need_repaint = activate_closest_hole(p.x(), p.y());
      if(need_repaint) { emit itemChanged(); }
    }

    // select closest hole
    bool left_clicked_now = state.left_button_pressing && !old_state.left_button_pressing;
    if(left_clicked_now && state.ctrl_pressing) {
      Selected_holes_set::iterator active_it = selected_holes.find(active_hole);
      if(active_it == selected_holes.end()) {
        selected_holes.insert(active_hole);
      }
      else { selected_holes.erase(active_it); }
      emit itemChanged();
    }
    return false;
  }

private:
  // finds closest polyline from polyline_data_list and makes it active_hole
  bool activate_closest_hole(int x, int y) {
    if(polyline_data_list.empty()) { return false; }

    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    qglviewer::Camera* camera = viewer->camera();

    Polyline_data_list::const_iterator min_it = polyline_data_list.begin();    
    const qglviewer::Vec& pos_it = camera->projectedCoordinatesOf(min_it->position);
    float min_dist = std::pow(pos_it.x - x, 2) + std::pow(pos_it.y - y, 2);

    for(Polyline_data_list::const_iterator it = polyline_data_list.begin(); it != polyline_data_list.end(); ++it)
    {
      const qglviewer::Vec& pos_it = camera->projectedCoordinatesOf(it->position);
      float dist = std::pow(pos_it.x - x, 2) + std::pow(pos_it.y - y, 2);
      if(dist < min_dist) {
        min_dist = dist;
        min_it = it;
      }
    }

    if(min_it == active_hole) { 
      return false;
    }
    active_hole = min_it;
    return true;
  }
  // clears internal data except poly_item
  void clear() {
    for(Polyline_data_list::const_iterator it = polyline_data_list.begin(); it != polyline_data_list.end(); ++it) {
      delete it->polyline;
    }
    polyline_data_list.clear();
    selected_holes.clear();
    active_hole = polyline_data_list.end();
  }
  
  Polyline_data_list::const_iterator active_hole;
  Mouse_keyboard_state state;
public:
  Selected_holes_set selected_holes;
  Scene_polyhedron_item* poly_item;
  Polyline_data_list polyline_data_list;
}; // end class Scene_edges_item
///////////////////////////////////////////////////////////////////////////////////////////////////

class Polyhedron_demo_hole_filling_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  bool applicable() const { return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())); }
  void print_message(QString message) { messages->information(message); }
  QList<QAction*> actions() const { return QList<QAction*>() << actionHoleFilling; }
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* m);

public slots:
  void hole_filling_action() { dock_widget->show(); }
  void on_Fill_all_holes_button();
  void on_Visualize_holes_button();
  void on_Fill_selected_holes_button();
  void item_about_to_be_destroyed(Scene_item*);
  void dock_widget_visibility_changed(bool visible);
private:
  QMainWindow* mw;
  Scene_interface* scene;
  Messages_interface* messages;
  QAction* actionHoleFilling;

  QDockWidget* dock_widget;
  Ui::HoleFilling* ui_widget;

  typedef std::map<Scene_polyhedron_item*, Scene_polylines_collection*> Polyhedron_item_hole_map;
  Polyhedron_item_hole_map polyhedron_item_hole_map;

  void fill(Polyhedron& polyhedron, Polyhedron::Halfedge_handle halfedge);
}; // end Polyhedron_demo_hole_filling_plugin

void Polyhedron_demo_hole_filling_plugin::init(QMainWindow* mainWindow,
                                      Scene_interface* scene_interface,
                                      Messages_interface* m)
{
  mw = mainWindow;
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

  connect(dock_widget, SIGNAL(visibilityChanged(bool)), this, SLOT(dock_widget_visibility_changed(bool)) );
  connect(ui_widget->Visualize_holes_button,  SIGNAL(clicked()), this, SLOT(on_Visualize_holes_button()));  
  connect(ui_widget->Fill_selected_holes_button,  SIGNAL(clicked()), this, SLOT(on_Fill_selected_holes_button())); 
  connect(ui_widget->Fill_all_holes_button,  SIGNAL(clicked()), this, SLOT(on_Fill_all_holes_button()));

  Scene* scene_casted = dynamic_cast<Scene*>(scene_interface);
  if(scene_casted) 
  { connect(scene_casted, SIGNAL(itemAboutToBeDestroyed(Scene_item*)), this, SLOT(item_about_to_be_destroyed(Scene_item*))); }
}

// removes items from polyhedron_item_hole_map
// if item is Scene_polyhedron_item then takes its Scene_polylines_collection and erases it from scene 
// (this will cause another call to this handler)
// else (item is Scene_polylines_collection) removes it from the map
void Polyhedron_demo_hole_filling_plugin::item_about_to_be_destroyed(Scene_item* scene_item) {
  Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene_item);
  if(poly_item) {
    Polyhedron_item_hole_map::iterator hole_it = polyhedron_item_hole_map.find(poly_item);
    if(hole_it != polyhedron_item_hole_map.end()) {
      scene->erase( scene->item_id(hole_it->second) );
    }
  }

  Scene_polylines_collection* polyline_item = qobject_cast<Scene_polylines_collection*>(scene_item);
  if(polyline_item) {
    polyhedron_item_hole_map.erase(polyline_item->poly_item);
  }
}
// removes Scene_polylines_collection items on visibility = false
void Polyhedron_demo_hole_filling_plugin::dock_widget_visibility_changed(bool visible)
{
  if(visible) { return; }

  for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries();
    i < end; ++i)
  {
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(i));
    if(poly_item) {
      Polyhedron_item_hole_map::iterator hole_it = polyhedron_item_hole_map.find(poly_item);
      if(hole_it != polyhedron_item_hole_map.end()) {
        scene->erase( scene->item_id(hole_it->second) );
      }
    }
  }
}
// creates a Scene_polylines_collection and associate it with active Scene_polyhedron_item
void Polyhedron_demo_hole_filling_plugin::on_Visualize_holes_button() {
  typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
  typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;

  int item_id = scene->mainSelectionIndex();
  Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(item_id));
  if(!poly_item) {
    print_message("Error: there is no selected polyhedron item!");
    return;
  }
  if(polyhedron_item_hole_map.find(poly_item) != polyhedron_item_hole_map.end()) {
    print_message("Error: selected polyhedron item already has an associated hole item!");
    return;
  }
  Scene_polylines_collection* polylines_collection = new Scene_polylines_collection(poly_item);
  polylines_collection->get_holes(*poly_item->polyhedron());

  if(polylines_collection->polyline_data_list.empty()) {
    print_message("There is no holes in selected polyhedron!");
    delete polylines_collection;
    return;
  }
  else {
    poly_item->setFlatMode(); // for better visualization
    polyhedron_item_hole_map[poly_item] = polylines_collection;
    int item_id = scene->addItem(polylines_collection);
    scene->setSelectedItem(item_id);
    mw->installEventFilter(polylines_collection);
    polylines_collection->setName(tr("%1-hole visualizer").arg(poly_item->name()));
    // poly_item->setColor(QColor(poly_item->color().red(), poly_item->color().green(), poly_item->color().blue(), 100));
  }
}
// fills selected holes on active Scene_polylines_collection
void Polyhedron_demo_hole_filling_plugin::on_Fill_selected_holes_button() {
  int item_id = scene->mainSelectionIndex();
  Scene_polylines_collection* polyline_item = qobject_cast<Scene_polylines_collection*>(scene->item(item_id));
  if(!polyline_item) {
    print_message("Error: there is no selected holes item!");
    return;
  }

  bool create_new = ui_widget->Create_new_polyhedron_check_box->checkState() == Qt::Checked;
  
  for(Scene_polylines_collection::Selected_holes_set::iterator it = polyline_item->selected_holes.begin();
    it != polyline_item->selected_holes.end(); ++it) {
      fill(*(polyline_item->poly_item->polyhedron()), (*it)->halfedge);
  }
  polyline_item->get_holes(*polyline_item->poly_item->polyhedron());

  scene->itemChanged(polyline_item->poly_item);
  if(polyline_item->polyline_data_list.empty()) {
    scene->erase( scene->item_id(polyline_item) );
  }
  else { scene->itemChanged(polyline_item); }
};
// fills all holes and removes associated Scene_polylines_collection if any
void Polyhedron_demo_hole_filling_plugin::on_Fill_all_holes_button() {
  typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
  typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;

  int item_id = scene->mainSelectionIndex();
  Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(item_id));
  if(!poly_item) {
    print_message("Error: there is no selected polyhedron item!");
    return;
  }

  bool create_new = ui_widget->Create_new_polyhedron_check_box->checkState() == Qt::Checked;
  bool fair = ui_widget->Triangulate_refine_fair_radio_button->isChecked();
  bool refine = fair || ui_widget->Triangulate_refine_radio_button->isChecked();
  double alpha = ui_widget->Density_control_factor_spin_box->value();

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
      fill(poly, it);
      it = poly.halfedges_begin();
    }
    else { ++it; }
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
  // remove assoc Scene_polylines_collection item if any
  if(!create_new) {
    Polyhedron_item_hole_map::iterator hole_it = polyhedron_item_hole_map.find(poly_item);
    if(hole_it != polyhedron_item_hole_map.end()) {
      scene->erase( scene->item_id(hole_it->second) );
    }
  }
}
// helper function for filling holes
void Polyhedron_demo_hole_filling_plugin::fill
  (Polyhedron& poly, Polyhedron::Halfedge_handle it) {

  bool fair = ui_widget->Triangulate_refine_fair_radio_button->isChecked();
  bool refine = fair || ui_widget->Triangulate_refine_radio_button->isChecked();
  double alpha = ui_widget->Density_control_factor_spin_box->value();

  if(!fair && !refine) {
    CGAL::triangulate_hole(poly, it);
  }
  else if(!fair) {
    CGAL::triangulate_and_refine_hole(poly, it, alpha);
  }
  else {
    if(ui_widget->Scale_dependent_weight_radio_button->isChecked())
      CGAL::triangulate_refine_and_fair_hole(poly, it, alpha,
        CGAL::Fairing_scale_dependent_weight<Polyhedron>());
    if(ui_widget->Uniform_weight_radio_button->isChecked())
      CGAL::triangulate_refine_and_fair_hole(poly, it, alpha, 
        CGAL::Fairing_uniform_weight<Polyhedron>());
    else
      CGAL::triangulate_refine_and_fair_hole(poly, it, alpha,
        CGAL::Fairing_cotangent_weight<Polyhedron>());
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_hole_filling_plugin, Polyhedron_demo_hole_filling_plugin)

#include "Polyhedron_demo_hole_filling_plugin.moc"
