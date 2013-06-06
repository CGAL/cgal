#undef NDEBUG
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

#include <boost/function_output_iterator.hpp>
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

  Scene_polylines_collection(Scene_polyhedron_item* poly_item, QMainWindow* mainWindow)
    : poly_item(poly_item) 
  {
    get_holes();
    active_hole = polyline_data_list.end();

    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    viewer->installEventFilter(this);

    mainWindow->installEventFilter(this);

    connect(poly_item, SIGNAL(itemChanged()), this, SLOT(poly_item_changed())); 
  }
  ~Scene_polylines_collection() {
    clear();
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
    
    for(Polyline_data_list::const_iterator it = polyline_data_list.begin(); it != polyline_data_list.end(); ++it) {
      if(it == active_hole) { ::glLineWidth(7.f); }
      else                  { ::glLineWidth(3.f); }

      if(selected_holes.find(it) != selected_holes.end()) 
      { it->polyline->setRbgColor(255, 0, 0); }
      else 
      { it->polyline->setRbgColor(0, 0, 255); }

      it->polyline->draw_edges();
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
  // find holes in polyhedron and construct a internal polyline for each
  void get_holes() {
    typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
    typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;

    clear();

    Polyhedron& poly = *poly_item->polyhedron();
    for(Halfedge_iterator it = poly.halfedges_begin(); it != poly.halfedges_end(); ++it)
    { it->id() = 0; }

    int counter = 0;
    for(Halfedge_iterator it = poly.halfedges_begin(); it != poly.halfedges_end(); ++it){
      if(it->is_border() && it->id() == 0){
        polyline_data_list.push_back(Polyline_data());
        Polyline_data& polyline_data = polyline_data_list.back();
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

public slots:
  void poly_item_changed() {
    get_holes();
    emit itemChanged();
  }

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

  template<class SceneType>
  SceneType* get_selected() {
    int item_id = scene->mainSelectionIndex();
    SceneType* scene_item = qobject_cast<SceneType*>(scene->item(item_id));
    if(!scene_item) {
      // no selected SceneType, if there is only one in list use it, otherwise error
      int counter = 0;
      for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end && counter < 2; ++i) {
        if(SceneType* tmp = qobject_cast<SceneType*>(scene->item(i))) { 
          scene_item = tmp;
          counter++; 
        }
      }
      if(counter != 1) { return NULL; }
    }
    return scene_item;
  }

  Scene_polylines_collection* get_polylines_collection(Scene_polyhedron_item* poly_item) {
    // did not use a map to assoc Scene_polyhedron_item with Scene_polylines_collection to prevent crowded code
    for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i) {
      Scene_polylines_collection* polylines_collection = qobject_cast<Scene_polylines_collection*>(scene->item(i));
      if(polylines_collection && polylines_collection->poly_item == poly_item) 
      { return polylines_collection; }
    }
    return NULL;
  }

public slots:
  void hole_filling_action() { dock_widget->show(); }
  void on_Fill_all_holes_button();
  void on_Visualize_holes_button();
  void on_Fill_selected_holes_button();
  void item_about_to_be_destroyed(Scene_item*);
  void dock_widget_visibility_changed(bool visible);
  void item_changed_polylines_collection();

private:
  struct Nop_functor {
    template<class T>
    void operator()(const T& t) const {}
  };
  typedef boost::function_output_iterator<Nop_functor> Nop_out;

  QMainWindow* mw;
  Scene_interface* scene;
  Messages_interface* messages;
  QAction* actionHoleFilling;

  QDockWidget* dock_widget;
  Ui::HoleFilling* ui_widget;

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

void Polyhedron_demo_hole_filling_plugin::item_about_to_be_destroyed(Scene_item* scene_item) {
  Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene_item);
  if(poly_item) {
    scene->erase( scene->item_id( get_polylines_collection(poly_item) ) );
  }
}
// removes Scene_polylines_collection items on visibility = false
void Polyhedron_demo_hole_filling_plugin::dock_widget_visibility_changed(bool visible)
{
  if(visible) { return; }
  // remove all Scene_polylines_collection items
  for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries();
    i < end; ++i)
  {
    Scene_polylines_collection* polylines_collection = qobject_cast<Scene_polylines_collection*>(scene->item(i));
    if(polylines_collection) {
      scene->erase( scene->item_id(polylines_collection) );
    }
  }
}
// creates a Scene_polylines_collection and associate it with active Scene_polyhedron_item
void Polyhedron_demo_hole_filling_plugin::on_Visualize_holes_button() {
  typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
  typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;

  Scene_polyhedron_item* poly_item = get_selected<Scene_polyhedron_item>();
  if(!poly_item) {
    print_message("Error: please select a polyhedron item from Geometric Objects list!");
    return;
  }

  if(get_polylines_collection(poly_item)) {
    print_message("Error: selected polyhedron item already has an associated hole item!");
    return;
  }

  Scene_polylines_collection* polylines_collection = new Scene_polylines_collection(poly_item, mw);
  connect(polylines_collection, SIGNAL(itemChanged()), this, SLOT(item_changed_polylines_collection()));

  if(polylines_collection->polyline_data_list.empty()) {
    print_message("There is no holes in selected polyhedron!");
    delete polylines_collection;
    return;
  }
  else {
    poly_item->setFlatMode(); // for better visualization
    int item_id = scene->addItem(polylines_collection);
    scene->setSelectedItem(item_id);
    polylines_collection->setName(tr("%1-hole visualizer").arg(poly_item->name()));
    // poly_item->setColor(QColor(poly_item->color().red(), poly_item->color().green(), poly_item->color().blue(), 100));
  }
}
// fills selected holes on active Scene_polylines_collection
void Polyhedron_demo_hole_filling_plugin::on_Fill_selected_holes_button() {
  Scene_polylines_collection* polyline_item = get_selected<Scene_polylines_collection>();
  if(!polyline_item) {
    print_message("Error: please select a hole visualizer from Geometric Objects list!");
    return;
  }

  bool create_new = ui_widget->Create_new_polyhedron_check_box->checkState() == Qt::Checked;
  
  for(Scene_polylines_collection::Selected_holes_set::iterator it = polyline_item->selected_holes.begin();
    it != polyline_item->selected_holes.end(); ++it) {
      fill(*(polyline_item->poly_item->polyhedron()), (*it)->halfedge);
  }
  scene->itemChanged(polyline_item->poly_item);
  polyline_item->poly_item_changed();
};

// fills all holes and removes associated Scene_polylines_collection if any
void Polyhedron_demo_hole_filling_plugin::on_Fill_all_holes_button() {
  typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
  typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;

  Scene_polyhedron_item* poly_item = get_selected<Scene_polyhedron_item>();
  if(!poly_item) {
    print_message("Error: please select a polyhedron item from Geometric Objects list!");
    return;
  }

  bool create_new = ui_widget->Create_new_polyhedron_check_box->checkState() == Qt::Checked;
  double alpha = ui_widget->Density_control_factor_spin_box->value();
  int action_index = ui_widget->action_combo_box->currentIndex();

  // create new polyhedron item if required
  Polyhedron* poly_pointer;
  Scene_polyhedron_item* new_item = 0;
  if(create_new) {
    new_item = new Scene_polyhedron_item(*poly_item->polyhedron());
    QString param_exp = action_index == 2 ? "Refine + Fair" : action_index == 1 ? "Refine" : "Triangulate";
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
    Scene_polylines_collection* polylines_collection = get_polylines_collection(poly_item);
    if(polylines_collection) { polylines_collection->poly_item_changed();}
  }
}
// helper function for filling holes
void Polyhedron_demo_hole_filling_plugin::fill
  (Polyhedron& poly, Polyhedron::Halfedge_handle it) {

  int action_index = ui_widget->action_combo_box->currentIndex();
  double alpha = ui_widget->Density_control_factor_spin_box->value();

  if(action_index == 0) {
    CGAL::triangulate_hole(poly, it, Nop_out());
  }
  else if(action_index == 1) {
    CGAL::triangulate_and_refine_hole(poly, it, Nop_out(), Nop_out(), alpha);
  }
  else {
    int weight_index = ui_widget->weight_combo_box->currentIndex();

    if(weight_index == 0)
      CGAL::triangulate_refine_and_fair_hole(poly, it, Nop_out(), Nop_out(), alpha, 
        CGAL::Fairing_uniform_weight<Polyhedron>());
    if(weight_index == 1)
      CGAL::triangulate_refine_and_fair_hole(poly, it, Nop_out(), Nop_out(), alpha,
       CGAL::Fairing_cotangent_weight<Polyhedron>());
  }
}

void Polyhedron_demo_hole_filling_plugin::item_changed_polylines_collection() {
  Scene_polylines_collection* polylines_collection = qobject_cast<Scene_polylines_collection*>(this->sender());
  if(polylines_collection && polylines_collection->polyline_data_list.empty()) {
    scene->erase( scene->item_id(polylines_collection));
  }
}
Q_EXPORT_PLUGIN2(Polyhedron_demo_hole_filling_plugin, Polyhedron_demo_hole_filling_plugin)

#include "Polyhedron_demo_hole_filling_plugin.moc"
