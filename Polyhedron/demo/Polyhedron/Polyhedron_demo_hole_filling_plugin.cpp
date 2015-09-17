#undef NDEBUG
#define DEBUG_TRACE
#include <QtCore/qglobal.h>

#include "Messages_interface.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polylines_item.h"
#include "Scene.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "ui_Hole_filling_widget.h"
#include "Polyhedron_type.h"

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Timer.h>
#include <CGAL/iterator.h>

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

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include "Kernel_type.h"

#include <boost/function_output_iterator.hpp>
#include <QMap>

// Class for visualizing holes in a polyhedron
// provides mouse selection functionality
class Q_DECL_EXPORT Scene_hole_visualizer : public Scene_item
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

  Scene_hole_visualizer(Scene_polyhedron_item* poly_item, QMainWindow* mainWindow)
    : poly_item(poly_item), block_poly_item_changed(false)
  {
    get_holes();
    active_hole = polyline_data_list.end();

    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    viewer->installEventFilter(this);
    mainWindow->installEventFilter(this);

    connect(poly_item, SIGNAL(item_is_about_to_be_changed()), this, SLOT(poly_item_changed())); 
  }
  ~Scene_hole_visualizer() {
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
  Scene_hole_visualizer* clone() const {
    return 0;
  }
  QString toolTip() const {
    return tr("%1 with %2 holes").arg(name()).arg(polyline_data_list.size());
  }

  bool supportsRenderingMode(RenderingMode m) const {
    return (m == Wireframe);
  }
  void draw() const {}
  void draw_edges(Viewer_interface* viewer) const {
    
    for(Polyline_data_list::const_iterator it = polyline_data_list.begin(); it != polyline_data_list.end(); ++it) {
      if(it == active_hole) { viewer->glLineWidth(7.f); }
      else                  { viewer->glLineWidth(3.f); }

      if(selected_holes.find(it) != selected_holes.end()) 
      { it->polyline->setRbgColor(255, 0, 0); }
      else 
      { it->polyline->setRbgColor(0, 0, 255); }

      it->polyline->draw_edges(viewer);
    }
  }

  void select_deselect_all(bool select) {
    if(select) {
      for(Polyline_data_list::iterator it = polyline_data_list.begin(); it != polyline_data_list.end(); ++it)
      { selected_holes.insert(it); }
    }
    else {
      selected_holes.clear();
    }
    Q_EMIT itemChanged();
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
      if(need_repaint) { Q_EMIT itemChanged(); }
    }

    // select closest hole
    bool left_clicked_now = state.left_button_pressing && !old_state.left_button_pressing;
    if(left_clicked_now && state.ctrl_pressing) {
      Selected_holes_set::iterator active_it = selected_holes.find(active_hole);
      if(active_it == selected_holes.end()) {
        selected_holes.insert(active_hole);
      }
      else { selected_holes.erase(active_it); }
      Q_EMIT itemChanged();
    }
    return false;
  }

private:
  // find holes in polyhedron and construct a internal polyline for each
  void get_holes() {
    typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
    typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
    // save selected hole positions to keep selected holes selected
    // we just use center position of holes for identification which might not work good for advanced cases...
    std::vector<qglviewer::Vec> selected_hole_positions;
    for(Selected_holes_set::const_iterator it = selected_holes.begin(); it != selected_holes.end(); ++it) {
      selected_hole_positions.push_back((*it)->position);
    }

    clear();

    Polyhedron& poly = *poly_item->polyhedron();
    for(Halfedge_iterator it = poly.halfedges_begin(); it != poly.halfedges_end(); ++it)
    { it->id() = 0; }

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
    //keep previous selected holes selected
    for(Polyline_data_list::const_iterator it = polyline_data_list.begin(); it != polyline_data_list.end(); ++it) {
      if(std::find(selected_hole_positions.begin(), selected_hole_positions.end(), it->position) != selected_hole_positions.end()) {
        selected_holes.insert(it);
      }
    }
  }
  // finds closest polyline from polyline_data_list and makes it active_hole
  bool activate_closest_hole(int x, int y) {
    typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
    if(polyline_data_list.empty()) { return false; }

    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    qglviewer::Camera* camera = viewer->camera();

    Polyline_data_list::const_iterator min_it;
    double min_dist = (std::numeric_limits<double>::max)();
    Kernel::Point_2 xy(x,y);
    for(Polyline_data_list::const_iterator it = polyline_data_list.begin(); it != polyline_data_list.end(); ++it)
    {
#if 0
      /* use center of polyline to measure distance - performance wise */
      const qglviewer::Vec& pos_it = camera->projectedCoordinatesOf(it->position);
      float dist = std::pow(pos_it.x - x, 2) + std::pow(pos_it.y - y, 2);
      if(dist < min_dist) {
        min_dist = dist;
        min_it = it;
      }
#else
      /* use polyline points to measure distance - might hurt performance for large holes */
      Halfedge_around_facet_circulator hf_around_facet = it->halfedge->facet_begin();
      do {
        
        const Polyhedron::Traits::Point_3& p_1 = hf_around_facet->vertex()->point();
        const qglviewer::Vec& pos_it_1 = camera->projectedCoordinatesOf(qglviewer::Vec(p_1.x(), p_1.y(), p_1.z()));
        const Polyhedron::Traits::Point_3& p_2 = hf_around_facet->opposite()->vertex()->point();
        const qglviewer::Vec& pos_it_2 = camera->projectedCoordinatesOf(qglviewer::Vec(p_2.x(), p_2.y(), p_2.z()));
        Kernel::Segment_2 s(Kernel::Point_2(pos_it_1.x, pos_it_1.y), Kernel::Point_2(pos_it_2.x, pos_it_2.y));

        double dist = CGAL::squared_distance(s, xy);
        if(dist < min_dist) {
          min_dist = dist;
          min_it = it;
        }
      } while(++hf_around_facet != it->halfedge->facet_begin());
#endif
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
  bool block_poly_item_changed;

public Q_SLOTS:
  void poly_item_changed() {
    if(block_poly_item_changed) { return; }
    get_holes();
    Q_EMIT itemChanged();
  }

}; // end class Scene_hole_visualizer
///////////////////////////////////////////////////////////////////////////////////////////////////

class Polyhedron_demo_hole_filling_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
  bool applicable(QAction*) const { return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())); }
  void print_message(QString message) { messages->information(message); }
  QList<QAction*> actions() const { return QList<QAction*>() << actionHoleFilling; }

  using Polyhedron_demo_plugin_helper::init;
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* m);

  Scene_hole_visualizer* get_hole_visualizer(Scene_polyhedron_item* poly_item) {
      return visualizers[poly_item];
  }

public Q_SLOTS:
  void hole_filling_action() { 
    dock_widget->show();
    dock_widget->raise();
    if(scene->numberOfEntries() < 2) { on_Visualize_holes_button(); }
  }
  void on_Select_all_holes_button();
  void on_Deselect_all_holes_button();
  void on_Visualize_holes_button();
  void on_Fill_selected_holes_button();
  void on_Create_polyline_items_button();
  void on_Accept_button();
  void on_Reject_button();
  void item_about_to_be_destroyed(Scene_item*);
  void hole_visualizer_changed();
  void dock_widget_closed();
  void on_Select_small_holes_button();
protected:
  bool eventFilter(QObject *, QEvent *event) {
    if(event->type() == QEvent::Close) {
      dock_widget_closed();
    }
    return false;
  }

  void change_poly_item_by_blocking(Scene_polyhedron_item* poly_item, Scene_hole_visualizer* collection) {
    if(collection) collection->block_poly_item_changed = true;
    poly_item->invalidate_buffers();
    scene->itemChanged(poly_item);
    if(collection) collection->block_poly_item_changed = false;
  }
private:
  Messages_interface* messages;
  QAction* actionHoleFilling;

  QDockWidget* dock_widget;
  Ui::HoleFilling ui_widget;

  //Maintains a reference between all the visualizers and their poly_item
  // to ease the management of the visualizers
  QMap<Scene_polyhedron_item*, Scene_hole_visualizer*> visualizers;
  // hold created facet for accept reject functionality
  std::vector<Polyhedron::Facet_handle> new_facets; 
  Scene_polyhedron_item* last_active_item; // always keep it NULL while not active-reject state
 
  bool fill(Polyhedron& polyhedron, Polyhedron::Halfedge_handle halfedge);
  bool self_intersecting(Polyhedron& polyhedron);
  void accept_reject_toggle(bool activate_accept_reject) {
    if(activate_accept_reject) {
      ui_widget.Accept_button->setVisible(true);
      ui_widget.Reject_button->setVisible(true);

      Q_FOREACH( QWidget* w, ui_widget.dockWidgetContents->findChildren<QWidget*>() )
      { w->setEnabled(false); }

      ui_widget.Accept_button->setEnabled(true);
      ui_widget.Reject_button->setEnabled(true);
    }
    else {
      ui_widget.Accept_button->setVisible(false);
      ui_widget.Reject_button->setVisible(false);

      Q_FOREACH( QWidget* w, ui_widget.dockWidgetContents->findChildren<QWidget*>() )
      { w->setEnabled(true); }
    }
  }

  static QString no_selected_hole_visualizer_error_message() {
    return "Error: please select a hole visualizer from Geometric Objects list."
      "Use 'Visualize Holes' button to create one by selecting the polyhedron item!";
  }

}; // end Polyhedron_demo_hole_filling_plugin

void Polyhedron_demo_hole_filling_plugin::init(QMainWindow* mainWindow,
                                      Scene_interface* scene_interface,
                                      Messages_interface* m)
{
  last_active_item = NULL;

  mw = mainWindow;
  scene = scene_interface;
  messages = m;

  actionHoleFilling = new QAction(tr("Hole Filling"), mw);
  connect(actionHoleFilling, SIGNAL(triggered()), this, SLOT(hole_filling_action()));

  dock_widget = new QDockWidget("Hole Filling", mw);
  dock_widget->setVisible(false);
  dock_widget->installEventFilter(this);

  ui_widget.setupUi(dock_widget);
  ui_widget.Accept_button->setVisible(false);
  ui_widget.Reject_button->setVisible(false);

  add_dock_widget(dock_widget);
  
  connect(ui_widget.Visualize_holes_button,  SIGNAL(clicked()), this, SLOT(on_Visualize_holes_button()));  
  connect(ui_widget.Fill_selected_holes_button,  SIGNAL(clicked()), this, SLOT(on_Fill_selected_holes_button())); 
  connect(ui_widget.Select_all_holes_button,  SIGNAL(clicked()), this, SLOT(on_Select_all_holes_button()));
  connect(ui_widget.Deselect_all_holes_button,  SIGNAL(clicked()), this, SLOT(on_Deselect_all_holes_button()));
  connect(ui_widget.Create_polyline_items_button,  SIGNAL(clicked()), this, SLOT(on_Create_polyline_items_button()));
  connect(ui_widget.Accept_button,  SIGNAL(clicked()), this, SLOT(on_Accept_button()));
  connect(ui_widget.Reject_button,  SIGNAL(clicked()), this, SLOT(on_Reject_button()));
  connect(ui_widget.Select_small_holes_button,  SIGNAL(clicked()), this, SLOT(on_Select_small_holes_button()));

  if(Scene* scene_casted = dynamic_cast<Scene*>(scene_interface)) 
  { connect(scene_casted, SIGNAL(itemAboutToBeDestroyed(Scene_item*)), this, SLOT(item_about_to_be_destroyed(Scene_item*))); }
}

void Polyhedron_demo_hole_filling_plugin::item_about_to_be_destroyed(Scene_item* scene_item) {
  Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene_item);
  if(poly_item) {
    // erase assoc polylines item
    scene->erase( scene->item_id( get_hole_visualizer(poly_item) ) );
    visualizers.remove(poly_item);
    // close accept-reject dialog if it is open
    if(last_active_item == poly_item) {
      on_Accept_button();
    }
  }
  else {
    Scene_hole_visualizer* visu_item = qobject_cast<Scene_hole_visualizer*>(scene_item);
    if(visu_item) {
        visualizers.remove(visu_item->poly_item);
    }
  }
}
// removes Scene_hole_visualizer items
void Polyhedron_demo_hole_filling_plugin::dock_widget_closed() {
  // remove all Scene_hole_visualizer items
  for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries();
    i < end; ++i)
  {
    Scene_hole_visualizer* hole_visualizer = qobject_cast<Scene_hole_visualizer*>(scene->item(i));
    if(hole_visualizer) {
      scene->erase( scene->item_id(hole_visualizer) );
    }
  }
  on_Accept_button(); 
}
// creates a Scene_hole_visualizer and associate it with active Scene_polyhedron_item
void Polyhedron_demo_hole_filling_plugin::on_Visualize_holes_button() {
  Scene_polyhedron_item* poly_item = get_selected_item<Scene_polyhedron_item>();
  if(!poly_item) {
    print_message("Error: please select a polyhedron item from Geometric Objects list!");
    return;
  }

  if(get_hole_visualizer(poly_item)) {
    print_message("Error: selected polyhedron item already has an associated hole visualizer!");
    return;
  }

  Scene_hole_visualizer* hole_visualizer = new Scene_hole_visualizer(poly_item, mw);
  visualizers[poly_item] = hole_visualizer;
  connect(hole_visualizer, SIGNAL(itemChanged()), this, SLOT(hole_visualizer_changed()));

  if(hole_visualizer->polyline_data_list.empty()) {
    print_message("There is no hole in selected polyhedron item!");
    visualizers.remove(poly_item);
    delete hole_visualizer;
    return;
  }
  else {
    // poly_item->setFlatMode(); // for better visualization
    int item_id = scene->addItem(hole_visualizer);
    scene->setSelectedItem(item_id);
    hole_visualizer->setName(tr("%1-hole visualizer").arg(poly_item->name()));
  }
}
// fills selected holes on active Scene_hole_visualizer
void Polyhedron_demo_hole_filling_plugin::on_Fill_selected_holes_button() {
  // get active polylines item
  Scene_hole_visualizer* hole_visualizer = get_selected_item<Scene_hole_visualizer>();
  if(!hole_visualizer) {
    print_message(no_selected_hole_visualizer_error_message());
    return;
  }
  if(hole_visualizer->selected_holes.empty()) {
    print_message("Error: there is no selected holes in hole visualizer!");
    return;
  }

  QApplication::setOverrideCursor(Qt::WaitCursor);
  // fill selected holes  
  int counter = 0;
  int filled_counter = 0;
  for(Scene_hole_visualizer::Selected_holes_set::iterator it = hole_visualizer->selected_holes.begin();
    it != hole_visualizer->selected_holes.end(); ++it, ++counter) {
      print_message(tr("Hole %1:").arg(counter));
      if( fill(*(hole_visualizer->poly_item->polyhedron()), (*it)->halfedge) ) { ++filled_counter;}
  }

  if(filled_counter > 0) {
    change_poly_item_by_blocking(hole_visualizer->poly_item, hole_visualizer);
    last_active_item = hole_visualizer->poly_item;
    accept_reject_toggle(true);
  }
  print_message(tr("%1 of %2 holes are filled!").arg(filled_counter).arg(counter));
  QApplication::restoreOverrideCursor();
}

// fills all holes and removes associated Scene_hole_visualizer if any
void Polyhedron_demo_hole_filling_plugin::on_Select_all_holes_button() {
  Scene_hole_visualizer* hole_visualizer = get_selected_item<Scene_hole_visualizer>();
  if(!hole_visualizer) {
    print_message(no_selected_hole_visualizer_error_message());
    return;
  }
  hole_visualizer->select_deselect_all(true);
}

void Polyhedron_demo_hole_filling_plugin::on_Select_small_holes_button() {
  Scene_hole_visualizer* hole_visualizer = get_selected_item<Scene_hole_visualizer>();
  if(!hole_visualizer) {
    print_message(no_selected_hole_visualizer_error_message());
    return;
  }

  std::size_t threshold = ui_widget.vertices_threshold_spin_box->value();
  typedef Scene_hole_visualizer::Polyline_data_list::const_iterator const_iterator;
  for(const_iterator it = hole_visualizer->polyline_data_list.begin();
                     it != hole_visualizer->polyline_data_list.end(); ++it)
  {
    if(it->polyline->polylines.front().size() <= threshold+1)
      hole_visualizer->selected_holes.insert(it);
  }
  scene->itemChanged(hole_visualizer);
}

void Polyhedron_demo_hole_filling_plugin::on_Deselect_all_holes_button() {
  Scene_hole_visualizer* hole_visualizer = get_selected_item<Scene_hole_visualizer>();
  if(!hole_visualizer) {
    print_message(no_selected_hole_visualizer_error_message());
    return;
  }
  hole_visualizer->select_deselect_all(false);
}

// Simply create polyline items and put them into scene - nothing related with other parts of the plugin
void Polyhedron_demo_hole_filling_plugin::on_Create_polyline_items_button(){
  Scene_hole_visualizer* hole_visualizer = get_selected_item<Scene_hole_visualizer>();
  if(!hole_visualizer) {
    print_message(no_selected_hole_visualizer_error_message());
    return;
  }
  if(hole_visualizer->selected_holes.empty()) {
    print_message("Error: there is no selected holes in hole visualizer!");
    return;
  }
  int counter = 0;
  for(Scene_hole_visualizer::Selected_holes_set::iterator it = hole_visualizer->selected_holes.begin();
    it != hole_visualizer->selected_holes.end(); ++it) {
      Scene_polylines_item* polyline_item = new Scene_polylines_item();
      polyline_item->polylines = (*it)->polyline->polylines;
      polyline_item->setName(QString("selected hole %1").arg(counter++));
      scene->addItem(polyline_item);
  }
}
void Polyhedron_demo_hole_filling_plugin::on_Accept_button() {
  if(last_active_item == NULL) { return; }

  accept_reject_toggle(false);
  if(Scene_hole_visualizer* hole_visualizer = get_hole_visualizer(last_active_item)) 
  { hole_visualizer->poly_item_changed();}

  new_facets.clear();
  last_active_item = NULL;
}
void Polyhedron_demo_hole_filling_plugin::on_Reject_button() {
  if(last_active_item == NULL) { return; }

  accept_reject_toggle(false);
  for(std::vector<Polyhedron::Facet_handle>::iterator it = new_facets.begin(); it != new_facets.end(); ++it) {
   last_active_item->polyhedron()->erase_facet((*it)->halfedge());
  }
  change_poly_item_by_blocking(last_active_item, get_hole_visualizer(last_active_item));

  new_facets.clear();
  last_active_item = NULL;
}
// To delete Scene_hole_visualizer when it becomes empty
void Polyhedron_demo_hole_filling_plugin::hole_visualizer_changed() {
  Scene_hole_visualizer* hole_visualizer = qobject_cast<Scene_hole_visualizer*>(this->sender());
  if(hole_visualizer && hole_visualizer->polyline_data_list.empty()) {
    scene->erase( scene->item_id(hole_visualizer));
    visualizers.remove(hole_visualizer->poly_item);
  }
}
// helper function for filling holes
bool Polyhedron_demo_hole_filling_plugin::fill
  (Polyhedron& poly, Polyhedron::Halfedge_handle it) {

  int action_index = ui_widget.action_combo_box->currentIndex();
  double alpha = ui_widget.Density_control_factor_spin_box->value();
  bool use_DT = ui_widget.Use_delaunay_triangulation_check_box->isChecked();
  unsigned int continuity = ui_widget.Continuity_spin_box->value();

  CGAL::Timer timer; timer.start();
  std::vector<Polyhedron::Facet_handle> patch;
  if(action_index == 0) {
    CGAL::Polygon_mesh_processing::triangulate_hole(poly,
             it, std::back_inserter(patch),
             CGAL::Polygon_mesh_processing::parameters::use_delaunay_triangulation(use_DT));
  }
  else if(action_index == 1) {
    CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(poly,
             it, std::back_inserter(patch), CGAL::Emptyset_iterator(),
             CGAL::Polygon_mesh_processing::parameters::density_control_factor(alpha).
             use_delaunay_triangulation(use_DT));
  }
  else {
    int weight_index = ui_widget.weight_combo_box->currentIndex();

    bool success;
    if(weight_index == 0) {
      success = CGAL::cpp11::get<0>(CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(poly,
              it, std::back_inserter(patch), CGAL::Emptyset_iterator(),
              CGAL::Polygon_mesh_processing::parameters::weight_calculator
                (CGAL::internal::Uniform_weight_fairing<Polyhedron>(poly)).
              density_control_factor(alpha).
              fairing_continuity(continuity).
              use_delaunay_triangulation(use_DT)));
    }
    else {
      success = CGAL::cpp11::get<0>(CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(poly,
              it, std::back_inserter(patch), CGAL::Emptyset_iterator(),
              CGAL::Polygon_mesh_processing::parameters::weight_calculator(CGAL::internal::Cotangent_weight_with_voronoi_area_fairing<Polyhedron>(poly)).
              density_control_factor(alpha).
              fairing_continuity(continuity).
              use_delaunay_triangulation(use_DT)));
    }

    if(!success) { print_message("Error: fairing is not successful, only triangulation and refinement are applied!"); }
  }
  print_message(QString("Took %1 sec.").arg(timer.time()));

  if(patch.empty()) {
    print_message(tr("Warning: generating patch is not successful! %1")
      .arg(use_DT ? "Please try without 'Use 3D Delaunay Triangulation'!" : ""));
    return false;
  }

  // Self intersection test
  if(ui_widget.Skip_self_intersection_check_box->checkState() == Qt::Checked) {
    timer.reset();

    typedef std::vector<std::pair<Polyhedron::Facet_const_handle, Polyhedron::Facet_const_handle> > Intersected_facets;
    Intersected_facets intersected_facets;
    CGAL::Polygon_mesh_processing::self_intersections(poly,
      std::back_inserter(intersected_facets),
      CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, poly)));

    print_message(QString("Self intersecting test: finding intersecting triangles in %1 sec.").arg(timer.time()));
    timer.reset();
    // this part might need speed-up
    bool intersected = false;
    for(Intersected_facets::iterator it = intersected_facets.begin(); 
      it != intersected_facets.end() && !intersected; ++it) {
      for(std::vector<Polyhedron::Facet_handle>::iterator it_patch = patch.begin(); 
        it_patch != patch.end() && !intersected; ++it_patch) {
        if(it->first == (*it_patch) || it->second == (*it_patch)) {
          intersected = true;
        }
      }
    }
    print_message(QString("Self intersecting test: iterate on patch in %1 sec.").arg(timer.time()));
    if(intersected) {
      for(std::vector<Polyhedron::Facet_handle>::iterator it = patch.begin(); it != patch.end(); ++it) {
        poly.erase_facet((*it)->halfedge());
      }
      print_message("Self intersecting patch is generated, and it is removed.");
      return false;
    }
    else { print_message("No Self intersection found, patch is valid."); }
  }
  // save facets for accept-reject 
  new_facets.insert(new_facets.end(), patch.begin(), patch.end());
  return true;
}

// Q_EXPORT_PLUGIN2(Polyhedron_demo_hole_filling_plugin, Polyhedron_demo_hole_filling_plugin)

#include "Polyhedron_demo_hole_filling_plugin.moc"
