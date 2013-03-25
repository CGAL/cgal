#include "Scene_edit_polyhedron_item_2.h"
#include "Kernel_type.h"
#include "Polyhedron_type.h"

#include <boost/foreach.hpp>
#include <algorithm>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include "Property_maps_for_edit_plugin.h"
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <QVariant>

#include <QObject>
#include <QMenu>
#include <QAction>

#include <QEvent>
#include <QKeyEvent>
#include <QMouseEvent>

#include <CGAL/gl_render.h>

#include "ui_Deform_mesh_2.h"

Scene_edit_polyhedron_item_2::Scene_edit_polyhedron_item_2(Scene_polyhedron_item* poly_item)
  : poly_item(poly_item), deform_mesh(*(poly_item->polyhedron()), Vertex_index_map(), Edge_index_map())
  , show_roi(true), ui_widget(NULL), frame(new qglviewer::ManipulatedFrame())
{
  // bind vertex picking 
  connect(poly_item, SIGNAL(selected_vertex(void*)), this, SLOT(vertex_has_been_selected(void*)));
  poly_item->enable_facets_picking(true);

  // interleave events of viewer (there is only one viewer) 
  QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
  viewer->installEventFilter(this);
    
  // create an empty handle group for starting
  create_handle_group();
   
  // start QObject's timer for continous effects 
  // (like selecting vertices as 'painting', deforming mesh while mouse not moving)
  startTimer(0);
}

Scene_edit_polyhedron_item_2::~Scene_edit_polyhedron_item_2()
{
  delete frame;
  while(is_there_any_handle_group())
  {
    delete_handle_group();
  }
}

struct Get_vertex_handle : public CGAL::Modifier_base<Polyhedron::HDS>
{
  Polyhedron::Vertex* vertex_ptr;
  vertex_descriptor vh;
  void operator()(Polyhedron::HDS& hds) {
    vh = hds.vertex_handle(vertex_ptr);
  }
};
/////////////////////////////////////////////////////////
/////////// Most relevant functions lie here ///////////
void Scene_edit_polyhedron_item_2::deform()
{
  if(handles.empty()) { return; }

  Deform_mesh::Handle_group hgb, hge;
  for(boost::tie(hgb, hge) = deform_mesh.handle_groups(); hgb != hge; ++hgb)
  {
    Handle_group_data& hd = get_data(hgb);
    qglviewer::ManipulatedFrame* hgb_frame = hd.frame;
    qglviewer::Vec translation = hgb_frame->position() - hd.initial_center;

    deform_mesh.rotate(hgb, 
      Point(hd.initial_center.x, hd.initial_center.y, hd.initial_center.z),
      hgb_frame->orientation(),
      translation);
  }
  deform_mesh.deform();

  emit mesh_deformed(this);
}
void Scene_edit_polyhedron_item_2::vertex_has_been_selected(void* void_ptr) {

  Polyhedron* poly = poly_item->polyhedron();
  // get vertex descriptor
  Get_vertex_handle get_vertex_handle;  
  get_vertex_handle.vertex_ptr = static_cast<Polyhedron::Vertex*>(void_ptr);
  poly->delegate(get_vertex_handle);
  vertex_descriptor clicked_vertex = get_vertex_handle.vh;
  // use clicked_vertex, do what you want 
  bool is_roi = ui_widget->ROIRadioButton->isChecked();
  bool is_insert = ui_widget->InsertRadioButton->isChecked();
  int k_ring = ui_widget->BrushSpinBox->value();
  process_selection(clicked_vertex, k_ring, is_roi, is_insert);
}
void Scene_edit_polyhedron_item_2::timerEvent(QTimerEvent *event)
{
  if(state.ctrl_pressing && (state.left_button_pressing || state.right_button_pressing) ) {
    deform(); 
  }
}
bool Scene_edit_polyhedron_item_2::eventFilter(QObject *target, QEvent *event)
{
  bool need_repaint = false;
  Mouse_keyboard_state old_state = state;
  ////////////////// TAKE EVENTS /////////////////////
  // key events
  if(event->type() == QEvent::KeyPress || event->type() == QEvent::KeyRelease) 
  {
    QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
    Qt::KeyboardModifiers modifiers = keyEvent->modifiers();

    state.ctrl_pressing = modifiers.testFlag(Qt::ControlModifier);
    state.shift_pressing = modifiers.testFlag(Qt::ShiftModifier);
  }
  // mouse events
  if(event->type() == QEvent::MouseButtonPress || event->type() == QEvent::MouseButtonRelease)
	{
    QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
    if(mouse_event->button() == Qt::LeftButton) {
      state.left_button_pressing = event->type() == QEvent::MouseButtonPress;
    }
    if(mouse_event->button() == Qt::RightButton) {
      state.right_button_pressing = event->type() == QEvent::MouseButtonPress;
    }    
  }
  ////////////////// //////////////// /////////////////////

  // check state changes between old and current state
  bool ctrl_pressed_now = state.ctrl_pressing && !old_state.ctrl_pressing;
  bool ctrl_released_now = !state.ctrl_pressing && old_state.ctrl_pressing;
  if(ctrl_pressed_now || ctrl_released_now) 
  {// update manipulated frame autoselection because they are not updated by QGLViewer until mouse moves
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    activate_closest_manipulated_frame(viewer->mapFromGlobal(QCursor::pos())); 

    need_repaint = true;
  }

  // use mouse move event for paint-like selection
  // specificly placed in here (not in timer function) for preventing unneccessary ray-casting 
  // (while mouse not moving, ray casting is pointless since the same vertex will be returned and processed)
  if(event->type() == QEvent::MouseMove)    
  { 
    // paint with mouse move event
    if(state.shift_pressing && state.left_button_pressing)
    {    
      QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
      QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
      qglviewer::Camera* camera = viewer->camera();

      bool found = false;
      const qglviewer::Vec& point = camera->pointUnderPixel(mouse_event->pos(), found);
      if(found)
      {
        const qglviewer::Vec& orig = camera->position();
        const qglviewer::Vec& dir = point - orig;
        select(orig.x, orig.y, orig.z, dir.x, dir.y, dir.z); // it will cause 'selected_vertex' signal where our slot will handle selection 

        need_repaint = true;
      }
    }// end shift_pressing
  }//end MouseMove

  if(need_repaint) { emit mesh_repaint_needed(this); }

  return false;
}
#include "opengl_tools.h"
void Scene_edit_polyhedron_item_2::draw() const {
  poly_item->direct_draw();
  CGAL::GL::Point_size point_size; point_size.set_point_size(5);
  CGAL::GL::Color color; color.set_rgb_color(0, 1.f, 0);
  if(show_roi) {
    // draw ROI
    ::glBegin(GL_POINTS);
    for(std::set<vertex_descriptor>::const_iterator it = roi.begin(); it != roi.end(); ++it)
    {
      if(handles.find(*it) == handles.end())
      {
        const Kernel::Point_3& p = (*it)->point();
        ::glVertex3d(p.x(), p.y(), p.z());
      }
    }
    ::glEnd();
  }
  // draw handles
  Deform_mesh::Const_handle_group hgb, hge;
  for(boost::tie(hgb, hge) = deform_mesh.handle_groups(); hgb != hge; ++hgb)
  {
    // draw axis using manipulated frame assoc with handle_group
    const Handle_group_data& hgb_data = get_data(hgb);
    if(hgb_data.frame->grabsMouse())
    {      
      // draw axis
      ::glPushMatrix();
        ::glMultMatrixd(hgb_data.frame->matrix());
        QGLViewer::drawAxis(0.1f);
      ::glPopMatrix();
      // draw bbox
      color.set_rgb_color(1.0f, 0, 0);
      ::glPushMatrix();
        ::glTranslated(hgb_data.frame->position().x, hgb_data.frame->position().y, hgb_data.frame->position().z);
        ::glMultMatrixd(hgb_data.frame->orientation().matrix());
        ::glTranslated(-hgb_data.initial_center.x, -hgb_data.initial_center.y, -hgb_data.initial_center.z);        
        draw_bbox(hgb_data.bbox);
      ::glPopMatrix();
    }

    if(hgb == active_group) { color.set_rgb_color(1.0f, 0, 0); }
    else                    { color.set_rgb_color(0, 0, 1.0f); }
    ::glBegin(GL_POINTS);
      Deform_mesh::Const_handle_iterator hb, he;
      for(boost::tie(hb, he) = deform_mesh.handles(hgb); hb != he; ++hb)
      {           
        const Kernel::Point_3& p = (*hb)->point();
        ::glVertex3d(p.x(), p.y(), p.z());
      }
    ::glEnd();
  }
}
//////////////////////////////////////////////////////////

/////////////// from trivial_plugin //////////////////////
void Scene_edit_polyhedron_item_2::draw_bbox(const Scene_interface::Bbox& bb ) const {
  ::glBegin(GL_LINES);
  gl_draw_edge(bb.xmin, bb.ymin, bb.zmin,
                bb.xmax, bb.ymin, bb.zmin);
  gl_draw_edge(bb.xmin, bb.ymin, bb.zmin,
                bb.xmin, bb.ymax, bb.zmin);
  gl_draw_edge(bb.xmin, bb.ymin, bb.zmin,
                bb.xmin, bb.ymin, bb.zmax);
    
  gl_draw_edge(bb.xmax, bb.ymin, bb.zmin,
                bb.xmax, bb.ymax, bb.zmin);
  gl_draw_edge(bb.xmax, bb.ymin, bb.zmin,
                bb.xmax, bb.ymin, bb.zmax);
    
  gl_draw_edge(bb.xmin, bb.ymax, bb.zmin,
                bb.xmax, bb.ymax, bb.zmin);
  gl_draw_edge(bb.xmin, bb.ymax, bb.zmin,
                bb.xmin, bb.ymax, bb.zmax);
    
  gl_draw_edge(bb.xmin, bb.ymin, bb.zmax,
                bb.xmax, bb.ymin, bb.zmax);
  gl_draw_edge(bb.xmin, bb.ymin, bb.zmax,
                bb.xmin, bb.ymax, bb.zmax);
    
  gl_draw_edge(bb.xmax, bb.ymax, bb.zmax,
                bb.xmin, bb.ymax, bb.zmax);
  gl_draw_edge(bb.xmax, bb.ymax, bb.zmax,
                bb.xmax, bb.ymin, bb.zmax);
  gl_draw_edge(bb.xmax, bb.ymax, bb.zmax,
                bb.xmax, bb.ymax, bb.zmin);
  ::glEnd();
}
void Scene_edit_polyhedron_item_2::gl_draw_edge(double px, double py, double pz,
                          double qx, double qy, double qz) const
{
  ::glVertex3d(px,py,pz);
  ::glVertex3d(qx,qy,qz);
}
/////////////////////////////////////////////////////////////

void Scene_edit_polyhedron_item_2::changed()
{
  poly_item->changed();
  Scene_item::changed();
//  last_pos = current_position();
}
Scene_polyhedron_item* Scene_edit_polyhedron_item_2::to_polyhedron_item() const {
  return poly_item;
}
qglviewer::ManipulatedFrame* Scene_edit_polyhedron_item_2::manipulatedFrame() {
  return frame;
}

Polyhedron* Scene_edit_polyhedron_item_2::polyhedron()       
{ return poly_item->polyhedron(); }
const Polyhedron* Scene_edit_polyhedron_item_2::polyhedron() const 
{ return poly_item->polyhedron(); }
QString Scene_edit_polyhedron_item_2::toolTip() const
{
  if(!poly_item->polyhedron())
    return QString();

  return QObject::tr("<p>Polyhedron <b>%1</b> (mode: %5, color: %6)</p>"
                     "<p>Number of vertices: %2<br />"
                     "Number of edges: %3<br />"
                     "Number of facets: %4</p>")
    .arg(this->name())
    .arg(poly_item->polyhedron()->size_of_vertices())
    .arg(poly_item->polyhedron()->size_of_halfedges()/2)
    .arg(poly_item->polyhedron()->size_of_facets())
    .arg(this->renderingModeName())
    .arg(this->color().name());
}
bool Scene_edit_polyhedron_item_2::isEmpty() const {
  return poly_item->isEmpty();
}
Scene_edit_polyhedron_item_2::Bbox Scene_edit_polyhedron_item_2::bbox() const {
  return poly_item->bbox();
}

void Scene_edit_polyhedron_item_2::setVisible(bool b) {
  poly_item->setVisible(b);
  Scene_item::setVisible(b);
}
void Scene_edit_polyhedron_item_2::setColor(QColor c) {
  poly_item->setColor(c);
  Scene_item::setColor(c);
}
void Scene_edit_polyhedron_item_2::setName(QString n) {
  Scene_item::setName(n);
  n.replace(" (edit)", "");
  poly_item->setName(n);
}
void Scene_edit_polyhedron_item_2::setRenderingMode(RenderingMode m) {
  poly_item->setRenderingMode(m);
  Scene_item::setRenderingMode(m);
}
Scene_edit_polyhedron_item_2* Scene_edit_polyhedron_item_2::clone() const {
  return 0;
}
void Scene_edit_polyhedron_item_2::select(
          double orig_x,
          double orig_y,
          double orig_z,
          double dir_x,
          double dir_y,
          double dir_z)
{
  Scene_item::select(orig_x,
                     orig_y,
                     orig_z,
                     dir_x,
                     dir_y,
                     dir_z);
  poly_item->select(orig_x,
                       orig_y,
                       orig_z,
                       dir_x,
                       dir_y,
                       dir_z);
}
#include "Scene_edit_polyhedron_item_2.moc"
