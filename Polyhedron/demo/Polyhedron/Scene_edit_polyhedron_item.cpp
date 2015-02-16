#include "opengl_tools.h"
#include "Scene_edit_polyhedron_item.h"
#include <boost/foreach.hpp>
#include <algorithm>


#include <CGAL/gl_render.h>

Scene_edit_polyhedron_item::Scene_edit_polyhedron_item
  (Scene_polyhedron_item* poly_item, 
  Ui::DeformMesh* ui_widget,
  QMainWindow* mw)
  : ui_widget(ui_widget), 
    poly_item(poly_item),
    deform_mesh(*(poly_item->polyhedron()), Deform_mesh::Vertex_index_map(), Deform_mesh::Hedge_index_map(), Array_based_vertex_point_map(&positions)),
    is_rot_free(true),
    own_poly_item(true),
    k_ring_selector(poly_item, mw, Scene_polyhedron_item_k_ring_selection::Active_handle::VERTEX, true),
    quadric(gluNewQuadric())
{
  mw->installEventFilter(this);
  gluQuadricNormals(quadric, GLU_SMOOTH);
  // bind vertex picking 
  connect(&k_ring_selector, SIGNAL(selected(const std::set<Polyhedron::Vertex_handle>&)), this,
    SLOT(selected(const std::set<Polyhedron::Vertex_handle>&)));

  poly_item->set_color_vector_read_only(true); // to prevent recomputation of color vector in changed()
  poly_item->update_vertex_indices();

  length_of_axis = bbox().diagonal_length() / 15.0;

  // interleave events of viewer (there is only one viewer) 
  QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
  viewer->installEventFilter(this);
    
  // create an empty group of control vertices for starting
  create_ctrl_vertices_group();
   
  // start QObject's timer for continuous effects 
  // (deforming mesh while mouse not moving)
  startTimer(0);

  // Required for drawing functionality
  positions.resize(num_vertices(*polyhedron())*3);
  normals.resize(positions.size());
  Polyhedron::Vertex_iterator vb, ve;
  std::size_t counter = 0;
  for(vb=polyhedron()->vertices_begin(), ve = polyhedron()->vertices_end();vb != ve; ++vb, ++counter) {
    positions[counter*3] = vb->point().x();
    positions[counter*3+1] = vb->point().y();
    positions[counter*3+2] = vb->point().z();

    const Polyhedron::Traits::Vector_3& n = 
      CGAL::Polygon_mesh_processing::compute_vertex_normal<Polyhedron::Traits>(vb, deform_mesh.halfedge_graph());

    normals[counter*3] = n.x();
    normals[counter*3+1] = n.y();
    normals[counter*3+2] = n.z();
  }

  tris.resize(polyhedron()->size_of_facets()*3);
  counter = 0;
  for(Polyhedron::Facet_handle fb = polyhedron()->facets_begin(); fb != polyhedron()->facets_end(); ++fb, ++counter) {
    tris[counter*3] =  static_cast<unsigned int>(fb->halfedge()->vertex()->id());
    tris[counter*3+1] = static_cast<unsigned int>(fb->halfedge()->next()->vertex()->id());
    tris[counter*3+2] = static_cast<unsigned int>(fb->halfedge()->prev()->vertex()->id());
  }

  edges.resize(polyhedron()->size_of_halfedges());
  counter = 0;
  for(Polyhedron::Edge_iterator eb = polyhedron()->edges_begin(); eb != polyhedron()->edges_end(); ++eb, ++counter) {
    edges[counter*2] = static_cast<unsigned int>(eb->vertex()->id());
    edges[counter*2+1] = static_cast<unsigned int>(eb->opposite()->vertex()->id());
  }
}

Scene_edit_polyhedron_item::~Scene_edit_polyhedron_item()
{
  while(is_there_any_ctrl_vertices_group())
  {
    delete_ctrl_vertices_group(false);
  }
  gluDeleteQuadric(quadric);
  if (own_poly_item) delete poly_item;
}

/////////////////////////////////////////////////////////
/////////// Most relevant functions lie here ///////////
void Scene_edit_polyhedron_item::deform()
{
  if(!is_there_any_ctrl_vertices()) { return; }

  for(Ctrl_vertices_group_data_list::iterator it = ctrl_vertex_frame_map.begin(); it != ctrl_vertex_frame_map.end(); ++it)
  { it->set_target_positions(); }
  deform_mesh.deform();

  poly_item->changed(); // now we need to call poly_item changed to delete AABB tree 
  emit itemChanged();
}

void Scene_edit_polyhedron_item::timerEvent(QTimerEvent* /*event*/)
{ // just handle deformation - paint like selection is handled in eventFilter()
  if(state.ctrl_pressing && (state.left_button_pressing || state.right_button_pressing)) {
    if(!ui_widget->ActivatePivotingCheckBox->isChecked()) {
      deform();
    }
    else {
      emit itemChanged(); // for redraw while Pivoting (since we close signals of manipulatedFrames while pivoting, 
                          // for now redraw with timer)
    }
  }
}
bool Scene_edit_polyhedron_item::eventFilter(QObject* /*target*/, QEvent *event)
{
  // This filter is both filtering events from 'viewer' and 'main window'
  Mouse_keyboard_state_deformation old_state = state;
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

  if(!poly_item->visible()) { return false; } // if not visible just update event state but don't do any action

  // check state changes between old and current state
  bool ctrl_pressed_now = state.ctrl_pressing && !old_state.ctrl_pressing;
  bool ctrl_released_now = !state.ctrl_pressing && old_state.ctrl_pressing;
  if(ctrl_pressed_now || ctrl_released_now || event->type() == QEvent::HoverMove) 
  {// activate a handle manipulated frame
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    const QPoint& p = viewer->mapFromGlobal(QCursor::pos());
    bool need_repaint = activate_closest_manipulated_frame(p.x(), p.y());

    if(need_repaint) { emit itemChanged(); }
  }

  return false;
}

#include "opengl_tools.h"
void Scene_edit_polyhedron_item::draw_edges() const {

  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_DOUBLE, 0, positions.data());
  glDrawElements(GL_LINES, (GLsizei) edges.size(), GL_UNSIGNED_INT, edges.data());
  glDisableClientState(GL_VERTEX_ARRAY); 

  if(rendering_mode == Wireframe) {
    draw_ROI_and_control_vertices();
  }
}
void Scene_edit_polyhedron_item::draw() const {
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);

  glVertexPointer(3, GL_DOUBLE, 0, positions.data());
  glNormalPointer(GL_DOUBLE, 0, normals.data());
  glDrawElements(GL_TRIANGLES, (GLsizei) tris.size(), GL_UNSIGNED_INT, tris.data());

  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);

  CGAL::GL::Color color;
  color.set_rgb_color(0, 0, 0);
  draw_edges();

  draw_ROI_and_control_vertices();
}

void Scene_edit_polyhedron_item::draw_ROI_and_control_vertices() const {
  GLboolean enable_back_lighting = glIsEnabled(GL_LIGHTING);
  glDisable(GL_LIGHTING);

  CGAL::GL::Color color;
  CGAL::GL::Point_size point_size; point_size.set_point_size(5);
  color.set_rgb_color(0, 1.f, 0);
  // draw ROI
  if(ui_widget->ShowROICheckBox->isChecked()) {
    BOOST_FOREACH(vertex_descriptor vd, deform_mesh.roi_vertices())
    {
      if(!deform_mesh.is_control_vertex(vd))
        gl_draw_point( vd->point() );
    }
  }
  // draw control vertices related things
  QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();

  for(Ctrl_vertices_group_data_list::const_iterator hgb_data = ctrl_vertex_frame_map.begin(); hgb_data != ctrl_vertex_frame_map.end(); ++hgb_data)
  {
    if(hgb_data->frame == viewer->manipulatedFrame())
    {      
      // draw axis
      ::glPushMatrix();
      ::glMultMatrixd(hgb_data->frame->matrix());
      QGLViewer::drawAxis(length_of_axis);
      ::glPopMatrix();
      // draw bbox
      if(!ui_widget->ActivatePivotingCheckBox->isChecked())
      {
        color.set_rgb_color(1.0f, 0, 0);
        ::glPushMatrix();
        ::glTranslated(hgb_data->frame->position().x, hgb_data->frame->position().y, hgb_data->frame->position().z);
        ::glMultMatrixd(hgb_data->frame->orientation().matrix());
        ::glTranslated(-hgb_data->frame_initial_center.x, -hgb_data->frame_initial_center.y, -hgb_data->frame_initial_center.z);        
        draw_bbox(hgb_data->bbox);
        ::glPopMatrix();
      }
    }
    // draw control vertices
    if(hgb_data == active_group) { color.set_rgb_color(1.0f, 0, 0); }
    else                    { color.set_rgb_color(0, 0, 1.0f); }
    for(std::vector<vertex_descriptor>::const_iterator hb = hgb_data->ctrl_vertices_group.begin(); hb != hgb_data->ctrl_vertices_group.end(); ++hb)
    {  gl_draw_point( (*hb)->point() );
    }
  }

  if(enable_back_lighting) { glEnable(GL_LIGHTING); }
}
void Scene_edit_polyhedron_item::gl_draw_point(const Point& p) const
{
  if(!ui_widget->ShowAsSphereCheckBox->isChecked()) {
    ::glBegin(GL_POINTS);
      ::glVertex3d(p.x(), p.y(), p.z());
    ::glEnd();
  } 
  else {
    GLint shading;
    ::glGetIntegerv(GL_SHADE_MODEL, &shading);
    ::glShadeModel(GL_SMOOTH);

    ::glPushMatrix();
      ::glTranslated(p.x(), p.y(), p.z());
      ::gluSphere(quadric, length_of_axis/15, 8, 8);
    ::glPopMatrix();

    ::glShadeModel(shading);
  }
}
//////////////////////////////////////////////////////////

/////////////// from trivial_plugin //////////////////////
void Scene_edit_polyhedron_item::draw_bbox(const Scene_interface::Bbox& bb ) const {
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
void Scene_edit_polyhedron_item::gl_draw_edge(double px, double py, double pz,
                          double qx, double qy, double qz) const
{
  ::glVertex3d(px,py,pz);
  ::glVertex3d(qx,qy,qz);
}
/////////////////////////////////////////////////////////////

void Scene_edit_polyhedron_item::changed()
{ update_normals(); }

Scene_polyhedron_item* Scene_edit_polyhedron_item::to_polyhedron_item() {
  Scene_polyhedron_item* poly_item_tmp = poly_item;
  poly_item->set_color_vector_read_only(false);
  own_poly_item=false;
  return poly_item_tmp;
}

Polyhedron* Scene_edit_polyhedron_item::polyhedron()       
{ return poly_item->polyhedron(); }
const Polyhedron* Scene_edit_polyhedron_item::polyhedron() const 
{ return poly_item->polyhedron(); }
QString Scene_edit_polyhedron_item::toolTip() const
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
bool Scene_edit_polyhedron_item::isEmpty() const {
  return poly_item->isEmpty();
}
Scene_edit_polyhedron_item::Bbox Scene_edit_polyhedron_item::bbox() const {
  return poly_item->bbox();
}

void Scene_edit_polyhedron_item::setVisible(bool b) {
  poly_item->setVisible(b);
  Scene_item::setVisible(b);
  if(!b) {
    (*QGLViewer::QGLViewerPool().begin())->setManipulatedFrame(NULL);
  }
}
void Scene_edit_polyhedron_item::setColor(QColor c) {
  poly_item->setColor(c);
  Scene_item::setColor(c);
}
void Scene_edit_polyhedron_item::setName(QString n) {
  Scene_item::setName(n);
  n.replace(" (edit)", "");
  poly_item->setName(n);
}
void Scene_edit_polyhedron_item::setRenderingMode(RenderingMode m) {
  poly_item->setRenderingMode(m);
  Scene_item::setRenderingMode(m);
}
Scene_edit_polyhedron_item* Scene_edit_polyhedron_item::clone() const {
  return 0;
}
void Scene_edit_polyhedron_item::select(
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

bool Scene_edit_polyhedron_item::keyPressEvent(QKeyEvent* e)
{
  //setting/unsetting rotation constraints
  if (e->key()==Qt::Key_R && !state.ctrl_pressing)
  {
    is_rot_free = !is_rot_free;
    rot_constraint.setRotationConstraintType( is_rot_free?
        qglviewer::AxisPlaneConstraint::FREE:
        qglviewer::AxisPlaneConstraint::AXIS);
    return true;
  }
  return false;
}

#include "Scene_edit_polyhedron_item.moc"
