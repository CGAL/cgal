#include "Scene_edit_polyhedron_item_2.h"
#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include "Custom_manipulated_frame.h"

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
#include <CGAL/gl_render.h>

#include "ui_Deform_mesh_2.h"

Scene_edit_polyhedron_item_2::Scene_edit_polyhedron_item_2(Scene_polyhedron_item* poly_item)
  : poly_item(poly_item), deform_mesh(*(poly_item->polyhedron()), Vertex_index_map(), Edge_index_map())
  , active_group(deform_mesh.create_handle_group())
{
  frame = new CustomManipulatedFrame();
  frame->setProperty("item", QVariant::fromValue<QObject*>(this));

  connect(poly_item, SIGNAL(selected_vertex(void*)), this, SLOT(vertex_has_been_selected(void*)));
  poly_item->enable_facets_picking(true);

  connect(frame, SIGNAL(modified()), this, SIGNAL(modified()));    
}

Scene_edit_polyhedron_item_2::~Scene_edit_polyhedron_item_2()
{ }

struct Get_vertex_handle : public CGAL::Modifier_base<Polyhedron::HDS>
{
  Polyhedron::Vertex* vertex_ptr;
  vertex_descriptor vh;
  void operator()(Polyhedron::HDS& hds) {
    vh = hds.vertex_handle(vertex_ptr);
  }
};

void Scene_edit_polyhedron_item_2::vertex_has_been_selected(void* void_ptr) {
  Polyhedron* poly = poly_item->polyhedron();

  Get_vertex_handle get_vertex_handle;  
  get_vertex_handle.vertex_ptr = static_cast<Polyhedron::Vertex*>(void_ptr);

  poly->delegate(get_vertex_handle);
  vertex_descriptor clicked_vertex = get_vertex_handle.vh;

  bool is_roi = ui_widget->ROIRadioButton->isChecked();
  bool is_insert = ui_widget->InsertRadioButton->isChecked();
  int k_ring = ui_widget->BrushSpinBox->value();
  process_selection(clicked_vertex, k_ring, is_roi, is_insert);
}

#include "opengl_tools.h"
void Scene_edit_polyhedron_item_2::draw() const {
  poly_item->direct_draw();
  CGAL::GL::Point_size point_size; point_size.set_point_size(5);
  CGAL::GL::Color color; color.set_rgb_color(0, 1.f, 0);
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
  // draw handles
  Deform_mesh::Const_handle_group hgb, hge;
  for(boost::tie(hgb, hge) = deform_mesh.handle_groups(); hgb != hge; ++hgb)
  {
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
Polyhedron* Scene_edit_polyhedron_item_2::polyhedron()       
{ return poly_item->polyhedron(); }
const Polyhedron* Scene_edit_polyhedron_item_2::polyhedron() const 
{ return poly_item->polyhedron(); }
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
