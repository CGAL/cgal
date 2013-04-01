#include "Scene_polyhedron_transform_item.h"
#include "Kernel_type.h"
#include "Polyhedron_type.h"

Scene_polyhedron_transform_item::Scene_polyhedron_transform_item(const qglviewer::Vec& pos,const Scene_polyhedron_item* poly_item_,const Scene_interface* scene_interface):
  poly_item(poly_item_),
  scene(scene_interface),
  manipulable(false),
  can_clone(true),
  frame(new ManipulatedFrame()),
  poly(poly_item->polyhedron()),
  center_(pos) { frame->setPosition(pos); }

void Scene_polyhedron_transform_item::draw() const{
  glPushMatrix();
  glMultMatrixd(frame->matrix());
  direct_draw_edges();
  //Scene_item_with_display_list::draw();
  glPopMatrix();  
}
    
void Scene_polyhedron_transform_item::direct_draw_edges() const {
  typedef Kernel::Point_3		Point;
  typedef Polyhedron::Edge_const_iterator	Edge_iterator;

  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINES);
  Edge_iterator he;
  for(he = poly->edges_begin();
      he != poly->edges_end();
      he++)
  {
    const Point& a = he->vertex()->point();
    const Point& b = he->opposite()->vertex()->point();
    ::glVertex3d(a.x()-center_.x,a.y()-center_.y,a.z()-center_.z);
    ::glVertex3d(b.x()-center_.x,b.y()-center_.y,b.z()-center_.z);
  }
  ::glEnd();
  ::glEnable(GL_LIGHTING);
}
QString Scene_polyhedron_transform_item::toolTip() const {
  return QObject::tr("<p>Affine transformation of <b>%1</b></p>"
                     "<p>Keep <b>Ctrl</b> pressed and use the arcball to define an affine transformation.<br />"
                     "Press <b>S</b> to apply the affine transformation to a copy of <b>%1</b>.</p>")
    .arg(getBase()->name());
}
bool Scene_polyhedron_transform_item::keyPressEvent(QKeyEvent* e){
  if (e->key()==Qt::Key_S){
    emit stop();
    return true;
  }
  return false;
}

Scene_polyhedron_transform_item::Bbox
Scene_polyhedron_transform_item::bbox() const {
  const Kernel::Point_3& p = *(poly->points_begin());
  CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
  for(Polyhedron::Point_const_iterator it = poly->points_begin();
      it != poly->points_end();
      ++it) {
    bbox = bbox + it->bbox();
  }
  return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
              bbox.xmax(),bbox.ymax(),bbox.zmax());
}

#include "Scene_polyhedron_transform_item.moc"

