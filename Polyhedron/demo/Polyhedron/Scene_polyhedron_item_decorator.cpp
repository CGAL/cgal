#include "Scene_polyhedron_item_decorator.h"


Scene_polyhedron_item_decorator::Scene_polyhedron_item_decorator
  (Scene_face_graph_item* poly_item, bool delete_item)
  : poly_item(poly_item), delete_poly_item(delete_item)
{ }

Scene_polyhedron_item_decorator::~Scene_polyhedron_item_decorator()
{
  if(delete_poly_item) { delete poly_item; }
}

Scene_polyhedron_item_decorator*
Scene_polyhedron_item_decorator::clone() const {
  return 0;
}

QString
Scene_polyhedron_item_decorator::toolTip() const
{
  if(!poly_item->polyhedron())
    return QString();

  return QObject::tr("<p>Selection <b>%1</b> (mode: %5, color: %6)</p>"
                     "<p>Number of vertices: %2<br />"
                     "Number of edges: %3<br />"
                     "Number of faces: %4</p>")
    .arg(this->name())
    .arg(num_vertices(*poly_item->polyhedron()))
    .arg(num_edges(*poly_item->polyhedron()))
    .arg(num_faces(*poly_item->polyhedron()))
    .arg(this->renderingModeName())
    .arg(this->color().name());
}

Face_graph*
Scene_polyhedron_item_decorator::polyhedron()
{ return poly_item->polyhedron(); }

const Face_graph*
Scene_polyhedron_item_decorator::polyhedron() const
{ return poly_item->polyhedron(); }

bool
Scene_polyhedron_item_decorator::isEmpty() const {
  return poly_item->isEmpty();
}

void
Scene_polyhedron_item_decorator::compute_bbox() const {
  setBbox(poly_item->bbox());
}


void
Scene_polyhedron_item_decorator::
invalidateOpenGLBuffers()
{
  poly_item->invalidateOpenGLBuffers();
  Scene_item::invalidateOpenGLBuffers();
  compute_bbox();
}

void
Scene_polyhedron_item_decorator::select(double orig_x,
                                   double orig_y,
                                   double orig_z,
                                   double dir_x,
                                   double dir_y,
                                   double dir_z)
{
  poly_item->select(orig_x,
                    orig_y,
                    orig_z,
                    dir_x,
                    dir_y,
                    dir_z);
}

Scene_face_graph_item* Scene_polyhedron_item_decorator::polyhedron_item() const {
  return poly_item;
}

void Scene_polyhedron_item_decorator::set_polyhedron_item(Scene_face_graph_item* poly_item) {
  this->poly_item = poly_item;
}

