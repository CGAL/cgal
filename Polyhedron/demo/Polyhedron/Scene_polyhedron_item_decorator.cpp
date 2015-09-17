#include "Scene_polyhedron_item_decorator.h"
#include "Polyhedron_type.h"

Scene_polyhedron_item_decorator::Scene_polyhedron_item_decorator
  (Scene_polyhedron_item* poly_item, bool delete_item)
  :Scene_item(), poly_item(poly_item), delete_poly_item(delete_item)
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

void Scene_polyhedron_item_decorator::draw() const {
  poly_item->draw();
}

void Scene_polyhedron_item_decorator::draw_edges() const {
  poly_item->draw_edges();
}

Polyhedron* 
Scene_polyhedron_item_decorator::polyhedron()       
{ return poly_item->polyhedron(); }

const Polyhedron* 
Scene_polyhedron_item_decorator::polyhedron() const 
{ return poly_item->polyhedron(); }

bool
Scene_polyhedron_item_decorator::isEmpty() const {
  return poly_item->isEmpty();
}

Scene_polyhedron_item_decorator::Bbox
Scene_polyhedron_item_decorator::bbox() const {
  return poly_item->bbox();
}


void
Scene_polyhedron_item_decorator::
invalidate_buffers()
{
  poly_item->invalidate_buffers();
  Scene_item::invalidate_buffers();
}

void 
Scene_polyhedron_item_decorator::select(double orig_x,
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

Scene_polyhedron_item* Scene_polyhedron_item_decorator::polyhedron_item() const {
  return poly_item;
}

void Scene_polyhedron_item_decorator::set_polyhedron_item(Scene_polyhedron_item* poly_item) {
  this->poly_item = poly_item;
}

