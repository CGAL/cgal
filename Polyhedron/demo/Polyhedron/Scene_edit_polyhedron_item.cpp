#include "Scene_edit_polyhedron_item.h"
#include "Kernel_type.h"
#include "Polyhedron_type.h"

#include <QVariant>
#include <list>

#include <QObject>
#include <QMenu>
#include <QAction>
#include <CGAL/gl_render.h>

struct Scene_edit_polyhedron_item_priv {
  Scene_polyhedron_item* poly_item;
};

Scene_edit_polyhedron_item::Scene_edit_polyhedron_item(Scene_polyhedron_item* poly_item)
  : d(new Scene_edit_polyhedron_item_priv)
{
  d->poly_item = poly_item;
}

Scene_edit_polyhedron_item::~Scene_edit_polyhedron_item()
{
  delete d;
}

Scene_edit_polyhedron_item* 
Scene_edit_polyhedron_item::clone() const {
  return 0;
}

QString 
Scene_edit_polyhedron_item::toolTip() const
{
  if(!d->poly_item->polyhedron())
    return QString();

  return QObject::tr("<p>Polyhedron <b>%1</b> (mode: %5, color: %6)</p>"
                     "<p>Number of vertices: %2<br />"
                     "Number of edges: %3<br />"
                     "Number of facets: %4</p>")
    .arg(this->name())
    .arg(d->poly_item->polyhedron()->size_of_vertices())
    .arg(d->poly_item->polyhedron()->size_of_halfedges()/2)
    .arg(d->poly_item->polyhedron()->size_of_facets())
    .arg(this->renderingModeName())
    .arg(this->color().name());
}

void Scene_edit_polyhedron_item::draw() const {
  d->poly_item->draw();
}

Polyhedron* 
Scene_edit_polyhedron_item::polyhedron()       { return d->poly_item->polyhedron(); }
const Polyhedron* 
Scene_edit_polyhedron_item::polyhedron() const { return d->poly_item->polyhedron(); }

bool
Scene_edit_polyhedron_item::isEmpty() const {
  return d->poly_item->isEmpty();
}

Scene_edit_polyhedron_item::Bbox
Scene_edit_polyhedron_item::bbox() const {
  return d->poly_item->bbox();
}


void
Scene_edit_polyhedron_item::
changed()
{
  d->poly_item->changed();
  Scene_item::changed();
}

void 
Scene_edit_polyhedron_item::select(double orig_x,
                                   double orig_y,
                                   double orig_z,
                                   double dir_x,
                                   double dir_y,
                                   double dir_z)
{
  d->poly_item->select(orig_x,
                       orig_y,
                       orig_z,
                       dir_x,
                       dir_y,
                       dir_z);
}

Scene_polyhedron_item* Scene_edit_polyhedron_item::to_polyhedron_item() const {
  return d->poly_item;
}

#include "Scene_edit_polyhedron_item.moc"
