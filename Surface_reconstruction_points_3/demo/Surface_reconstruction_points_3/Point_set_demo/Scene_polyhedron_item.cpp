#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#include <CGAL/IO/Polyhedron_iostream.h>

#include <QObject>
#include <CGAL/gl_render.h>

Scene_polyhedron_item::Scene_polyhedron_item()
  : Scene_item_with_display_list(),
    poly(new Polyhedron)
{
}

Scene_polyhedron_item::Scene_polyhedron_item(Polyhedron* const p)
  : Scene_item_with_display_list(),
    poly(p)
{
}

Scene_polyhedron_item::Scene_polyhedron_item(const Polyhedron& p)
  : Scene_item_with_display_list(),
    poly(new Polyhedron(p))
{
}

// Scene_polyhedron_item::Scene_polyhedron_item(const Scene_polyhedron_item& item)
//   : Scene_item_with_display_list(item),
//     poly(new Polyhedron(*item.poly))
// {
// }

Scene_polyhedron_item::~Scene_polyhedron_item()
{
  delete poly;
}

Scene_polyhedron_item*
Scene_polyhedron_item::clone() const {
  return new Scene_polyhedron_item(*poly);
}

// Loads polyhedron from .OFF file
bool
Scene_polyhedron_item::load(std::istream& in)
{
  in >> *poly;
  return in && !isEmpty();
}

// Write polyhedron to .OFF file
bool
Scene_polyhedron_item::save(std::ostream& out) const
{
  out << *poly;
  return out;
}

QString
Scene_polyhedron_item::toolTip() const
{
  if(!poly)
    return QString();

  return QObject::tr("<p><b>%1</b> (mode: %5, color: %6)</p>"
                     "<p>Number of vertices: %2<br />"
                     "Number of edges: %3<br />"
                     "Number of facets: %4</p>")
    .arg(this->name())
    .arg(poly->size_of_vertices())
    .arg(poly->size_of_halfedges()/2)
    .arg(poly->size_of_facets())
    .arg(this->renderingModeName())
    .arg(this->color().name());
}


bool Scene_polyhedron_item::supportsRenderingMode(RenderingMode m) const {
  return m != PointsPlusNormals && m != Splatting;
}

// Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
void Scene_polyhedron_item::direct_draw() const {
  gl_render_facets(*poly);
}

Polyhedron*
Scene_polyhedron_item::polyhedron()       { return poly; }
const Polyhedron*
Scene_polyhedron_item::polyhedron() const { return poly; }

bool
Scene_polyhedron_item::isEmpty() const {
  return (poly == 0) || poly->empty();
}

Scene_polyhedron_item::Bbox
Scene_polyhedron_item::bbox() const {
  const Point& p = *(poly->points_begin());
  CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
  for(Polyhedron::Point_iterator it = poly->points_begin();
      it != poly->points_end();
      ++it) {
    bbox = bbox + it->bbox();
  }
  return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
              bbox.xmax(),bbox.ymax(),bbox.zmax());
}

#include "Scene_polyhedron_item.moc"
