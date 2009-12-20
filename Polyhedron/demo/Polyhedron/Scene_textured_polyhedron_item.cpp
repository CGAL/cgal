#include "Scene_textured_polyhedron_item.h"
#include "Textured_polyhedron_type.h"
#include <CGAL/IO/Polyhedron_iostream.h>

#include <QObject>
#include <CGAL/gl_render.h>

typedef EPIC_kernel::Point_3 Point;

Scene_textured_polyhedron_item::Scene_textured_polyhedron_item()
  : Scene_item_with_display_list(),
    poly(new Textured_polyhedron)
{
  texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
}

Scene_textured_polyhedron_item::Scene_textured_polyhedron_item(Textured_polyhedron* const p)
  : Scene_item_with_display_list(),
    poly(p)
{
  texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
}

Scene_textured_polyhedron_item::Scene_textured_polyhedron_item(const Textured_polyhedron& p)
  : Scene_item_with_display_list(),
    poly(new Textured_polyhedron(p))
{
  texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
}

// Scene_textured_polyhedron_item::Scene_textured_polyhedron_item(const Scene_textured_polyhedron_item& item)
//   : Scene_item_with_display_list(item),
//     poly(new Textured_polyhedron(*item.poly))
// {
// }

Scene_textured_polyhedron_item::~Scene_textured_polyhedron_item()
{
  delete poly;
}

Scene_textured_polyhedron_item* 
Scene_textured_polyhedron_item::clone() const {
  return new Scene_textured_polyhedron_item(*poly);
}

// Load textured_polyhedron from .OFF file
bool
Scene_textured_polyhedron_item::load(std::istream& in)
{
  in >> *poly;
  return in && !isEmpty();
}

// Write textured_polyhedron to .OFF file
bool 
Scene_textured_polyhedron_item::save(std::ostream& out) const
{
  out << *poly;
  return out;
}

QString 
Scene_textured_polyhedron_item::toolTip() const
{
  if(!poly)
    return QString();

  return QObject::tr("<p>Textured polyhedron <b>%1</b> (mode: %5, color: %6)</p>"
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

// Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
void Scene_textured_polyhedron_item::direct_draw() const {
  glTexImage2D(GL_TEXTURE_2D,
	       0,
	       GL_RGB,
	       texture.GetWidth(),
	       texture.GetHeight(),
	       0,
	       GL_RGB,
	       GL_UNSIGNED_BYTE, 
	       texture.GetData());
  glEnable(GL_TEXTURE_2D);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

  poly->gl_draw_textured_triangles(true, true, 1.0);
  glDisable(GL_TEXTURE_2D);
}

Textured_polyhedron* 
Scene_textured_polyhedron_item::textured_polyhedron()       { return poly; }
const Textured_polyhedron* 
Scene_textured_polyhedron_item::textured_polyhedron() const { return poly; }

bool
Scene_textured_polyhedron_item::isEmpty() const {
  return (poly == 0) || poly->empty();
}

Scene_textured_polyhedron_item::Bbox
Scene_textured_polyhedron_item::bbox() const {
  const Point& p = *(poly->points_begin());
  CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
  for(Textured_polyhedron::Point_iterator it = poly->points_begin();
      it != poly->points_end();
      ++it) {
    bbox = bbox + it->bbox();
  }
  return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
              bbox.xmax(),bbox.ymax(),bbox.zmax());
}

#include "Scene_textured_polyhedron_item.moc"
