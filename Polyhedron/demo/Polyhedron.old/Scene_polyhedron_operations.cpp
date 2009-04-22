#include "Scene.h"
#include "Polyhedron_type.h"
#include <CGAL/IO/Polyhedron_iostream.h>

Polyhedron* Scene::new_polyhedron() 
{
  return new Polyhedron;
}

Polyhedron* Scene::copy_polyhedron(Polyhedron* poly)
{
  return new Polyhedron(*poly);
}

void Scene::destroy_polyhedron(Polyhedron* poly)
{
  delete poly;
}

bool Scene::load_polyhedron(Polyhedron* poly, std::istream& in)
{
  in >> *poly;
  return in;
}

bool Scene::save_polyhedron(Polyhedron* poly, std::ostream& out)
{
  out << *poly;
  return out;
}

QString Scene::polyhedronToolTip(int index) const
{
  Polyhedron* poly = polyhedron(index);
  if(!poly)
    return QString();

  return tr("<p><b>%1</b> (mode: %5, color: %6)</p>"
            "<p>Number of vertices: %2<br />"
            "Number of edges: %3<br />"
            "Number of facets: %4</p>")
    .arg(polyhedronName(index))
    .arg(poly->size_of_vertices())
    .arg(poly->size_of_halfedges()/2)
    .arg(poly->size_of_facets())
    .arg(polyhedra[index].rendering_mode == Wireframe ? "wireframe" : "fill")
    .arg(polyhedra[index].color.name());
}
