#include "Scene.h"
#include "Textured_polyhedron_type.h"

Textured_polyhedron* Scene::new_tex_polyhedron() 
{
  return new Textured_polyhedron;
}

Textured_polyhedron* Scene::copy_tex_polyhedron(Textured_polyhedron* poly)
{
  return new Textured_polyhedron(*poly);
}

void Scene::destroy_tex_polyhedron(Textured_polyhedron* poly)
{
  delete poly;
}

QString Scene::texPolyhedronToolTip(int index) const
{
  Textured_polyhedron* poly = texPolyhedron(index);
  if(!poly)
    return QString();

  return tr("<p><b>%1</b> (mode: %5, color: %6)</p>"
            "<p>Number of vertices: %2<br />"
            "Number of edges: %3<br />"
            "Number of facets: %4</p>")
    .arg("textured polyhedron")
    .arg(poly->size_of_vertices())
    .arg(poly->size_of_halfedges()/2)
    .arg(poly->size_of_facets())
    .arg(polyhedra[index].rendering_mode == Wireframe ? "wireframe" : "fill")
    .arg(polyhedra[index].color.name());
}


