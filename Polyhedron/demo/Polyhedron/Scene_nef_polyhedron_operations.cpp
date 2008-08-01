#include "Scene.h"
#include "Nef_type.h"

Nef_polyhedron* Scene::new_nef_polyhedron() 
{
  return new Nef_polyhedron;
}

Nef_polyhedron* Scene::copy_nef_polyhedron(Nef_polyhedron* poly)
{
  return new Nef_polyhedron(*poly);
}

void Scene::destroy_nef_polyhedron(Nef_polyhedron* poly)
{
  delete poly;
}

QString Scene::nefPolyhedronToolTip(int index) const
{
  Nef_polyhedron* poly = nefPolyhedron(index);
  if(!poly)
    return QString();

  return tr("<p><b>%1</b> (mode: %5, color: %6)<br />"
	    "<i>Nef_3 polyhedron</i></p>"
            "<p>Number of vertices: %2<br />"
            "Number of edges: %3<br />"
            "Number of facets: %4<br />"
	    "number of volumes: %7</p>")
    .arg(polyhedronName(index))
    .arg(poly->number_of_vertices())
    .arg(poly->number_of_edges())
    .arg(poly->number_of_facets())
    .arg(polyhedra[index].rendering_mode == Wireframe ? "wireframe" : "fill")
    .arg(polyhedra[index].color.name())
    .arg(poly->number_of_volumes());
}
