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

void Scene::destroy(Polyhedron* poly)
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

Scene::Bbox Scene::bbox()
{
  if(polyhedra.empty()) {
    Bbox bbox;
    bbox.xmin = bbox.ymin = bbox.zmin = 0.0;
    bbox.xmax = bbox.ymax = bbox.zmax = 1.0;
    return bbox;
  }
  else
  {
    Point p = polyhedra.begin()->polyhedron_ptr->vertices_begin()->point();
    CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    for(Polyhedra::iterator 
          poly_it = polyhedra.begin(),
          poly_end = polyhedra.end();
        poly_it != poly_end; ++poly_it) {
      for(Polyhedron::Vertex_iterator
            v = poly_it->polyhedron_ptr->vertices_begin(),
            v_end = poly_it->polyhedron_ptr->vertices_end();
          v != v_end; ++v)
      {
        bbox = bbox + v->point().bbox();
      }
    }
    Bbox result;
    result.xmin = bbox.xmin();
    result.ymin = bbox.ymin();
    result.zmin = bbox.zmin();
    result.xmax = bbox.xmax();
    result.ymax = bbox.ymax();
    result.zmax = bbox.zmax();
    return result;
  }
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
